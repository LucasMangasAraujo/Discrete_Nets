"""
Script defining the NetworkClass.
"""

import numpy as np
import networkx as nx
from collections import defaultdict
from scipy.interpolate import UnivariateSpline
import random

class NetworkClass:
    """
    A class to assist querying information from discrete networks. 
    This is the parent class, defining methods that are general for
    all types of networks
    """
    
    
    def __init__(self, data_file, dump_file, input_file):
        """
        Class constructor
        
        Inputs:
            data_file (str): name of LAMMPS data file
            dump_file (str): name of LAMMPS dump file
            input_file (str): name of LAMMPS inpute file
            
        """
        self.data_file = data_file
        self.dump_file = dump_file
        self.input_file = input_file
    
    
    
    def detect_tooLong_chains(self, model, cut_off = 0.95):
        """
        Detect ids of initially too long chains.
        
        Inputs:
            model (string):type of bond behaviour.
                    model = '1': Gaussian
                    model = '2': FJC
                    model = '3': Breakable extensible FJC
                    model = '4': Breakable FJC
            cut_off (float, default: 0.95): fraction of Nb defining a tool 
                                                elongated chain.
            
        Outputs:
            tooLong_ids: set with the ids of the chains that are above 
                            the cutoff
        """
        # Get equilrated node positions
        Nodes, _  = self.get_nodes_and_bonds()
        
        # Get Bonds and bond coefficients
        Bonds, Coeffs = self.get_bonds_and_coeffs(model)
        
        # Scan network to detect the broken chains
        tooLong_ids = set()
        for idx, (bond_type, n1, n2) in Bonds.items():
            ## Compute distance
            dist = np.linalg.norm(Nodes[n1] - Nodes[n2])
            
            ## Asses if failure occured
            chain_coeffs = Coeffs[bond_type]
            if NetworkClass.is_too_long(dist, model, chain_coeffs, cut_off):
                tooLong_ids.add(idx)
        
        
        return tooLong_ids
        
    
    
    @staticmethod
    def is_too_long(r, model, chain_coeffs, cut_off = 0.95):
        """
        Detect if chain is broken
        
        Inputs:
            r(float): chain end-to-end distance
            
            model (string):type of bond behaviour.
                    model = '1': Gaussian
                    model = '2': FJC
                    model = '3': Breakable extensible FJC
                    model = '4': Breakable FJC
                    
            chain_coeffs (tuple): chain parameters
            
        Output
            True if chain is too long.
            False if chain is not too long.
        """
        bKuhn, NKuhn = chain_coeffs[:2]
        Nb = bKuhn * NKuhn
        
        return (r / Nb) >= cut_off
    
    
    
    def rewrite_data_file(self, data_file, tooLong_ids, model):
        """
        Rewrite the data file removing bonds that are too long initially
        
        Inputs:
            data_file (string): LAMMPS data file.
            tooLong_ids (set): ids of chains initially too long
        """
        
        
        # Get equilrated node positions
        Nodes, _  = self.get_nodes_and_bonds()
        
        # Get Bonds and bond coefficients
        Bonds, Coeffs = self.get_bonds_and_coeffs(model)
        bond_ids = set(Bonds.keys())
        nTypes = len(tuple(Coeffs.keys()))
        
        # Generate set with the ids of the intact chains
        intact_ids = bond_ids - tooLong_ids
        intact_Bonds = {}
        if nTypes > 1: intact_Coeffs = {}
        
        temp = 0
        for idx in intact_ids:
            temp +=1
            
            ## Check if there is any source of polydipersity
            if nTypes > 1:
                ## Store previous bond type
                bond_type = Bonds[idx][0]
                intact_Coeffs[temp] = Coeffs[bond_type]
                intact_Bonds[temp] = temp, Bonds[idx][1], Bonds[idx][2]
            else:
                intact_Bonds[temp] = Bonds[idx]
        
        # Store lines from the original data file
        with open(data_file, "r") as f:
            lines = f.readlines()
            
        
        # Based on the new dict of bonds, rewrite data
        with open(data_file, "w") as f:
            ## Reuse some information 
            for line in lines:
                if 'atoms' in line:
                    f.write('%d atoms\n' %len(Nodes))
                elif 'bonds' in line:
                    f.write('%d bonds\n' %len(intact_Bonds));
                elif nTypes > 1 and ('bond types' in line):
                    f.write(f'{len(intact_Coeffs)} bond types\n')
                elif nTypes > 1 and ('Bond Coeffs' in line):
                    ## Write header of the coefficients and skip one line
                    f.write(line)
                    f.write('\n')
                    
                    ## Write the coefficients of each bond
                    for idx, coeffs in intact_Coeffs.items():
                        f.write(f'{idx} {coeffs[0]} {coeffs[1]}\n')
                    f.write('\n')
                    
                    ## Write the header of the atoms section and break
                    f.write("Atoms\n")
                    break
                    
                else:
                    f.write(line);
                if 'Atoms' in line:
                    break;
                    
            
            ## Write the nodes
            f.write('\n')
            for idx in Nodes.keys():
                f.write('%d 1 1 %g %g %g\n' %(idx,Nodes[idx][0],Nodes[idx][1],Nodes[idx][2]))
            
            ## Velocities
            f.write('\n')
            # Bonds IDs and pairs of nodes connected by each bond
            f.write('Velocities\n\n') 
            for idx in Nodes:
                f.write('%d 0 0 0\n' %idx)
            
            ## Now the bonds
            f.write('\n')
            ## Bonds IDs and pairs of nodes connected by each bond
            f.write('Bonds\n\n') 
            for idx in intact_Bonds.keys():
                f.write('%d %d %d %d\n' %(idx, intact_Bonds[idx][0], intact_Bonds[idx][1], intact_Bonds[idx][2]));
        
        
        
        
        return
    
    def compute_macro_stretch(self, loading, ref_lengths):
        """
        Obtain stretch based on the box dimension and load type
        
        Inputs:
            loading (int): type of loading. See runinc in sim_executor 
                           for details.
           ref_lengths (dict): lengths of the reference configuration.
       
       Outputs:
            macro_stretch (float): stretch charactherising the load.
        """
        
        # Get the box dimension 
        box_boundaries, box_lengths = self.get_box()
        
        if loading < 3:
            Lx0 = ref_lengths['x']
            Lx = box_lengths['x']
            stretch = Lx / Lx0
        
        
        return macro_stretch
    
    def get_avg_preStretch(self, model):
        """
        Return the rms pre-stretch
        
        Inputs:
            model (str): type of chain model used. See sim_exectutor.py
                         for details.
        
        Outputs:
            lambda0 (float): root-mean-square pre-stretch
        """
        
        # Get the pre-stretch distributio
        preStretch_distr = self.get_preStretch_distr(model)
        
        # Take the average 
        avg_squared_preStretch = np.mean(np.square(preStretch_distr))
        
        return np.sqrt(avg_squared_preStretch)
    
    def get_rms_r0(self):
        """
        Inputs:
            None
            
        Outputs:
            root mean square end-to-end distance
        """
        # Calculate distances in the network
        distances = np.array(tuple(self.get_distances().values()))
        
        # Obtain root mean square values
        rms_r0 = np.sqrt(np.mean(distances **2))
        
        return rms_r0
    
    
    
    def get_preStretch_distr(self, model):
        """
        Get the pre-stretch distribution.
        """
        # Extract Nodes, Bonds, and BondsCoeffs
        Nodes, _ = self.get_nodes_and_bonds()
        Bonds, BondCoeffs = self.get_bonds_and_coeffs(model)
        
        
        # Loop over the Bonds dict
        preStretch_distr = []
        for idx, (bond_type, n1, n2) in Bonds.items():
            v = Nodes[n1] - Nodes[n2]
            r0 = np.linalg.norm(v)
            if int(model) in [2, 4] :
                b, N = BondCoeffs[bond_type][0], BondCoeffs[bond_type][1]
                lambda0 = r0 / (np.sqrt(N)* b)
            elif model == '6':
                kappa, _, _, _ = BondCoeffs[bond_type]
                lambda0 = r0 * np.sqrt(2 * kappa / 3)
                
            
            
            preStretch_distr.append(lambda0)
        
        
        return np.array(preStretch_distr)
    
    
    
    
    def get_computational_params(self, params, dim = 3):
        """
        Get computational params used in the simulation.
        Inputs:
            params (tuple): parameters in with physical units when relevant. 
                            The order is the following:
                                bKuhn: Kuhn length in nm
                                NKuhn: Number of Kuhn segments in the chain.
                                nub3: Normalised chain density in bKuhn3 units.
            dim (int, default = 3): problem dimention.
        
        Outputs:
            computational_params (tuple): parameters in computational units.
        """
        # Extract information DN structure
        Nodes, Bonds = self.get_nodes_and_bonds()
        Boundary = self.get_boundary()
        
        # Unpack input params 
        bKuhn, NKuhn, nub3 = params
        
        # Normalise using the density of crosslinks (subtracting bounary nodes)
        crosslinks = len(Nodes) - len(Boundary)
        upsilonb3 = nub3 / 2
        computational_bKuhn = np.power( upsilonb3 / crosslinks, 1/dim)
        
        # Assemple computational_params tuple
        computational_params = (computational_bKuhn, NKuhn)
        
        return computational_params
    
    
    def create_DN_graph(self):
        """
        Turn DN into Graph.
        
        Inputs:
            None
            
        Outputs:
            G (networkx Graph): networkx graph object.
        """
        # Get Network structure
        Nodes, Bonds = self.get_nodes_and_bonds()
        
        # Create Graph
        G = nx.Graph();
        G.add_nodes_from(Nodes);
        for idx, (n1,n2) in Bonds.items():
            dist = np.linalg.norm(Nodes[n1] - Nodes[n2])
            G.add_edge(n1,n2, weigth = dist);
        
        return G

        
        
    def get_boundary(self, initial_Nodes = None):
        """
        Get bounday nodes ids
        Inputs:
            initial_Nodes (dict, default = None): Coords dict in the ref config 
                                                    in format (x, y, z)
        
        Outputs:
            Boudary (tuple): ids of boundary nodes
        """
        
        
        # Read main file to find idx of the boundary nodes
        try:
            with open(self.input_file, "r") as f:
                ## Read file until group of boundary nodes is found
                key = f.readline()
                while 'group' not in key:
                    key = f.readline()
                
                ## Loop over split line and store ids of boundary nodes
                data = key.strip('\n').split(" ")
                Boundary = []
                for idx in data:
                    try:
                        Boundary.append(int(idx))
                    except ValueError:
                        continue
                
        except FileNotFoundError:
            print("Warning: input file not found. Using coordinates in initial_Nodes instead.")
            Boundary = []
            ref_G = self.create_DN_graph()
            for idx, coord in initial_Nodes.items():
                is_boundary = np.any(np.isclose(coord, 0)) or np.any(np.isclose(coord, 1))
                if is_boundary:
                    ## Check if might be an internal node
                    if ref_G.degree[idx] == 1:
                        Boundary.append(idx)
                
        
        return tuple(Boundary)
    
    
    def get_distances(self):
        """
        Get end-to-end distance of the chains at a given configuration.
        
        Inputs: 
            None
            
        Outputs:
            distances (dict): distances of regular chains.
        """
        
        # Get Node coordinates and their positions
        Nodes, Bonds = self.get_nodes_and_bonds()
        
        # Scan and store results
        distances = {}
        for idx, bond in Bonds.items():
            n1, n2 = bond
            vector = Nodes[n1] - Nodes[n2]
            distances[idx] = np.linalg.norm(vector)
        
        
        return distances


    def get_forces(self, model):
        """
        Get chain forces
        
        Inputs: 
            model (str): chain model used
            
        Outputs:
            fbkT (ndarray): array of normalised chain forces
            chain_forces (dict): same information but with in dict form.
        """
        # Get coefficients and bonds
        Bonds, Coeffs = self.get_bonds_and_coeffs(model)
        Nodes, _ = self.get_nodes_and_bonds()
        
        # Loop
        fbkT, bond_forces = [], {}
        for idx, (bond_type, n1, n2) in Bonds.items():
            dist = np.linalg.norm(Nodes[n1] - Nodes[n2])
            if model == "2":
                bKuhn, NKuhn = Coeffs[bond_type]
                rNb = dist / (NKuhn * bKuhn)
                fbkT.append(NetworkClass.invLangevin(rNb))
            elif model == "4":
                bKuhn, NKuhn, critical_rNb = Coeffs[bond_type]
                rNb = dist / (NKuhn * bKuhn)
                fbkT.append(NetworkClass.invLangevin(rNb))
                
            chain_forces[idx] = NetworkClass.invLangevin(rNb)
        
        
        return np.array(fbkT), chain_forces



    def get_bonds_and_coeffs(self, model):
        """
        Get bonds and their coefficients
        
        Inputs:
            model (str): chain model used
            
        Outputs:
            Bonds (dict): bonds and bond type in format (node1, node2, bond_type)
            BondCoeffs(dict): coefficients of each bond in formater (bond_ids, coeffs)
        """
        
        # Read data file
        with open(self.data_file,"r") as f:
            ## Get the Bond types and their coefficients
            key = f.readline()
            while "Bond Coeffs" not in key:
                key = f.readline()
            
            BondCoeffs = {}
            f.readline()
            data = f.readline().strip("\n").split(" ")
            while len(data) > 1:
                if model == '2':
                    BondCoeffs[int(data[0])] = float(data[1]), float(data[2])
                elif model == '4':
                    BondCoeffs[int(data[0])] = float(data[1]), float(data[2]), float(data[3])
                elif model == '6':
                    BondCoeffs[int(data[0])] = float(data[1]), float(data[2]), float(data[3]), float(data[4])
                
                data = f.readline().strip("\n").split(" ")
            
            
            ## Keep readin until bond section is reached
            key = f.readline()
            while "Bonds" not in key:
                key = f.readline()
            f.readline()
            
            ## Read bonds
            Bonds = {}
            data = f.readline().strip("\n").split(" ")
            while len(data) > 1:
                Bonds[int(data[0])] = int(data[1]), int(data[2]), int(data[3])
                data = f.readline().strip("\n").split(" ")
        
        
        return Bonds, BondCoeffs



    def get_nodes_and_bonds(self):
        """
        Get the bonds and node coordinates of the network.
        Inputs:
            None
            
        Outputs:
            Nodes (dict): dictionary containing node positions
            Bonds (dict): dictionary containing bonds
        """
        # Get data file containing network 
        filename = self.data_file
        
        # Initialize dict
        Bonds = {}
        Nodes = {}
        
        with open(filename, 'r') as f:
            key = f.readline().strip()

            while 'Atoms' not in key:
                key = f.readline().strip()

            f.readline()  # Skip header

            data = f.readline().split()
            while len(data) > 1:
                idx = int(data[0])
                x = float(data[3])
                y = float(data[4])
                z = float(data[5])
                Nodes[idx] = np.array([x, y, z]) 

                data = f.readline().split()
            
            key = f.readline().strip()

            while 'Bonds' not in key:
                key = f.readline().strip()

            f.readline()  # Skip header
            
            data = f.readline().split()
            while len(data) > 1:
                idx = int(data[0])
                n1 = int(data[2])
                n2 = int(data[3])
                Bonds[idx] = [n1, n2] 

                data = f.readline().split()
        
        return Nodes, Bonds
    
    
    def get_box(self):
        """
        Get current box bounds.
        
        Inputs:
            
        Outputs:
            box_boundaries (dict): coords defininf the planes of each boundary.
            box_lengths (dict): dimensions of the simulation box.
        """
        # Initialise dicts
        box_boundaries = {}
        box_lengths = {}
        
        # Read file
        with open(self.data_file, "r") as f:
            key = f.readline()
            ## kepp reading until relevant section is reached
            while "xlo" not in key:
                key = f.readline()
            ## Read x length
            data = key.strip("\n").split(" ")
            box_lengths['x'] = float(data[1]) - float(data[0])
            box_boundaries['x'] = float(data[0]), float(data[1])
            
            ## Read y length
            data = f.readline().strip("\n").split(" ")
            box_lengths['y'] = float(data[1]) - float(data[0])
            box_boundaries['y'] = float(data[0]), float(data[1])
            
            ## Finally z length
            data = f.readline().strip("\n").split(" ")
            box_lengths['z'] = float(data[1]) - float(data[0])
            box_boundaries['z'] = float(data[0]), float(data[1])
            
        
        
        return box_boundaries, box_lengths
    
    
    
    def calculate_stress(self, dim):
        """
        Calculate  cauchy stress using virtual work principle.
        
        Inputs:
            dim (int): problem dimension
        
        Outputs:
            S (ndarray): 3x1 array with the principal stresses
        """
        S = np.zeros(3)
        
        with open(self.dump_file, 'r') as f:
            key = f.readline().strip()
            
            while 'id' not in key:
                key = f.readline().strip()
            
            data = f.readline().split()
            while len(data) > 1:
                if dim == 3:
                    x = float(data[2])
                    y = float(data[3])
                    z = float(data[4])
                    fx = -float(data[5])  # minus sign for reaction force
                    fy = -float(data[6])
                    fz = -float(data[7])

                    S[0] += fx * x
                    S[1] += fy * y
                    S[2] += fz * z

                else:
                    x = float(data[2])
                    y = float(data[3])
                    fx = -float(data[4])  
                    fy = -float(data[5])

                    S[0] += fx * x
                    S[1] += fy * y
                    
                data = f.readline().split()

        return S
    
    
    @staticmethod
    def calculate_nominal_stress(dim, F, cauchy_stress):
        """
        Calculate nominal stress stress using the virtual work principle.
        
        Inputs:
            dim (int): problem dimension
            F (ndarray): deformation gradient
        
        Outputs:
            S (ndarray): 3x1 array with the principal stresses
        """
        # Assemble F^{-T} in the principal space
        F_minusT = np.array([1/stretch for stretch in F])
        
        # Use regular relation to obtain the nominal stress
        nominal = np.zeros(3)
        nominal = cauchy_stress * F_minusT
        
        return nominal
    
    
    
    @staticmethod
    def stress_units(stress_array, bKuhn, T = 298):
        """
        Go from b3/kT units to kPa
        
        Inputs:
            stress_array (ndarray): rubbery components of stress in b^3/kT units
            bKuhn (float): Kuhn length in nm.
            T (float, default = 298 K): temperature in Kelvin.
            
        Ouputs:
            stress_kPa (ndarray): stress array in kPa.
            
        """
        # Declare Boltzmann constant
        kB = 1.380649e-23
        kT = kB * T ## temperatur in energy units
        
        # Render stress array with J/nm3 units
        stress_J_over_nm3 = stress_array * kT / np.power(bKuhn, 3)
        
        # Convert to kPa
        stress_kPa = stress_J_over_nm3 * 1e24
        return stress_kPa
    
    
    @staticmethod
    def invLangevin(x):
        '''
        '''
        
        #Pade approximation of the inverse Langevin function
        #cf. Cohen 1991

        invL = x*(3.-(x**2))/(1.-(x**2))

        return invL
    

class FracNetworkClass(NetworkClass):
    """
    A class for networks where detersministic scissions are allowed.
    Inherented from the NetworkClass.
    
    """
    
    
    def get_chain_stretches(self, initial_Nodes, F):
        """
        Get chain streches for the surviving chains
        """
        Nodes, Bonds = self.get_nodes_and_bonds()
        
        ## Loop
        chain_stretch_array = []
        affine_chain_stretch_array = []
        for idx, (n1, n2) in Bonds.items():
            v = Nodes[n1] - Nodes[n2]
            v0 = initial_Nodes[n1] - initial_Nodes[n2]
            chain_stretch = np.linalg.norm(v) / np.linalg.norm(v0)
            chain_stretch_array.append(chain_stretch)
            
            ## Compute affine stretch now
            affine_chain_stretch = np.linalg.norm(F*v0) / np.linalg.norm(v0)
            affine_chain_stretch_array.append(affine_chain_stretch)
            
        
        return np.array(chain_stretch_array), np.array(affine_chain_stretch_array)
    
    
    
    
    
    def detect_broken_chains(self, model, extra_params):
        """
        Detect ids of broken chains once after mechanical 
        equilibrium was reached.
        """
        # Get equilrated node positions
        Nodes, _  = self.get_nodes_and_bonds()
        
        # Get Bonds and bond coefficients
        Bonds, Coeffs = self.get_bonds_and_coeffs(model)
        
        # Scan network to detect the broken chains
        broken_ids = set()
        for idx, (bond_type, n1, n2) in Bonds.items():
            ## Compute distance
            dist = np.linalg.norm(Nodes[n1] - Nodes[n2])
            
            ## Asses if failure occured
            chain_coeffs = Coeffs[bond_type]
            if FracNetworkClass.is_broken(dist, model, chain_coeffs, extra_params):
                broken_ids.add(idx)
        
        
        return broken_ids
        
    
    
    
    def remove_ineffective_clusters(self, G):
        """
        Remover clusters that are coiled from the Graph.
        
        Inputs:
            G (networkx graph obj): self-explanatory
            
        Outputs: 
            None
        """
        # Get current positions of the nodes
        Nodes, _ = self.get_nodes_and_bonds()
        
        # Get connected components of the graph
        connected_components = list(nx.connected_components(G))
        
        # Removed the connected components that are coiled
        ineffective_clusters_ids = []
        for i, component in enumerate(connected_components):
            subG = G.subgraph(component)
            sub_edges = list(subG.edges)
            subG_distances = []
            
            for bond in sub_edges:
                n1, n2 = bond
                dist = np.linalg.norm(Nodes[n1] - Nodes[n2])
                subG_distances.append(dist)
            
            ## check if all distances are close to zero
            is_coiled = np.all(np.array(subG_distances) < 1e-6)
            if is_coiled:
                ineffective_clusters_ids.append(i)
        
        # With ids of the ineffective clusters, remove them from the graph
        for i in ineffective_clusters_ids:
            G.remove_nodes_from(connected_components[i])
        
        return