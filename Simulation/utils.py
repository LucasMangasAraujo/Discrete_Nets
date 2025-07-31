"""
    This script contains functions that are called by the control code
    during discrete network simulations.
"""
from math import *
import numpy as np
import sys, os

# Import NetworkClass (new_feature)
from network_class import NetworkClass


def writeMain(simfile,posfile,Boundary,dim,model):
    """ 
    Writes the main input file in LAMMPS.
    
    simfile : name of the input files in LAMMPS.
    posfile : name of the file for postprocessing the results
            from the energy minimisation.
    Boundary : list containing the IDs of the nodes at the boundary 
            of the RVE.
    dim : integer indicating the dimension of the problem.
    model : string indicating the type of bond behaviour.
            model = '1': Gaussian
            model = '2': FJC
            model = '3': Breakable extensible FJC
            model = '4': Breakable FJC
    
    The function returns None.
    """
    
    min_algo='fire' #algorithm for minimization
    dmax = 0.05      #how much a single atom can move during line search

    if model == '1':
        bond_style = 'harmonic'
        
    elif model == '2' or model == '4':
        bond_style = 'langevin'
        
    elif model == '3':
        bond_style = 'Xlangevin'
        
        
    else:
        print('unknown bond type: %s' %model)
        exit()

    f = open(simfile,'w')


    f.write('#Main input file for LAMMPS\n')

    f.write('units\tlj\n')
    f.write('dimension\t%d\n' %dim)
    if dim == 3:
        f.write('boundary\tf f f \n')
    else:
        f.write('boundary\tf f p\n')
    
    f.write('atom_style\tbond\n')
    f.write('bond_style\t%s\n' %(bond_style))
    f.write('atom_modify\tsort 0 0\n')
    f.write('pair_style\tnone\n\n')

    f.write('read_data\t%s\n\n' %posfile)

    f.write('reset_timestep\t0\n')
    f.write('timestep\t0.0001\n')
    f.write('neighbor\t0.1 nsq\n')
    f.write('thermo\t1\n')
    if dim == 3:
        f.write('thermo_style\tcustom etotal press pxx pyy pzz pxy pxz pyz\n')
    else:
        f.write('thermo_style\tcustom etotal press pxx pyy pxy\n')
    
    f.write('min_style\t%s\n' %(min_algo))
    f.write('min_modify\tdmax %s\n\n' %(dmax))

    f.write('group\tboundary id ')
    for i in range(len(Boundary)):
        f.write('%s ' %(Boundary[i]))
    f.write('\n\n')

    # Step 1: deform the box affinely 
    #delta values: change in box boundaries at the end of run  
    #Note: actual mode of deformation applied here is not important as these lines will be replaced
    #by run.py on the go
    if dim == 3:
        f.write('fix 1 all deform 1 x delta 0 0 y volume z volume remap x units box\n')
    else:
        f.write('fix 1 all deform 1 x delta 0 0 y volume remap x units box\n')

    #need a run to apply the fix deform command above
    f.write('run 1\n\n')

    # Step 2: Apply zero force on boundary nodes (prevent their motion) and minimize energy
    f.write('fix\t2 boundary setforce 0 0 0\n')
    f.write('minimize\t0 1e-16 1000 10000\n\n')

    # Step 3: remove the zero-force constraint on the boundary
    f.write('unfix 2\n\n')

    # Define computation to calculate forces
    if dim == 3:
        f.write('compute\t1 boundary property/atom fx fy fz\n')
        f.write('dump\t1 boundary custom 1 test.res id type x y z c_1[1] c_1[2] c_1[3]\n')

    else:
        f.write('compute\t1 boundary property/atom fx fy\n')
        f.write('dump\t1 boundary custom 1 test.res id type x y c_1[1] c_1[2]\n')
        
    f.write('dump_modify\t1 sort id\n')

    #run dummy step (0 increment) to perform the dump operation and write test.res
    f.write('run\t0\n\n')   
    
    #write new atom positions
    f.write('write_data\t%s\n\n' %posfile)

    f.close()
    

    return


def writePositions(filename, Nodes, Bonds, Boundary, BondTypes, model, params):

    """
    Writes a file containing the architecture of the discrete network 
    and the parameters of each chain in it in a way that LAMMPS can 
    read it and proceed with the energy minimisation process.
    
    
    filename : name of the file that will be generated
    Nodes : dictionary whose keys are the IDs of the nodes,
            and the values are a list with node coordinates
    Bonds: dictionary whose keys are the bond IDs and the,
            values are a list containing the pair of nodes 
            connected.
    Boundary: list of strings with the IDs of the boundary nodes
    BondTypes: dictionary whose keys are the bond IDs and the,
            values are the chain lengths.
    model : string indicating the type of bond behaviour.
            model = '1': Gaussian
            model = '2': FJC
            model = '3': Extensible FJC
    params: list containing chain parameters other than 
            the chain length.
            
    
    The function returns None.
    
    """
    
    # Calculate the number of nodes, bonds, and boundary nodes
    Natoms = len(Nodes);
    Nbonds = len(Bonds);
    Nboundary = len(Boundary);
    NbondTypes = len(BondTypes);
    
    # Check for polydispersity
    chain_lengths = np.array(list(BondTypes.values()));
    polydispersity_flag = not np.all(chain_lengths == chain_lengths[0])
    if not polydispersity_flag: NbondTypes = 1;
    
    # Open file and write on it
    with open(filename, 'w') as f:
        
        #Header
        f.write('LAMMPS data file for the initial network geometry\n\n');

        #Number of nodes and bonds
        f.write('%d atoms\n' %Natoms);
        f.write('1 atom types\n');
        f.write('%d bonds\n' %Nbonds);
        f.write('%d bond types\n\n' %NbondTypes);

        #Box dimensions
        f.write('-0.1 1.1 xlo xhi\n');
        f.write('-0.1 1.1 ylo yhi\n');
        f.write('-0.1 1.1 zlo zhi\n\n');

        #Masses
        f.write('Masses\n\n1 1\n\n');

        #Bond coefficients
        f.write('Bond Coeffs\n\n');
        bKuhn = params[0]; ## Kuhn length
        
        if polydispersity_flag:
            for idx, N in BondTypes.items():
                
                if model == '1': ## Gaussian chain (harmonic)
                    kappa = (3./2.) * (1 / ( N * pow(bKuhn, 2) )); ## Bond stiffness in the Gaussian regime
                    f.write('%d %g %g\n'%(idx, kappa, 0.)); ## zero rest length
                    
                elif model == '2': ## FJC with Non-Gaussian statistics
                    f.write('%d %g %g\n' %(idx, bKuhn, N));
                
                elif model == '3': ## Extensible FJC
                    bKuhn, Eb, critical_eng = tuple(params);
                    f.write('%d %g %g %g %g\n' %(idx, bKuhn, N, Eb, critical_eng));
                
            
        else:
            N = chain_lengths[0];
            
            if model == '1': ## Gaussian chain (harmonic)
                kappa = (3./2.) * (1 / ( N * pow(bKuhn, 2) )); ## Bond stiffness in the Gaussian regime
                f.write('1 %g %g\n'%(kappa, 0.)); ## zero rest length
                
            elif model == '2': ## Kuhn length and Kuhn segments 
                f.write('1 %g %g\n' %(bKuhn, N));
            
            elif model == '3': ## Extensible FJC
                bKuhn, Eb, critical_eng = tuple(params);
                f.write('1 %g %g %g %g\n' %(bKuhn, N, Eb, critical_eng));
            
            
        f.write('\n\n');
        
        #Atoms ids and positions
        f.write('Atoms\n\n')
        for idx in Nodes:
            f.write('%d 1 1 %g %g %g\n' %(idx,Nodes[idx][0],Nodes[idx][1],Nodes[idx][2]))
        
        f.write('\n')
        # Bonds IDs and pairs of nodes connected by each bond
        f.write('Bonds\n\n') 

        # Check for polydispersity
        for idx in Bonds:
            if(polydispersity_flag): ## Each bond has its own type
                f.write('%d %d %g %g\n' % (idx,idx,Bonds[idx][0],Bonds[idx][1]) );
            else: ## there is one bond type only
                f.write('%d 1 %g %g\n' %(idx,Bonds[idx][0],Bonds[idx][1]));

        f.write('\n')

    return


def readBondTypes(input,Nbonds):
    '''
    Reads the chain lengths and store them in dictionary called BondTypes
    
    input: pointer of the file being read.
    Nbonds: number of bonds in the network.
    
    The function returns
    1) A dictionary whose keys are the bond IDs and the,
        values are the chain lengths: BondTypes.
    2) A string with previous line read: key.
    '''
    
    # Initialise the dict
    BondTypes = {};
    
    # Read and store
    for i in range(1,Nbonds + 1):
        key = input.readline().strip('\n');
        data = key.split();
        BondTypes[i] = float(data[1]);
    
    return BondTypes, key


def readBoundary(input):

    """
    Read the list of boundary nodes from a file with pointer input
    and store them in an arra.
    
    input: pointer of the file being read.
    
    The function returns:
    1) A list of strings representing the IDs of boundary nodes:
        data.
    2) A string containing the last read line of the geometry file
        key.
    
    """
    key = input.readline().strip(' \n') #read the next line in the input file
    if ("," in key):
        data = key.split(',')
    else:
        data = key.split(' ');

    return data,key


def readBonds(input):

    """
    Read the list of bonds from a file with pointer input
    and store them in a dictionary.
    
    
    input: pointer of the file being read.
    
    The function returns:
    1) A dictionary containing the list of bonds: bond_dict.
    2) A string containing the last read line of the geometry file:
        key.
    
    """

    again = True;
    bond_dict = {};

    while again is True:
        key = input.readline().strip('\n')  #read the next line in the input file
        if '$' in key:
            again = False;
        else:
            data = key.split(',');
            bond_dict[int(data[0])] = [int(data[1]),int(data[2])]

    return bond_dict,key


def readNodes(input):

    """
    Read the list of node from a file with pointer input
    and store them in a dictionary.
    
    input: pointer of the file being read.
    
    The function returns:
    1) A dictionary containing the list of nodes coordinates:
        node_dict.
    2) A string containing the last read line of the geometry file:
        key.
    """

    #read nodes
    again = True
    node_dict = {}

    while again is True:

        #read the next line in the input file
        key = input.readline().strip('\n')

        if '$' in key:
            again = False
        else:
            data = key.split(',')
            x = float(data[1])
            y = float(data[2])
            z = float(data[3])
            node_dict[int(data[0])] = [float(x),float(y),float(z)]

    return node_dict,key


def readGeometry(filename):

    """ 
    Read network geometry file.
    
    geomfile : file or full path to file containing the 
                network architecture.
    
    The function returns:
    1) A dictionary containing lists of nodes coordinates:
        Nodes.
    2) A dictionary containing lists of bonds: Bonds.
    3) A list of strings representing the boundary
        nodes IDs: Boundary.
    4) A dictionary containing the chain lengths of 
        each bond: BondTypes.
    """

    f = open(filename,'r')

    again = True
    key = f.readline().strip('\n')  #read the next line in the input file

    while again is True:

        if 'nodes' in key:
            Nodes,key = readNodes(f)

        elif 'bonds' in key:
            Bonds, key = readBonds(f)

        elif 'boundary' in key:
            Boundary, key = readBoundary(f)
            key = f.readline().strip('\n') # Start Reading Again $
        elif 'BondTypes' in key:
            BondTypes, key = readBondTypes(f,len(Bonds));
        else:
            again = False

    f.close()

    return Nodes,Bonds,Boundary,BondTypes


def readParams(filename):

    """
    Read the simulation parameters file 'filename'.
    
    filename: string representing the input file for the
                simulation.
    
    The function returns
    1) An integer indicating the problem dimension : dim.
    2) A string to access the file containing the DN:
        geomfile.
    3) A float indicating the chain density used to 
        generate the DN: chain_density.
    4) A string indicating the chain behaviour: spring_type.
    5) A list containing the chain parameters other than 
        the chain length(s): params.
    6) An integer indicating the type of loading: loading.
    7) A float indicating the target final stretch: max_stretch.
    8) An integer with the number of stretch increments: Ninc.
    9) A string representing the folder in which the results
        of the simulation will be stored: results_folder.
    10) A boolean indicating if current configs should be 
        stored: visual_flag.
    
    """
    
    f = open(filename,'r')
    params = []
    
    f.readline()
    dim = int(f.readline())
    
    
    f.readline()
    geomfile = f.readline().strip('\n')
    
    f.readline()
    chain_density = float(f.readline())

    f.readline()
    spring_type = f.readline().strip('\n')

    f.readline()
    data = f.readline().split()

    for i in range(len(data)):
        params.append(float(data[i]))

    f.readline()
    loading = int(f.readline().strip('\n'))

    f.readline()
    data = f.readline().split()
    max_stretch = float(data[0])
    Ninc = int(data[1])
    
    f.readline()
    results_folder = f.readline().strip('\n')
    
    
    f.readline()
    data = f.readline().strip('\n');
    if data == 'y':
        visual_flag = True;
    else:
        visual_flag = False;
    
    f.close()

    
    return dim, geomfile, chain_density , spring_type, params, loading, max_stretch, Ninc, results_folder, visual_flag


def createFolder(folder_name):
    """
    Check if folder exists. If it does not, create it.
    
    folder_name: string representing the path of the folder.
    
    The function returns None.
    """
    
    # Remove trailing slash if it exists
    if folder_name.endswith('/'):
        folder_name = folder_name[:-1]
    
    # Check if the folder exists, if not, create it
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
        print(f"Folder '{folder_name}' created.")
    else:
        print(f"Folder '{folder_name}' already exists.")
    
    return


def printChainPar(model, chain_params, polydispersity_flag):
    """
    Print normalised chain parameters.
    
    model : string indicating the type of bond behaviour.
            model = '1': Gaussian
            model = '2': FJC
            model = '3': Extensible FJC
    chain_params: tuple with the chain parameters
    polydispersity_flag: boolean that is true if the DN
                        is polydisperse.
    
    The function returns none.
    """
    
    if model == '1': 
        print('spring type: harmonic (Gaussian)');
        if not polydispersity_flag:
            print('b = %g, N = %g' %chain_params)
        else:
            print('polydispersed DN')
            print('b = %g' %chain_params[0])
            
    elif model == '2': 
        print('spring type: FJC');
        if not polydispersity_flag:
            print('b = %g, N = %g' %chain_params)
        else:
            print('polydispersed DN')
            print('b = %g' %chain_params[0])
            
    elif model == '3': 
        print('spring type: Extensible FJC');
        if not polydispersity_flag:
            print('b = %g, N = %g, Eb = %g, critical_eng = %g' %chain_params)
        else:
            tmp = chain_params[0], chain_params[2], chain_params[3],
            print('polydispersed DN')
            print('b = %g, Eb = %g, critical_eng = %g' %tmp);
            
    
    return


def moveDatFile(original_file, dat_folder, inc, loading):
    """
    Copy current dat file and move the copy to destination 
    folder.
    
    original_file: name of the dat file to be copied.
    dat_folder: name of the folder the copied gile will
                be moved to.
    inc: integer indicating the current increment number.
    loading: integer indicating the loading type.
            loading = 1: uniaxial tension
            loading = 2: biaxial tension
            loading = 3: pure shear
    
    The function returns None.
    """
    
    # Initialise copy file name 
    if loading == 1:
        dat_file = 'uniaxial_inc' + str(inc) + '.dat'
    elif loading == 2:
        dat_file = 'biaxial_inc' + str(inc) + '.dat'
    else:
        dat_file = 'pshear_inc' + str(inc) + '.dat'
    
    # Check if destination folder exists
    createFolder(dat_folder);
    
    # Copy and move the copied file
    os.system('cp test.dat %s' %dat_file)
    os.system('mv %s %s' %(dat_file, dat_folder))
    
    return



def remove_initially_tooLong(model):
    """
    Scan non-equilibrated DN to detect too long chains
    
    Inputs:
        model (str): type of bond behaviour.
                    model = '1': Gaussian
                    model = '2': FJC
                    model = '3': Breakable extensible FJC
                    model = '4': Breakable FJC
                    
    Outputs:
        found_tooLong (bool): Flag returning True if too elongated bonds 
                                were found.
    """
    
    # Print initial message
    print('Scanning network to detect too long chains...')
    
    # Create DN object
    DN = NetworkClass ("test.dat", "", "main.in")
    
    # Detect too elongated chains
    tooLong_ids = DN.detect_tooLong_chains(model, cut_off = 0.95)
    if len(tooLong_ids) > 0:
        print("Detected bonds that are too elongated at the start!!")
        
        ## Remove from DNs the elongated chains
        print("Removing these bonds from the network...")
        DN.rewrite_data_file("test.dat", tooLong_ids, model)
        print("Done!")
        
        ## Set found_tooLong to True
        found_tooLong = True
        
    else:
        ## Print message and set found_tooLong = False
        print("No bonds are too elongated at the start.")
        found_tooLong = False
        
    
    print("Scanning complete!")
    
    return found_tooLong



def reduce_LAMMPS_timestep(main_file):
    """
    Reduce the time step of LAMMPS integrator.
    
    Inputs:
        main_file: name of the original LAMMPS input file.
        
    Outputs:
        None
    """
    
    # Name of temporary main_file
    temp_file = "small_dt.in"
    
    # Read all lines from the original main_file
    with open(main_file, "r") as f:
        lines = f.readlines()
        
    # Write the temporari file
    with open(temp_file, "w+") as f:
        for line in lines:
            if not "timestep" in line:
                f.write(line)
            elif "reset_timestep" in line:
                f.write(line)
            else:
                data = line.strip("\n").split("\t")
                deltaT = float(data[1]) / 10
                f.write("timestep\t%g\n" %deltaT)
        
    return





if __name__ == "__main__":
    main()
