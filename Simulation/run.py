from math import *
import numpy as np
import sys, os

# Import utils
import utils

# Import NetworkClass (new_feature)
from network_class import NetworkClass

'''
    Control code of the DN simulations. 
'''

def main(): 
    
    #standalone script: read parameters from file
    inputfile = 'inputs.txt';
    
    # Read simulation parameters
    dim, geomfile, chain_density, model, chain_params, loading, max_stretch, Ninc, results_folder, visual_flag = utils.readParams(inputfile)
    
    # run simulation
    print(100 * "=")
    print('Starting simulation...')
    runsim(dim, geomfile, chain_density, model, chain_params, loading, max_stretch, Ninc, results_folder, visual_flag);
    print('Completed!')
    print(100 * '=')
    
    return



def runsim(dim, geomfile, chain_density, model, chain_params, loading, max_stretch, Ninc, results_folder, 
            visual_flag = False):
    """
    Run the DN simulation and generate files with simulation results.
    
    dim : integer indicating the dimension of the problem.
    geomfile: string representing the path to the txt file 
                containing the DN architecture.
    chain_density: float representing the as-generated chain 
                    density.
    model: string indicating the type of bond behaviour.
            model = '1': Gaussian
            model = '2': FJC
            model = '3': Breakable extensible FJC
            model = '4': Breakable FJC
    chain_params: list containing chain parameters other than 
            the chain length.
    loading: integer indicating the type of load.
            loading = 1: uniaxial tension
            loading = 2: biaxial tension
            loading = 3: pure shear
    max_stretch: float indicating the target final stretch
    Ninc: integer with the number of stretch increments
    results_folder: tring representing the folder in which 
                    the results of the simulation will be stored.
    visual_flag: Optional boolean flag indicating if the current 
                    configurations should be stored for visualisation
                    Default: False.
    
    The function returns None.
    """
    # Generate input files for LAMMPS
    print('Generating input files for LAMMPS...')
    mainfile = 'main.in'
    posfile = 'test.dat'
    
    # Print inputs
    print('dim: %d' %dim)
    print('max stretch: %g' %max_stretch)
    print('number of incs: %d' %Ninc);
    
    # Read geometry file
    Nodes, Bonds, Boundary, BondTypes = utils.readGeometry(geomfile);
    Nnodes, Nbonds, Nboundary = len(Nodes), len(Bonds), len(Boundary);
    
    # Check if DN is polydispersed
    BondTypes = {idx: chain_params[1] for idx in BondTypes.keys()}
    chain_lengths = np.array(list(BondTypes.values()))
    polydispersity_flag = not np.all(chain_lengths == chain_lengths[0])
    
    # Normalise the Kuhn length by the RVE's length
    nub3 = chain_density * chain_params[0];
    bKuhn_normalised = pow(nub3 / (2 * (Nnodes - Nboundary ) ), 1/3); ## Additional boundary bonds
    
    
    # And assemble chain params array in dimensionless format
    chain_params_normalised = np.copy(chain_params);
    chain_params_normalised[0] = bKuhn_normalised;
    utils.printChainPar(model, tuple(chain_params_normalised), polydispersity_flag)
    if polydispersity_flag:
        mean_N = np.mean(chain_lengths);
        print("mean chain length: %g" %mean_N)
    
    
    # Print number of chains
    print('total number of nodes: %d' %Nnodes);
    print('total number of bonds: %d' %Nbonds)
    print('number of boundary nodes: %d' %Nboundary)
    
    # Write initial position file for LAMMPS
    utils.writePositions(posfile, Nodes, Bonds, Boundary , BondTypes, model, chain_params_normalised)
    
    # Create loading history
    stretches = np.linspace(1, max_stretch, Ninc)
    
    # Run minisation with breakable bonds to remove unrealistic chains
    if (model == "2" or polydispersity_flag):
        utils.writeMain(mainfile,posfile,Boundary,dim,"2");
        
    else:
        utils.writeMain(mainfile,posfile,Boundary,dim,model)
        
    
    # Initialise arrays
    Fall, Sall = [], [];
    
    # Run the simulation
    for i, stretch in enumerate(stretches):
        print(100 * '-')
        
        # Run increment
        inc = i; 
        print('###Inc %d' %inc);
        
        if i == 0: ## initial relaxation
            Nbonds_start = len(Bonds); 
            F = np.ones(3)
            print('Running initial relaxation...')
            err = runinc(loading,inc,0,dim);
            
            if err:
                ## Scan for bonds that are too elongated from the start
                print('Found error in initial relaxation!!')
                found_tooLong = utils.remove_initially_tooLong(model)
                
                breakpoint()
                if found_tooLong:
                    print('Running initial relaxation ...')
                    err = runinc(loading,inc,0,dim);
                else:
                    print("Convergence issues stem from unknown reasons.")
                    inc=inc-1
                    break
            else:
                print('Done!')
            
            
            # Obtain reordered Bonds dictionary
            Bonds, typeIDs = getUpdatedBonds('test.dat');
            Nbonds = len(Bonds);
            if Nbonds_start != Nbonds:
                print('The first minisation step broke %g bonds.' %(Nbonds_start - Nbonds));
            else:
                print('No broken bonds in the first step.')
            
            # Calculate initial squared end-to-end distance and pre-stretch
            r2, pre_stretch2 = getDist('test.dat', Bonds, BondTypes, typeIDs, bKuhn_normalised, polydispersity_flag);
            print('r02: %g, r0: %g' %(r2, sqrt(r2)))
            print('lambda02: %g, lambda0: %g' %(pre_stretch2, sqrt(pre_stretch2)))
            print('affine modulus in b^3/kT units: %g' %(nub3 * pre_stretch2))
            
            # Rewrite LAMMPS input file if needed
            if (model == "2" or polydispersity_flag):
                utils.writeMain(mainfile,posfile,Boundary,dim,model);
            
        else:
            inc_stretch = stretch - stretches[i - 1];
            F = defgrad(F,inc_stretch,loading,dim)
            print('Running increment...')
            err = runinc(loading,inc,inc_stretch,dim)
            if err:
                print('Found error, exit')
                inc=inc-1
                break
            else:
                print('done!')
            
        
        # Calculate the stress
        S = calculateStress('test.res',dim);
        S *= pow(bKuhn_normalised, 3); ## b^3/kT units
        
        print('F: %g, %g, %g' %(F[0],F[1],F[2]))
        print('Sb^3/kT: %g, %g, %g' %tuple(S))
        Fall.append(F[0]);
        Sall.append(S);
        
        # Move dat file if needed
        if visual_flag:
            print('Copying current configs...')
            utils.moveDatFile('test.dat', 'Current_geometries', inc, loading)
            print('Done!')
        
        print('Moving to next increment...')
        print(100* '-')
        
        
    print('Simulation completed!')
    
    # Write the txt file with the stress-stretch data
    print(100 * '*')
    print('Writing file with results...')
    utils.createFolder(results_folder);
    tup = np.array(Fall), np.array(Sall);
    matrix = np.column_stack(tup, )
    
    if loading == 1: tmp = '_uniaxial.txt';
    elif loading == 2: tmp = '_biaxial.txt';
    else: tmp = '_pshear.txt';
    
    
    if polydispersity_flag:
        path_to_file = results_folder + 'dataPoly' + tmp;
    else:
        path_to_file = results_folder + 'data' + tmp;
        
    
    print('stress and stretch data will be written to %s' %path_to_file.split('/')[-1]);
    header = 'lambda S1E S2E S3E'
    np.savetxt(path_to_file, matrix, delimiter = " ", header = header);
    print('Done!')
    print(100 * '*')
    
    return


def runinc(loading,inc,dl,dim):
    """ 
    Run the strech increment using LAMMPS.
    
    loading: integer indicating the type of load.
        loading = 1: uniaxial tension
        loading = 2: biaxial tension
        loading = 3: pure shear
    inc: integer indicating the icrement number.
    dl: float representing the increment of stretch.
    dim: integer indicating tre problem dimension.
    
    The funtion returns
    1) A boolean that is True if the LAMMPS simulation
        was aborted.
    """

    #displacement increment applied to the box of dimension 1.2 
    du = 1.2*dl
    
    #copy main file in temporary file 
    os.system('cp main.in main_tmp.in')

    if dim==3:
    
        #uniaxial loading
        if loading == 1:
            newline = "fix 1 all deform 1 x delta " + str(-du/2.) + " " + str(du/2)  + " y volume z volume remap x units box\n"

        #equi-biaxial loading
        elif loading == 2:
            newline = "fix 1 all deform 1 x delta " + str(-du/2.) + " " + str(du/2) + " y delta " + str(-du/2) + " " + str(du/2) + " z volume remap x units box\n"

        #pure shear
        elif loading == 3:
            newline = "fix 1 all deform 1 x delta " + str(-du/2) + " " + str(du/2) + " y volume z delta 0 0 remap x units box\n"
        
        # Free-Sweelling
        elif loading == 4:
            newline = "fix 1 all deform 1 x delta " + str(-du/2.) + " " + str(du/2) + " y delta " + str(-du/2) + " " + str(du/2) + " z delta " + str(-du/2) + " " + str(du/2) + " remap x units box\n";

        
        else:
            print('invalid loading in runinc: %d' %loading)
            exit()

    #2D simulation
    else:
        newline = "fix 1 all deform 1 x delta " + str(-du/2) + " " + str(du/2) + " y volume remap x units box\n"

    #rewrite the main file by replacing the line with the loading
    fin = open('main_tmp.in','r')
    fout = open('main.in','w')

    for line in fin:
        
        if 'delta' in line:
            fout.write(newline)
        else:
            fout.write(line)

    fin.close()
    fout.close()

    #run lammps
    os.system('~/.local/bin/lmp -in main.in > log')

    #check for error in log file
    err = checkerror('log.lammps')

    return err



def run_reduced_dt(main_file):
    """
    Run FIRE with smaller dt
    
    Inputs:
        main_file: original lammps input file
        
    Output:
        err: boolean informing the success or failure of the 
                     minimisation
    """
    
    # Generate neem input file
    utils.reduce_LAMMPS_timestep(main_file)
    
    # Call lammps
    os.system('~/.local/bin/lmp -in small_dt.in > log')
    
    # Check if smaller dt led to convergence
    err = checkerror('log.lammps')
    breakpoint()
    return err 


def checkerror(filename):
    """
    Access LAMMPS log file and check if errors ocurred.
    
    filename: name of the LAMMPS log file.
    
    The function returns True if an error was found.
    """
    
    fin = open(filename,'r')

    for line in fin:
        
        if 'ERROR' in line:
            print('Error message in log.lammps:')
            print(line)
            return True
    return False


def getUpdatedBonds(filename):
    """
    Get the bonds new ordering after the first energy 
    minisation in LAMMPS.
    
    filename: string representing the dat file generated
                by LAMMPS.
                
    The function returns:
    1) Updated dictionary of bonds: Bonds.
    2) Dictionary with the type of each bond: typeIDs.
    """
    
    # Initialize the outputs
    Bonds, typeIDs = {}, {}
    
    # Open file and read
    with open(filename,"r") as f:
        
        # Read until the bonds were reached
        key = f.readline().strip("\n");
        
        while "Bonds" not in key:
            key = f.readline().strip("\n");
            
        f.readline().split()

        data = f.readline().split()
        while len(data) > 1:
            idx = int(data[0])
            bondType = int(data[1]);
            n1 = int(data[2])
            n2 = int(data[3])
            typeIDs[idx] = bondType;
            Bonds[idx] = [n1,n2] 
            

            data = f.readline().split()
        
    
    return Bonds, typeIDs


def getDist(filename, Bonds, BondTypes, typeIDs, bKuhn_normalised, polydispersity_flag):
    
    '''
    calculate distances in a relaxed network from 'test.dat'.
    
    filename: output dat file from LAMMPS after relaxation.
    Bonds: dictionary containing a list of the network bonds.
    BondTypes: dictionary containing the chain lengths in the 
                network.
    typeIDs: list of integers of each bond type after LAMMPS
            reordered the bonds.
    bKuhn_normalised: normalised Kuhn length.
    polydispersity_flag: boolean indicating if the DN is polydispersed.
    
    The function returns
    1) A float representing the mean squared end-to-end distance (in a unit box):
        r2.
    2) A float representing the mean squared pre-stretch: pre_stretch_2.
    '''

    f = open(filename,'r')
    
    Nbond = len(Bonds)

    #create a dictionary with current positions of nodes
    Coords = {} #dictionary with current positions

    key = f.readline().strip('\n')

    while 'Atoms' not in key:
        key = f.readline().strip('\n')
    
    f.readline()
    
    data = f.readline().split()
    while len(data) > 1:
        idx = int(data[0])
        x = float(data[3])
        y = float(data[4])
        z = float(data[5])
        Coords[idx] = [x,y,z] 
    
        data = f.readline().split()
    
    
    # Check if network is polydispersed
    if not polydispersity_flag:
        N = BondTypes[1]; 
    
    #calculate relaxed distances 
    i=0
    r2 = 0.
    pre_stretch2 = 0.;
    

    for idx,bond in Bonds.items():
        # Calculate distance
        n1 = bond[0]
        n2 = bond[1]
        dx = Coords[n1][0]-Coords[n2][0]
        dy = Coords[n1][1]-Coords[n2][1]
        dz = Coords[n1][2]-Coords[n2][2]
        dist2 = dx**2 + dy**2 + dz**2
        r2 += dist2
        
        # Calculate pre-stretch
        if polydispersity_flag:
            N = BondTypes[typeIDs[idx]];
        
        pre_stretch2 += dist2 / (N * pow(bKuhn_normalised, 2))

        i += 1

    #average squared distance
    
    r2 /= len(Bonds);
    pre_stretch2 /= len(Bonds);
    
    return r2, pre_stretch2


def defgrad(F,dl,loading,dim):
    
    """
    Calculate the diagonalised deformation gradient
    
    F: numpy array containing the past deformation 
        gradient.
    dl: float with the stretch increment 
        in direction 1.
    loading: integer indicating the type of load.
        loading = 1: uniaxial tension
        loading = 2: biaxial tension
        loading = 3: pure shear
        
    dim: integer indicating problem dimension.
    
    The function returns
    1) A numpy array with the updated deformation
        gradient: F.
    """

    F[0] += dl

    if dim == 3:

        #uniaxial
        if loading == 1:
            F[1] = 1./sqrt(F[0])
            F[2] = F[1]

        #equibiaxial
        elif loading == 2:
            F[1] = F[0]
            F[2] = 1./(F[1]*F[1])       

        #pure shear
        elif loading == 3:
            F[1] = 1./F[0]
            F[2] = 1.
        
        # Free-Sweelling
        elif loading == 4:
            F[1] = F[0];
            F[2] = F[0];
        
        
        else:
            print("invalid loading type: %d" %loading)

    #2D: only one loading mode
    else:

        F[1] = 1./F[0]
        F[2] = 1

    return F


def calculateStress(filename,dim):

    """ Read test.res and calculate the stress from reaction forces on boundary nodes.
    
    filename: name of file containing the reaction forces in the boundary.
    dim: integer indicating problem dimension.
    
    The function returns:
    1) An array containing (true) stress components in L^3/kT units: S.
    """

    f = open(filename,'r')
    S = np.zeros(3)

    key = f.readline().strip('\n')

    while 'id' not in key:
        key = f.readline().strip('\n')
    
    data = f.readline().split()
    while len(data) > 1:

        if dim == 3:
            x = float(data[2])
            y = float(data[3])
            z = float(data[4])
            fx = -float(data[5])    #minus sign to take the reaction force
            fy = -float(data[6])
            fz = -float(data[7])
            
            
            S[0] += fx*x
            S[1] += fy*y
            S[2] += fz*z

        else:
            x = float(data[2])
            y = float(data[3])
            fx = -float(data[4])    #minus sign to take the reaction force
            fy = -float(data[5])

            S[0] += fx*x
            S[1] += fy*y

                
        data = f.readline().split()

    return S 



if __name__ == "__main__":
    main();
