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
    (dim, geomfile, chain_density, model, chain_params, loading, max_stretch, Ninc, 
        results_folder, visual_flag, rate_independent_scission) = utils.readParams(inputfile)
    
    
    # run simulation
    print(100 * "=")
    print('Starting simulation...')
    runsim(dim, geomfile, chain_density, model, chain_params, loading, max_stretch, Ninc, 
            results_folder, visual_flag, rate_independent_scission);
    print('Completed!')
    print(100 * '=')
    
    return



def runsim(dim, geomfile, chain_density, model, chain_params, loading, max_stretch, Ninc, results_folder, 
            visual_flag = False, rate_independent_scission = False):
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
    rate_independent_scission: Optional boolean informing whether 
                    rate-independent scissions are enabled
    The function returns None.
    """
    # Generate input files for LAMMPS
    print('Generating input files for LAMMPS...')
    mainfile = 'main.in'
    posfile = 'test.dat'
    
    # Print some inputs
    print('dim: %d' %dim)
    print('max stretch: %g' %max_stretch)
    print('number of incs: %d' %Ninc)
    
    # Read geometry file
    Nodes, Bonds, Boundary, BondTypes = utils.readGeometry(geomfile);
    Nnodes, Nbonds, Nboundary = len(Nodes), len(Bonds), len(Boundary);
    
    # Check if DN is polydispersed
    chain_lengths = np.array(list(BondTypes.values()))
    polydispersity_flag = not np.all(chain_lengths == chain_lengths[0])
    
    # Print number of chains
    print('total number of nodes: %d' %Nnodes);
    print('total number of bonds: %d' %Nbonds)
    print('number of boundary nodes: %d' %Nboundary)
    
    # Check if rate-independent scissions are active and if chain model is appropriate
    if rate_independent_scission:
        print('Rate-independent scissions are enabled')
        is_consistent = check_consistent_chain_model(model, rate_independent_scission)
        if not is_consistent: return
    
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
    
    
    # Write initial position file for LAMMPS
    utils.writePositions(posfile, Nodes, Bonds, Boundary , BondTypes, model, chain_params_normalised)
    
    # Create loading history
    stretches = np.linspace(1, max_stretch, Ninc)
    
    # Write LAMMPS input file
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
                
                if found_tooLong:
                    print('Running initial relaxation ...')
                    err = runinc(loading,inc,0,dim);
                else:
                    print("Convergence issues stem from unknown reasons.")
                    inc=inc-1
                    break
            else:
                print('Done!')
            
            
            # If convergence was reached, create DN object
            DN = NetworkClass ("test.dat", "test.res", "main.in")
            
            # Obtain reordered Bonds dictionary
            _, Bonds = DN.get_nodes_and_bonds()
            Nbonds = len(Bonds);
            if Nbonds_start != Nbonds:
                print('The first minisation step broke %g bonds.' %(Nbonds_start - Nbonds));
            else:
                print('No broken bonds in the first step.')
            
            
            # Calculate initial squared end-to-end distance and pre-stretch
            pre_stretch2 = np.square(DN.get_avg_preStretch(model))
            rms_r0 = DN.get_rms_r0()
            print(f'r02: {np.square(rms_r0)}, r0_rms: {rms_r0}')
            print(f'lambda02: {pre_stretch2}, lambda0: {sqrt(pre_stretch2)}')
            print(f'affine modulus in b^3/kT units: {nub3 * pre_stretch2}')
            
            # Store the pristine number of chains in case scissions are enabled
            if rate_independent_scission:
                Nbonds_pristine = Nbonds
                
            
        
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
            
        
        # Break chains if rate-independent scissions were enabled
        breakpoint()
        # Update DN object if not initial relaxation step
        if i != 0:
            DN = NetworkClass ("test.dat", "test.res", "main.in")
        
        # Calculate the stress
        S = DN.calculate_stress(dim)
        S *= pow(bKuhn_normalised, 3); ## b^3/kT units
        
        
        print('F: %g, %g, %g' %(F[0],F[1],F[2]))
        print('Sb^3/kT: %g, %g, %g' %tuple(S))
        Fall.append(F[0])
        Sall.append(S)
        
        
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
    
    if loading == 1: tmp = '_uniaxial.csv';
    elif loading == 2: tmp = '_biaxial.csv';
    else: tmp = '_pshear.csv';
    
    
    if polydispersity_flag:
        path_to_file = results_folder + 'dataPoly' + tmp;
    else:
        path_to_file = results_folder + 'data' + tmp;
        
    
    print('stress and stretch data will be written to %s' %path_to_file.split('/')[-1]);
    header = 'lambda S1E S2E S3E'
    np.savetxt(path_to_file, matrix, delimiter = ",", header = header);
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



def break_chains():
    '''
    Break chains in the network that have reached to the scission thredhold
    '''
    
    
    return

def check_consistent_chain_model(model, rate_independent_scission):
    '''
    Cross check if the chain model makes sense when rate-independent
    scissions are enabled.
    
    Inputs:
        model (str): chain model used.
            model = '1': Gaussian
            model = '2': FJC
            model = '3': Breakable extensible FJC
            model = '4': Breakable FJC
        
        rate_independent_scission (bool): True if chain scissions are enabled.
        
    Output:
        is_consistent (bool): True if chain model is appropriate.
    '''
    print('Checking if chain model selected is consistent with rate-independent scission...')
    if int(model) > 2:
        print('Chain model is consistent, proceed!')
        is_consistent = True
    else:
        print('Chains model is not consistent!!. Please, check inputs.txt.')
        print('Breaking simulation')
        is_consistent = False
    
    return is_consistent


if __name__ == "__main__":
    main();
