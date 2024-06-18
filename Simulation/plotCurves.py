"""
This script is used to plot the stress-stretch curves from DN simulations
"""


# Import packages for plotting
from math import *
import numpy as np
import sys, os
import seaborn as sns
import matplotlib.pyplot as plt 
import matplotlib.colors as mcolors # Colors available


# Import user-defined scripts
from utils import readParams, readGeometry

# Plot pre-sets
sns.set(rc={"figure.dpi":300})
sns.set_style("ticks")
plt.rc('text', usetex=True)  # enable LaTeX rendering
plt.rcParams["font.family"] = "serif"




def main():
    
    # Get the geometry file
    _, geomfile, _, _, _, loading, _, _, results_folder,_ = readParams('inputs.txt')
    
    
    # Decide whether the DN used is polydispersed
    _, _, _, BondTypes = readGeometry(geomfile);
    chain_lengths = np.array(list(BondTypes.values()))
    polydispersity_flag = not np.all(chain_lengths == chain_lengths[0])
    
    # Read the txt file with the results
    path_to_file = resultsFile(results_folder, loading, polydispersity_flag);
    data = np.loadtxt(path_to_file, skiprows=1)
    
    # Plot
    fig, ax = plt.subplots(1,1, figsize=(8, 8));
    stretch, stress = data[:,0], data[:,1] - data[:,-1];
    ax.scatter(stretch, stress, color = "mediumblue", s = 80, 
               label = 'DN results')
    
    ax.set_ylabel(r"$\frac{\sigma_1b^3}{kT}$", fontsize = 35);
    ax.set_xlabel(r"$\lambda_1$", fontsize = 35);
    ax.set_xlim(left = 1, right = ceil(stretch[-1]));
    ax.set_ylim(bottom = 0, top = 0.06)
    ax.tick_params(axis='both', which='major', labelsize = 26)
    ax.legend(fontsize = 20);
    
    plt.tight_layout()
    plt.savefig("Stress_stretch.svg")
    plt.tight_layout()
    plt.show()
    # breakpoint()
    return


def resultsFile(results_folder, loading, polydispersity_flag):
    """
    Return the name of the txt file with the data.
    
    results_folder: name of the folder where the text file is.
    loading: integer indicating the type of load.
    polydispersity_flag: Boolean informing if the DN is polydispersed
    
    The function returns:
    1) A string representing the path to the text file: path_to_file.
    """
    
    # Assemble 'sufix' 
    if loading == 1: tmp = '_uniaxial.txt';
    elif loading == 2: tmp = '_biaxial.txt';
    else: tmp = '_pshear.txt';
    
    if polydispersity_flag:
        path_to_file = results_folder + 'data' + tmp;
    else:
        path_to_file = results_folder + 'dataPoly' + tmp;
    
    
    return path_to_file


if __name__ == "__main__":
    main()