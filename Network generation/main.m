% Network generation algorithm in matlab
% -------------------------------------------------------------------------
% This is a control code that calls the functions required to generate
% discrete networks (DN). The main feature of this code is the independent
% tuning of chain density and end-to-end distance distribution.
% author: Lucas Mangas Araujo
% email: lucas.mangasaraujo@jesus.ox.ac.uk
% date: 19/06/2024
% -------------------------------------------------------------------------
% For more details of the code, refer to the README.md in GitHub repository
% at https://github.com/LucasMangasAraujo/Discrete_Nets/tree/main
clc; close all; clear all
%% User inputs
% -------------------------------------------------------------------------
% Set the network parameters
dim = 3; % Dimension in space
n = 10e3; % Expected number of chains in the netwoek
expected_bond_length = 15; % Average end-to-end distance with units
var_bond_length = 2.7; % Variance end-to-end distance with units
nu = 1e-3; % Crosslink density
compensation_factor = 1.35; % Relaxation compensation factor
file_name = "test.txt"; % text file containing the network topology
folder_name = "Networks/"; % Path of  the folder containing the text files
chain_length = 100; % Chain length in number of Kuhn segments

%    0 1 2 3 4 5
D = [0 0 0 0 1 0 ]; % Functionality of the crosslinks
D = D/sum(D); 
average_degree = sum((0:length(D)-1).*D); % Average degree

% Set polydisperse parameters
polydisperse_flag = true; % Polydispersity flag
distribution = "Log"; % Distribution ID
distribution_par = {chain_length, 30}; % Distribution parameters


% Adjust the prescribed average end-to-end distance with the compensation
expected_bond_length = expected_bond_length * compensation_factor;
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Calculate some other quantities based in the user inputs
density = 2 * nu / average_degree; % density of crosslinks
N = 2 * n/ average_degree; % Number of crosslinks to be seeded in the computational domain
physical_length = power(n / nu, 1/dim); % length of the RVE
mu_d = expected_bond_length / physical_length; % Dimensionless average
s_d = var_bond_length / physical_length; % Initialise the nodes matrix and store information; % Dimensionless standard deviation
if dim == 2
    theoretical_xi = power(2* n, 1/2) * mu_d; 
else
    theoretical_xi = power(4 * n, 1/3) * mu_d / sqrt(3); 
end
% -------------------------------------------------------------------------

% Print initialisation of the code
fprintf("--------------------------------------------------------------\n")
fprintf("Generating network in unit box with the following parameters:\n")
fprintf("mu_d:  %g\n", mu_d);
fprintf("s_d:   %g\n", s_d);
fprintf("For this parameters, the as-generated xi is %g\n", theoretical_xi)
fprintf("--------------------------------------------------------------\n")
fprintf("\n")
%% Fixed parameters
% -------------------------------------------------------------------------
% Futher parameters used in the algorithm that are kept fixed
std_cutoff = 1; % Cut off the distribution
spacing_factor = 0.5; % Controls the radius of the Poisson disk
Boundaries   = 'periodic'; % Generation boundary conditions
PointProcess = 'poisson'; % Type of domain seeding
% -------------------------------------------------------------------------

%%  Generate the network
% -------------------------------------------------------------------------
% Call the function that creates the network
fprintf("--------------------------------------------------------------\n")
fprintf("Generating the network...\n");

[V , G, BoundaryEdges] = createGraph(N, D, density, dim, spacing_factor, expected_bond_length, var_bond_length, std_cutoff, PointProcess,Boundaries);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Create boundary nodes

% Extract the nodes, connections and boundary links in matrix form
[nodes, bonds, boundary_bonds] = networkInformation(V, G.Edges, BoundaryEdges, dim);

% Compute the length of relevant arrays
Nboundary = length(boundary_bonds); % Number of periodic chains
Ninternal  = length(bonds) - Nboundary; % Number of "internal" chains

% Creat the boundary nodes
[new_nodes, new_bonds, boundary_nodes] = createBoundaryNodes(nodes, bonds, Ninternal, Nboundary, dim);
    
% Create array with the chain length of each chain in the network 
bond_types = ones(length(new_bonds), 1) * chain_length;

fprintf("Finished generating the network...\n");
fprintf("--------------------------------------------------------------\n")
fprintf("\n")
% -------------------------------------------------------------------------

%% Write the text file
% -------------------------------------------------------------------------
fprintf("--------------------------------------------------------------\n")
fprintf("Writing text file...\n");
fprintf("Folder name: %s\n", folder_name);
fprintf("File name:   %s\n", file_name);
writeTxtFile(folder_name, file_name, new_nodes, new_bonds, boundary_nodes, bond_types)
fprintf("Done!\n");
fprintf("--------------------------------------------------------------\n")
fprintf("\n");
% -------------------------------------------------------------------------
%% Create the polydisperse twin of the generated network
% -------------------------------------------------------------------------
if polydisperse_flag
    % Generate the chain length distribution
    fprintf("--------------------------------------------------------------\n")
    fprintf("Generating polydispersed network...\n");
    chain_lengths = polydispersity(distribution, distribution_par, length(new_bonds));
    fprintf("Done!\n");
    fprintf("--------------------------------------------------------------\n")
    fprintf("\n");
    
    % Write the polydisperse text file
    temp = file_name.split(".");
    poly_file = temp(1) + "_" + distribution + ".txt";
    writeTxtFile(folder_name, poly_file, new_nodes, new_bonds, boundary_nodes, chain_lengths)

end
% -------------------------------------------------------------------------

%% Terminate
fprintf("Network generation done!!\n")