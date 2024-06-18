%% createBoundaryNodes
% -------------------------------------------------------------------------
% This function creates the boundary nodes based on the list of periodic
% chains of network.
% 
% nodes: matrix containing the coordinates of the nodes in the unit box
% bonds: matrix containing the connections between the nodes
% Ninternal: Number of non-periodic chains
% Nboundary: Number of periodic chains
% dim: Dimension of the box
% 
% The function returns
% 1) A matrix containing coordinates of original and boundary nodes: new_nodes
% 2) A matrix containing the updated connections between nodes in the network: new_bonds
% 3) An array with the IDs of the boundary nodes: boundary_nodes
% -------------------------------------------------------------------------

function [new_nodes, new_bonds, boundary_nodes] = createBoundaryNodes(nodes, bonds, Ninternal, Nboundary, dim)
    % -------------------------------------------------------------------------
    % Each periodic link generates two boundary nodes
    Nnodes_new = length(nodes) + 2 * length(Nboundary);
    Nbonds_new = length(bonds) + 2 * length(Nboundary);
    Nboundary_nodes = 2 * length(Nboundary);
    
    % Preallocate the output arrays
    new_nodes = zeros(Nnodes_new, 3);
    new_bonds = zeros(Nbonds_new, 2);
    boundary_nodes = zeros(Nboundary_nodes, 1);
    
    % Pre-allocate local arrays
    boundary_nodes_coordinates = zeros(Nboundary_nodes, 3);
    created_bonds = zeros(Nboundary_nodes, 2);
    % -------------------------------------------------------------------------
    
    % -------------------------------------------------------------------------
    % Create the boundary nodes
    Nnodes = length(nodes); % As-generated number of nodes
    for i = (Ninternal + 1):length(bonds)
        id1 = bonds(i, 1); id2 = bonds(i,2); % IDs of the connected nodes
        a = nodes(id1, :)'; b = nodes(id2, :)'; % Coordinates of these nodes
        local_counter = 2 * (i - Ninternal); % Pair of boundary nodes
        boundary_id1 = local_counter - 1 + Nnodes;  boundary_id2 = local_counter + Nnodes; % IDs of the pair of boundary nodes
        
        boundary_nodes_coordinates(local_counter - 1: local_counter, :) = findBoundaryNodes(a, b, dim);
        created_bonds(local_counter - 1: local_counter, :) = [boundary_id1, id1; boundary_id2, id2];
        boundary_nodes(local_counter - 1: local_counter, :) =  [boundary_id1; boundary_id2];
        
    end
    
    % Update the nodes and bonds matrices
    new_nodes = [nodes; boundary_nodes_coordinates];
    new_bonds = [bonds; created_bonds];
    % -------------------------------------------------------------------------

end