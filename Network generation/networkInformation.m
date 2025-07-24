%% networkInformation
% -------------------------------------------------------------------------
% This function returns information about the network that is useful for
% writing text files and creating boundary nodes. 
% 
% V: Network nodes and their coordinates
% InternalEdges: table with non-periodic or "internal" connections
% BoundaryEdges: Boundary chains     
% dim: Space dimension
% 
% This function returns:
% 1) Matrix with x, y and z coordinates of the nodes: nodes
% 2) Matrix with all connections in the graph: bonds
% 3) Matrix with the periodic bonds: boundary_bonds
% -------------------------------------------------------------------------

function [nodes, bonds, boundary_bonds] = networkInformation(V, InternalEdges, BoundaryEdges, dim)
    % -------------------------------------------------------------------------
    % Initialise the nodes matrix and store information
    if dim == 3
        nodes = V;
    else
        nodes = [V zeros(length(V), 1)];
    end
    
    % Generate matrix with the boundary connections
    [rows, columns] = find(triu(BoundaryEdges)); % Find non-zero entries in upper triangular part
    boundary_bonds = [columns, rows];
    
    
    % Create matrix with all the bonds in the network
    bonds = [table2array(InternalEdges); boundary_bonds];
    
    % -------------------------------------------------------------------------
end