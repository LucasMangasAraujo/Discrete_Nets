%% writeTxtFile
% -------------------------------------------------------------------------
% This function writes the text file that can be used in DN simulations
% 
% folder: folder that will receive the txt file
% file: name of the destination text file
% nodes: array containing the nodes coordinates
% bonds: array witht connection between the nodes
% boundary_nodes: array with the IDs of the nodes at the boundary
% bond_types: array with the chain lenght of each chain in the network
% 
% The function has no returns
% -------------------------------------------------------------------------

function writeTxtFile(folder, file, nodes, bonds, boundary_nodes, bond_types)
    % -------------------------------------------------------------------------
    % Check the existance of the folder, and assamble the path to the file
    folderCheck(folder)
    file_path = folder + file;

    % Write the nodes
    f = fopen(file_path, "w+");
    fprintf(f,'$nodes\n');
    for i = 1:length(nodes)
        fprintf(f,'%d, %.7e, %.7e, %.7e \n',i, nodes(i,1), nodes(i, 2), nodes(i, 3));
    end
    
    % Write the connections
    fprintf(f,'$bonds\n');
    for i = 1:length(bonds)
        fprintf(f,'%d, %d, %d \n',i , bonds(i,2), bonds(i,1));
    end
    
    % Write the boundary nodes
    fprintf(f,'$boundary\n');
    
    for i = 1:length(boundary_nodes)
        fprintf(f,"%d ",boundary_nodes(i));
    end
    fprintf(f,"\n");
    
    % Write the bond types
    fprintf(f,"$BondTypes\n");
    
    % Print information
    for i = 1:length(bond_types)
        fprintf(f,'%d, %g \n',i,bond_types(i));
    end
    
    % Close the file
    fclose(f); 
    
    % -------------------------------------------------------------------------

end

%% Function for the creation of the destination folder
% -------------------------------------------------------------------------
% This function creates the folder that will contain the generated text
% file containing the architecture of the polymer network.
% folder_name: string with the full path of the folder
% 
% The function has no returns
% -------------------------------------------------------------------------

function folderCheck(folder_name)
% -------------------------------------------------------------------------
    % Check if the folder exists
    if ~isfolder(folder_name)
        % If the folder does not exist, create it
        mkdir(folder_name);
        fprintf('Folder "%s" created successfully.\n', folder_name);
    else
        fprintf('Folder "%s" already exists.\n', folder_name);
    end
% -------------------------------------------------------------------------
end