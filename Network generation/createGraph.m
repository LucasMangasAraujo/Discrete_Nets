%% createNetwork
% -------------------------------------------------------------------------
% This function creates a discrete network with prescribed end-to-end
% distance and crosslink connectivity distributions, and chain density. 
% Its main task is to call all the secondary functions needed for the 
% assembly process.
% 
% N: number of crosslinks in the network
% D: target crosslink connectivity distribution
% density: crosslink density
% dim: the dimension of the box, typically 2 or 3
% spacing_factor: the parameter used in the Poisson disk model, it defines the packing ratio
% expected_bond_length: target average of the end-to-end distance distribution
% var_bond_length: target standard deviation of the end-to-end distance distribution
% std_cutoff: the cut off on the edge distribution, for efficiency
% PointProcess: seeding method used to sample the unit box
% Boundaries: type of boundary conditions.
% 
% The function returns:
% 1) Array of nodes (crosslinks) coordinates in the unit box:  V
% 2) Graph object of the generated network: G
% 3) Sparse matrix representing the periodic chains: BoundaryEdges
% -------------------------------------------------------------------------

function [V, G, BoundaryEdges] = createGraph(N, D, density, dim, spacing_factor, expected_bond_length, var_bond_length, std_cutoff, PointProcess,Boundaries)
    % -------------------------------------------------------------------------
    % Call the function that creates the network

    [A V empirical_deree_distribution EdgeLength_list Unmathched_edges spacing box_scale x C] = GeometricGraph( N, D, density, dim, spacing_factor, expected_bond_length, var_bond_length, std_cutoff, PointProcess,Boundaries);
    disp('Target/Empirical degree distribtuion:')
    disp([ D; empirical_deree_distribution ])
    disp({'Number of matched edges:', sum(sum(A))/2})
    disp({'Number of unmatched edges:', Unmathched_edges})
    % -------------------------------------------------------------------------

    % -------------------------------------------------------------------------
    % Parameters for the visualisation of the network

    edge_alpha = 0.2; % Transparency of the edges
    marker_size = 6;
    % Show the obtained network and empirical bond length distribution
    [G, BoundaryEdges] = showNetwork(EdgeLength_list, x, box_scale , spacing , expected_bond_length, var_bond_length, V, A, edge_alpha, marker_size);
    % -------------------------------------------------------------------------

end


%% showNetwork

% -------------------------------------------------------------------------
% This function returns the graph object for further computations, and the
% list of periodic edges. It also shows graphically the generated network, 
% and compares target with empirical edge length distributions.
% 
% EdgeLength_list: Array of chains lengths in the unit box
% x: Array with discretisation of the distribution support
% box_scale: Length scale associated with the prescribed density
% spacing: Spacing between the vertices in the unit box
% spacing_factor: the parameter used in the Poisson disk model, it defines the packing ratio
% expected_bond_length: target average of the end-to-end distance distribution
% V: Array of nodes coordinates in the unit box:
% A: Adjacency matrix
% edge_alpha: Transparency of the edges
% marker_size: Size of the nodes
% 
% The function returns:
% 1) Graph object: G 
% 2) Sparse matrix representing the periodic chains: BoundaryEdges
% -------------------------------------------------------------------------

function [G, BoundaryEdges] =  showNetwork(EdgeLength_list, x, box_scale , spacing , expected_bond_length, var_bond_length, V, A, edge_alpha, marker_size)
    % -------------------------------------------------------------------------
    subplot(4,1,1:3)
    cla

    % Obtain the empirical distribution
    x = x * box_scale;
    spacing = spacing * box_scale;
    EdgeLength_list = EdgeLength_list * box_scale;
    
    % Computing the empirical end-to-end distance distribution
    h = hist( EdgeLength_list, x );
    h( end ) = 0;
    h = h/trapz( x, h );
    ch = cumtrapz(x,h);
    
    npdf = normpdf( x, expected_bond_length, var_bond_length );
    npdf( x<spacing ) = 0;
    npdf = npdf/trapz( x, npdf );
    ncdf = cumtrapz( x, npdf );
    emp_length = trapz( x, x.*h );
    % -------------------------------------------------------------------------

    % -------------------------------------------------------------------------
    % Plot the network and the distributions
    A( length( V ), length( V ) ) = 0;
    A = A|A';
    
    % D is a symmetric traceles matrix. D_ij = D_ji gives the distance between
    % points i and j
    D = dist( V' );
    
    % 
    BoundaryEdges = A & D>1/2;
    A = A&~BoundaryEdges;
    
    G = graph( A );
    plt = plot(G);
    axis square
    plt.NodeLabel='';
    
    plt.XData = V( :, 1 );
    plt.YData = V( :, 2 );
    
    xlim( [0 1] )
    ylim( [0 1] )
    
    xlabel("x", FontSize= 20, Interpreter="latex");
    ylabel("y", FontSize= 20, Interpreter="latex");
    zlabel("z", FontSize= 20, Interpreter="latex");
    
    if size(V,2)==3
        plt.ZData = V( :, 3 );
        view([45 45])
        camproj perspective
        zlim( [0 1] )
    end
    
    plt.EdgeColor = 'k';
    plt.Marker = '.';
    plt.LineWidth = 3;
    plt.EdgeAlpha = edge_alpha;
    plt.MarkerSize = marker_size;
    
    hold on
    set(gcf,'Color','w');
    
    subplot(4,1,4)
    
    
    
    tar=normpdf( x, expected_bond_length, var_bond_length );
    tar(x<spacing)=0;
    tar=tar/trapz(x,tar);
    
    
    ch=cumtrapz(x,h);
    
    ctar = normcdf( x, expected_bond_length, var_bond_length );
    ctar = (ctar-ctar(1));
    ctar=ctar/ctar(end);
    
    nrm = sqrt(trapz( x, (ch-ctar).^2 ) );
    
    plot( x, h, x, tar )
    xlim([ spacing*0.9 max(x) ] );
    xlabel('Bond length', FontSize= 16, Interpreter="latex");
    ylabel('Density', FontSize= 16, Interpreter="latex");
    
    % -------------------------------------------------------------------------

end