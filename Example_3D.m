dim = 3;
N   = 1000;

density = 200; % points per unit square

var_bond_length = 0.05; % in a unit box
expected_bond_length =  0.15; 
std_cutoff = 1;

Boundaries   = 'periodic';
PointProcess = 'poisson';
spacing_factor  = 0.78; %used in the poisson disk model;

%    0 1 2 3 4 5
D = [0 0 0 0 0 1 ];
D = D/sum(D); 

%for plotting:
edge_alpha = 0.2;
%make network:
[A V empirical_deree_distribution EdgeLength_list Unmathched_edges spacing box_scale x C] = GeometricGraph( N, D, density, dim, spacing_factor, expected_bond_length, var_bond_length, std_cutoff, PointProcess,Boundaries);

%rescale parameters for output:
x = x * box_scale;
spacing = spacing * box_scale;
EdgeLength_list = EdgeLength_list * box_scale;


disp('Target/Empirical degree distribtuion:')
disp([ D; empirical_deree_distribution ])
disp({'Number of matched edges:', sum(sum(A))/2})
disp({'Number of unmatched edges:', Unmathched_edges})



h = hist( EdgeLength_list, x );
h( end ) = 0;
h = h/trapz( x, h );
ch = cumtrapz(x,h);


npdf = normpdf( x, expected_bond_length, var_bond_length );
npdf( x<spacing ) = 0;
npdf = npdf/trapz( x, npdf );
ncdf = cumtrapz( x, npdf );
emp_length = trapz( x, x.*h );



subplot(4,1,1:3)



cla

A( length( V ), length( V ) ) = 0;
A = A|A';

D = dist( V' );
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
plt.MarkerSize = round(N/120);
hold on
set(gcf,'Color','w');
%plot(V(:,1),V(:,2),'r.','MarkerSize', 10);
%plt.NodeColor = [.9 0.3 0.3];

subplot(4,1,4)



tar=normpdf( x, expected_bond_length, var_bond_length );
tar(x<spacing)=0;
tar=tar/trapz(x,tar);


ch=cumtrapz(x,h);

ctar = normcdf( x, expected_bond_length, var_bond_length );
ctar = (ctar-ctar(1));
ctar=ctar/ctar(end);

nrm = sqrt(trapz( x, (ch-ctar).^2 ) )

plot( x, h, x, tar )
xlim([ spacing*0.9 max(x) ] );
xlabel('Bond length')
ylabel('Density')




