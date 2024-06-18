function [A V empirical_deree_distribution EdgeLength_list Unmathched_edges spacing box_length x C hh f_emp] = GeometricGraph( nPts, degree_distribution, density, dim, spacing_factor, expected_bond_length, var_bond_length, std_cutoff, Process, broundaries)
%  This function constructs a geometric network with a given degree
%  distribution in N dimensions. The coordinates of the vertices are
%  taken from the Poisson disk model.

%  nPts number of ponts.
%  DENSITY  is the number of points in a unit box
%  DIM   is the dimension of the box, typically 2 or 3
%  EXPECTED_BOND_LENGTH and VAR_BOND_LENGTH are the mean and variance of
%  the edge length distribution
%  SPACING_FACTOR is the parameter used in the Poisson disk model, it defines packing ratio

%
% The function returns:
% 1) Adjacency matrix A, 
% 2) the array of vertex coordinates  V
% 3) empirical_deree_distribution and an array of edge lengths.
% 4) number of edges the algorithm could not satisfy.

%% Parameters %%
% number of consecutive failures before termination.
max_tries = 300; 

%% Precomputation %%
% number of vertices



% spacing for Poisson disk model, for the given density in the unit cube
spacing  = ( 1/density )^( 1/dim ) * spacing_factor; 

Volume = nPts/density;

% size of the box
box_length = Volume^( 1/dim ); 


% scalling all the parameters to fit into unit box.
spacing = spacing/box_length;
expected_bond_length = expected_bond_length/box_length;
var_bond_length = var_bond_length/box_length;
%cutoff = cutoff/box_length;

% cut off on the edge distribution, for efficiency
cutoff = std_cutoff*sqrt(var_bond_length); %the standard deviation away from mean


%% Vertices
if strcmp(Process,'uniform')
    V = rand( nPts, dim );
    spacing=0;
elseif  strcmp(Process,'poisson')
    

    % N.B.: this goes into a very long loop at dense shpere packings
    V=0;

    while abs(size( V, 1 )/nPts-1)>0.04 % we allow the number of vertices be 4% away from the target number.
        sizeI = ones( 1, dim ) * 2;
        V   = poissonDisc( sizeI, spacing, nPts );
        V   = V - 1;
    end


else
    error('No point process specified.')
end

%% Distances between vertices in V
if strcmp(broundaries,'periodic') 
    [ D, i, j ] = dist_periodic( V, expected_bond_length+cutoff, spacing, 1 );
elseif  strcmp(broundaries,'nonperiodic') 
    [ D, i, j ] = dist_nonperiodic(V,expected_bond_length+cutoff, spacing);
else
    error('No boundary type specified.')
end

%% Edge length distribution, with importance sampling bias factor %%

x = linspace( max(min(expected_bond_length-cutoff,spacing),0), expected_bond_length+cutoff, 400);
h = normpdf( x, expected_bond_length, var_bond_length );
h = h/trapz( x, h );

f_emp = hist( D, x );
f_emp = f_emp / trapz( x, f_emp );
hf = h./f_emp;
hh = h;
hf( h==0 | isinf(hf)| isnan(hf) )=0;


C = max(hf);

p = interp1( x, hf, D );
p(isnan(p))=0;
p = p/sum(p);


%% Randomized sequential construction %%
EdgeLength_list=[];
A = sparse(length(V),length(V),0);

%normalisation to be sure
degree_distribution = degree_distribution/sum( degree_distribution );

%generate graphic degree sequence
%deg=[1 0]; % the while loop is not so relevant for large graphs.
%while ~isgraphic(deg)
    deg = weightedSample( length(V), 0:length(degree_distribution)-1, degree_distribution );
%end
deg=deg';
DEG=deg;

%number of edges
M = sum(deg) / 2;

Unmathched_edges = 0;
pp =  p .* deg( i ) .* deg( j ).*( 1 - DEG(i).*DEG(j)/4);

% Place edges one by one:
for k=1:M
   
    
    if all( pp==0 )
        %terminate, no good edges left   
        break
    end

    id = weightedSample( 1, 1:length(pp), pp );
    
    % do we have edge #id already?
    counter=0;
    %{0
    while  A( i( id ), j( id ) ) == 1
       id = weightedSample( 1, 1:length(pp), pp );
        counter = counter+1;
        if counter == max_tries %give up
            break
        end
    end
    %}
   
    
    if counter < max_tries
        
        % place an edge
        A( i( id ), j( id ) ) = 1;
        A( j( id ), i( id ) ) = 1;
        
        %list of used distances
        EdgeLength_list( end+1 ) = D( id );
        
        F = unique([find(i==i(id))' find(j==j(id))' find(j==i(id))' find(i==j(id))']);
        
        % update degree sequence
        deg( i(id) ) = deg( i(id) ) - 1;
        deg( j(id) ) = deg( j(id) ) - 1;

        p(id) = 0;
        
        pp(F) =  p(F) .* deg( i(F) ) .* deg( j(F) ) .*( 1 - DEG(i(F)).*DEG(j(F))/4);
        
    
    end

end;

Unmathched_edges = sum( deg )/2;
% empirical degree distribution
empirical_deree_distribution = histomap( sum( A ), 0:length( degree_distribution ) - 1 );
empirical_deree_distribution = empirical_deree_distribution / sum( empirical_deree_distribution );


end
%%


function s = weightedSample( n, P, W ) 
% INPUTs: number of draws from a discrete distribution (n)
%         possible values to pick from, (P)
%         set of normalized weights/probabilities, (W)
% OUTPUTs: s - set of n numbers drawn from P
%              according to the weights in W
    s = [];

    W = W(:)' / sum(W);

    unit = [0 cumsum(W)];
    
    while length( s )<n
      s = [ s P( find( unit<rand, 1, 'last' ) ) ];         
    end 
    
end


function B = isgraphic(seq)
% Check whether a sequence of number is graphical via Erdos-Gallai formula
% INPUTs: a sequence (vector) of numbers
% OUTPUTs: boolean, true or false

if not(isempty(find(seq<=0))) | mod(sum(seq),2)==1
    % there are non-positive degrees or their sum is odd
    B = false; return;
end

n=length(seq);
seq=-sort(-seq);  % sort in decreasing order

for k=1:n-1
    sum_dk = sum(seq(1:k));
    sum_dk1 = sum(min([k*ones(1,n-k);seq(k+1:n)]));
    
    if sum_dk > k*(k-1) + sum_dk1; B = false; return; end

end
B = true;
end

function h=histomap(X,b)
%compute simple histogram
    h=[];
    for i=b
        h(end+1)=sum(X==i);
    end;
end


function [pts] = poissonDisc(sizeI,spacing,nPts,showIter)

% Purpose:
% N-dimensional poisson disc sampling function. This can also be used to
% randomly sample k pts from N-dimensional space with a minimum separation
% distance.
%
% Inputs:
% sizeI -   [required] Size of volume from which points are to be 
%           sampled
% spacing - [required] Minimum separation distance between points
% nPts -    [Default is 0] if nPts = 0 For poisson disc sampling.
%           nPts = k to sample k-pts from N-dimensional space with
%           minimum separation distance
% showIter - [Default is 0] If showIter == 1, this option can be used to 
%            see how points are generated through each iteration. It can be
%            useful when code is taking a long time to generate points. 
%
% Output:
% pts - All eligible points
% M.Patel, 2016

%%%%%%% Initial parameters setup
% Parsing inputs and setting default values
if nargin == 3; showIter = 0; end
if nargin == 2; showIter = 0; nPts = 0; end

% Setting properties for iterations
ndim = length(sizeI);   % Number of Dimensions
k = 5;  % Number of 'dart' tries in each grid.
dartFactor = 4; %Select number of sample data in each iterations. Change it to
% reduce run time for code. Have to play around with number. 


%%%%%%% Making Grid read for iterations
%Make grid size such that there is just one pt in each grid
dm = spacing/sqrt(ndim);    % grize cell size [Bridson 2007]

%Make Grid
for i = 1:ndim
    sGrid{1,i} = 1:dm:sizeI(i);
end
[sGrid{:}] = ndgrid(sGrid{:});
sizeGrid = size(sGrid{1});

% Convert Grid points into a nx3 array;
for i = 1:ndim
    sGrid{i} = sGrid{i}(:);
end
sGrid = cell2mat(sGrid);

% Arrays to show eligible grids for dart throws and keeping score of darts
% thrown in a particular grid
emptyGrid = logical(ones(size(sGrid,1),1)); %Eligible Grids
nEmptyGrid = sum(emptyGrid);    %Number of eligible Grids
scoreGrid = zeros(size(emptyGrid)); %Score of darts thrown in Grid

% Darts to be thrown per iterations
% This hugely influences speed of the algorithm. Change dartFactor for it. 
if nPts == 0
    nPts = nEmptyGrid;
    ndarts = round(nEmptyGrid/dartFactor);
end
ndarts = round(nPts/dartFactor);

%%%%%%%%% Iterative process to generate points
% Initialize parameters
ptsCreated = 0;
pts = [];
iter = 0;

% Start Iterative process
tic
while ptsCreated<nPts & nEmptyGrid >0
   % if length(pts)>1
   %     plot(pts(:,1),pts(:,2),'o')
   % end
    %Thrown darts in eligible grids
    availGrid = find(emptyGrid == 1);   %Eligible grids for dart throw
    dataPts = min([nEmptyGrid,ndarts]); % Darts to be thrown
    p = datasample(availGrid,dataPts,'Replace',false); %Select grids for darts
    tempPts = sGrid(p,:) + dm*rand(length(p),ndim); %Dart throw!!!
    
    
    % Find good dart throws
    [~,D] = knnsearch([pts;tempPts],tempPts,'k',2); %Finding distance between all darts(pts)
    D = D(:,2); 

    withinI = logical(prod(bsxfun(@lt,tempPts,sizeI),2)); %Eligible pts should be withing sizeI 
    eligiblePts = withinI & D>spacing; %elgible pts should also have minimum separation distance
    
    scorePts = tempPts(~eligiblePts,:); %Keep score from bad dart throws :(
    tempPts = tempPts(eligiblePts,:);   % Save good dart throws :)
    
    
    %Update empty Grid
    emptyPts = floor((tempPts+dm-1)/dm);
    emptyPts = num2cell(emptyPts,1);
    emptyIdx = sub2ind(sizeGrid,emptyPts{:});
    emptyGrid(emptyIdx) = 0;
    
    %Update score pts
    scorePts = floor((scorePts+dm-1)/dm);
    scorePts = num2cell(scorePts,1);
    scoreIdx = sub2ind(sizeGrid,scorePts{:});
    scoreGrid(scoreIdx) = scoreGrid(scoreIdx) + 1;
    
    %Update emptyGrid if scoreGrid has exceeded k dart throws
    emptyGrid = emptyGrid & (scoreGrid<k);
    
    %Update quantities for next iterations
    nEmptyGrid = sum(emptyGrid);
    pts = [pts;tempPts];
    ptsCreated = size(pts,1);
    iter = iter+1;
    ttoc = toc;
    
    %Display iteration details
    if showIter == 1
        disp(sprintf('Iteration: %d    Points Created: %d    EmptyGrid:%d    Total Time: %0.3f',iter,ptsCreated,nEmptyGrid,ttoc));
    end
    
end

% Cut down pts if more points are generated
if size(pts,1)>nPts
    p = 1:size(pts,1);
    p = datasample(p,nPts,'Replace',false);
    pts = pts(p,:);
end

end


function [D,i,j] = dist_periodic(X,cutoff,spacing,box_scale)
x_size = box_scale;

D = 0;

for d = 1:size(X,2)

    x = X( :, d );

    dx = x - x';
    dx = triu(dx);
    dx( 1:length( dx ) + 1:end ) = 0;
    
    dx = dx(:);

    dx( dx > x_size * 0.5  ) = dx( dx >   x_size * 0.5 ) - x_size;
    dx( dx <= -x_size * 0.5) = dx( dx <= -x_size * 0.5 ) + x_size;
    
    D = D + dx.^2;
    
end
D = sqrt(D);

[ i j ] = meshgrid( 1:size(X,1) );
i = i(:);
j = j(:);

F = D ==0 | D > cutoff | D<spacing; % Discard distances below the Poisson distacne, larger 
i( F ) = [];
j( F ) = [];
D( F ) = [];


end



function [D,i,j] = dist_nonperiodic(V,cutoff,spacing)

D = dist( V' );
D = triu(D);
D( 1:length( D ) + 1:end ) = 0;
[ i j ] = meshgrid( 1:length(D) );
i = i(:);
j = j(:);
D = D(:);

% filtering out large distances for efficiancty
F = D ==0 | D > cutoff | D<spacing;
i( F ) = [];
j( F ) = [];
D( F ) = [];

end

