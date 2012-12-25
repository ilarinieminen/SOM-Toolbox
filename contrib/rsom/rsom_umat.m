function U = rsom_umat(sMap, Dists, mode)

%RSOM_UMAT computes the U-matrix defined by the SOM.
%
% U = rsom_umat(sMap, D, mode)
%
%  U = rsom_umat(sMap, D, 'mean');
%
% Input and output arguments: 
%   sMap       (struct) RSOM map struct
%   D          (matrix) dissimilarity matrix used for training
%   [mode]     (string) 'min','mean','median','max', default is 'median'
%
%   U          (matrix) u-Matrix 
%
% For more help, try 'type rsom_umat' or check out online documentation.
% See also RSOM_SHOW.

%%%%%%%%%%%%% DETAILED DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% rsom_umat
%
% PURPOSE
%
% Calculates the U-matrix.
%
% SYNTAX
%
%  U  = rsom_umat(sMap, D);
%  U  = rsom_umat(sMap, D, 'mean');
%
% DESCRIPTION
%
% Calculates the U-matrix which can be used for visualization.
% 
% REFERENCES
%
%   Barbara Hammer, Alexander Hasenfuss: Topographic Mapping of Large
%   Dissimilarity Data Sets. Neural Computation 22(9): 2229-2284 (2010)
%
%   Iivarinen, J., Kohonen, T., Kangas, J., Kaski, S., "Visualizing 
%   the Clusters on the Self-Organizing Map", in proceedings of
%   Conference on Artificial Intelligence Research in Finland,
%   Helsinki, Finland, 1994, pp. 122-126.
%
% REQUIRED INPUT ARGUMENTS
%
%   sMap       (struct) RSOM map struct
%   D          (matrix) dissimilarity matrix used for training
%
% OPTIONAL INPUT ARGUMENTS 
%
%  mode     (string) 'min','mean','median','max', default is 'median
%
% EXAMPLES
%
%   U  = rsom_umat(sMap, D, 'mean');
%
% SEE ALSO
%
%   rsom_show    Visualizes a calcualted U-matrix

% Contributed to SOM Toolbox vs2, December 7th, 2012 by Alexander Schulz
% Copyright (c) Alexander Schulz
% http://www.cis.hut.fi/projects/somtoolbox/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin < 3)
    mode = 'median';
end

M       = sMap.cCodebook;
lattice = sMap.topol.lattice;
y       = sMap.topol.msize(1);
x       = sMap.topol.msize(2);

ux = 2 * x - 1; 
uy = 2 * y - 1;
U  = zeros(uy, ux);

calc = sprintf('%s(a)',mode);

mask = ones(x*y,1);
if size(mask,2)>1, mask = mask'; end


%% u-matrix computation

% Calculate the distances between all neurons
neuronDists = determine_distance_relational_neurons_to_neurons(Dists, M, sMap.topol.msize);

% distances between map units

if strcmp(lattice, 'rect'), % rectangular lattice
  
  for j=1:y, for i=1:x,
      if i<x, 
    dx = neuronDists(j,i,j,i+1); % horizontal
	U(2*j-1,2*i) = sqrt(dx);
      end 
      if j<y, 
    dy = neuronDists(j,i,j+1,i); % vertical
	U(2*j,2*i-1) = sqrt(dy);
      end
      if j<y & i<x,	 
    dz1 = neuronDists(j,i,j+1,i+1); % diagonals
    dz2 = neuronDists(j+1, i, j, i+1);
	U(2*j,2*i) = (sqrt(dz1)+sqrt(dz2))/(2 * sqrt(2));
      end
    end
  end

elseif strcmp(lattice, 'hexa') % hexagonal lattice  %%%% TODO change %%%%

  for j=1:y, 
    for i=1:x,
      if i<x,
    dx = neuronDists(j,i,j,i+1); % horizontal
	U(2*j-1,2*i) = sqrt(dx);
      end
      
      if j<y, % diagonals
    dy = neuronDists(j,i,j+1,i);
	U(2*j,2*i-1) = sqrt(dy);	
	
	if rem(j,2)==0 & i<x,
      dz = neuronDists(j,i,j+1,i+1);
	  U(2*j,2*i) = sqrt(dz);
	elseif rem(j,2)==1 & i>1,
      dz = neuronDists(j,i,j+1,i-1);
	  U(2*j,2*i-2) = sqrt(dz);
	end
      end
    end
  end
  
end

% values on the units

if (uy == 1 | ux == 1),
  % in 1-D case, mean is equal to median 

  ma = max([ux uy]);
  for i = 1:2:ma,
    if i>1 & i<ma, 
      a = [U(i-1) U(i+1)]; 
      U(i) = eval(calc);
    elseif i==1, U(i) = U(i+1); 
    else U(i) = U(i-1); % i==ma
    end
  end    

elseif strcmp(lattice, 'rect')

  for j=1:2:uy, 
    for i=1:2:ux,
      if i>1 & j>1 & i<ux & j<uy,    % middle part of the map
	a = [U(j,i-1) U(j,i+1) U(j-1,i) U(j+1,i)];        
      elseif j==1 & i>1 & i<ux,        % upper edge
	a = [U(j,i-1) U(j,i+1) U(j+1,i)];
      elseif j==uy & i>1 & i<ux,       % lower edge
	a = [U(j,i-1) U(j,i+1) U(j-1,i)];
      elseif i==1 & j>1 & j<uy,        % left edge
	a = [U(j,i+1) U(j-1,i) U(j+1,i)];
      elseif i==ux & j>1 & j<uy,       % right edge
	a = [U(j,i-1) U(j-1,i) U(j+1,i)];
      elseif i==1 & j==1,              % top left corner
	a = [U(j,i+1) U(j+1,i)];
      elseif i==ux & j==1,             % top right corner
	a = [U(j,i-1) U(j+1,i)];
      elseif i==1 & j==uy,             % bottom left corner
	a = [U(j,i+1) U(j-1,i)];
      elseif i==ux & j==uy,            % bottom right corner
	a = [U(j,i-1) U(j-1,i)];
      else
	a = 0;
      end
      U(j,i) = eval(calc);
    end
  end

elseif strcmp(lattice, 'hexa')
  
  for j=1:2:uy, 
    for i=1:2:ux,
      if i>1 & j>1 & i<ux & j<uy,      % middle part of the map
	a = [U(j,i-1) U(j,i+1)];
	if rem(j-1,4)==0, a = [a, U(j-1,i-1) U(j-1,i) U(j+1,i-1) U(j+1,i)];
	else a = [a, U(j-1,i) U(j-1,i+1) U(j+1,i) U(j+1,i+1)]; end       
      elseif j==1 & i>1 & i<ux,        % upper edge
	a = [U(j,i-1) U(j,i+1) U(j+1,i-1) U(j+1,i)];
      elseif j==uy & i>1 & i<ux,       % lower edge
	a = [U(j,i-1) U(j,i+1)];
	if rem(j-1,4)==0, a = [a, U(j-1,i-1) U(j-1,i)];
	else a = [a, U(j-1,i) U(j-1,i+1)]; end
      elseif i==1 & j>1 & j<uy,        % left edge
	a = U(j,i+1);
	if rem(j-1,4)==0, a = [a, U(j-1,i) U(j+1,i)];
	else a = [a, U(j-1,i) U(j-1,i+1) U(j+1,i) U(j+1,i+1)]; end
      elseif i==ux & j>1 & j<uy,       % right edge
	a = U(j,i-1);
	if rem(j-1,4)==0, a=[a, U(j-1,i) U(j-1,i-1) U(j+1,i) U(j+1,i-1)];
	else a = [a, U(j-1,i) U(j+1,i)]; end
      elseif i==1 & j==1,              % top left corner
	a = [U(j,i+1) U(j+1,i)];
      elseif i==ux & j==1,             % top right corner
	a = [U(j,i-1) U(j+1,i-1) U(j+1,i)];
      elseif i==1 & j==uy,             % bottom left corner
	if rem(j-1,4)==0, a = [U(j,i+1) U(j-1,i)];
	else a = [U(j,i+1) U(j-1,i) U(j-1,i+1)]; end
      elseif i==ux & j==uy,            % bottom right corner
	if rem(j-1,4)==0, a = [U(j,i-1) U(j-1,i) U(j-1,i-1)];
	else a = [U(j,i-1) U(j-1,i)]; end
      else
	a=0;
      end
      U(j,i) = eval(calc);
    end
  end
end




function neuronDist = determine_distance_relational_neurons_to_neurons(Dist, cNeurons, msize)


y = msize(1);
x = msize(2);

neuronDist = zeros(y,x,y,x);

for j=1:y
    for i=1:x
        neuronDist(j,i,:,:) = reshape(determine_distance_relational_1neuron_to_neurons(Dist, cNeurons, j + (i-1)*y), [y,x]);
    end
end




function dists = determine_distance_relational_1neuron_to_neurons(Dist, cNeurons, ind)

% Calculate the distance from the neuron on the position ind to all others

n = size(cNeurons, 1);
dists = zeros(n,1);
neuron = cNeurons(ind,:);

for i=1:n
    dists(i) = cNeurons(i,:)*Dist*neuron' - 0.5*neuron*Dist*neuron' - 0.5*cNeurons(i,:)*Dist*cNeurons(i,:)';
end






