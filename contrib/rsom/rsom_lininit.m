function sMap = rsom_lininit(Distances, msize, varargin)

%RSOM_LININIT computes linear initialization for the RSOM. 
%
% sM = rsom_lininit(D, msize, [argID, value, ...])
%
%  sM = rsom_lininit(D, msize, 'lattice', 'hexa');
%
% Input and output arguments: 
%   D          (matrix) dissimilarity data for training, size nData x nData
%   msize      (vector) map size
%   [argID,    (string) See below.
%    value]    (varies) 
%
%   sM         (struct) RSOM map struct, the initialized map 
%
% Here are the valid argument IDs and corresponding values.
%   'lattice'     (string) map lattice: 'hexa' or 'rect'
%   'shape'       (string) map shape: 'sheet', 'cyl' or 'toroid'
%   'name'        (string) name of the RSOM struct
%
% For more help, try 'type rsom_lininit' or check out online documentation.
% See also RSOM_RANDINIT, RSOM_BATCHTRAIN.

%%%%%%%%%%%%% DETAILED DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% rsom_lininit
%
% PURPOSE
%
% Initializes the RSOM in the input space linearly.
%
% SYNTAX
%
%  sM = rsom_lininit(D, msize);
%  sM = rsom_lininit(...,'argID',value,...);
%
% DESCRIPTION
%
% The grid in the input space is initialized along the eigenvectors 
% corresponding to the largest eigenvalues. These are calculated by 
% classical MDS from the dissimilarity matrix D. If D is generated from
% euclidean data, squared distances should be provided i.e.
% (D)_ij = ||x_i - x_j||^2. 
% 
% REFERENCES
%
%   Barbara Hammer, Alexander Hasenfuss: Topographic Mapping of Large
%   Dissimilarity Data Sets. Neural Computation 22(9): 2229-2284 (2010)
%
% REQUIRED INPUT ARGUMENTS
%
%   D          (matrix) dissimilarity data for training, size nData x nData
%   msize      (vector) map size
%
% OPTIONAL INPUT ARGUMENTS 
%
%  argID (string) Argument identifier string (see below).
%  value (varies) Value for the argument (see below).
%
%  The optional arguments can be given as 'argID',value -pairs.
%  The valid IDs and corresponding values are listed below. 
%
%  Below is the list of valid arguments: 
%   'lattice'     (string) map lattice: 'hexa' or 'rect'
%   'shape'       (string) map shape: 'sheet', 'cyl' or 'toroid'
%   'name'        (string) name of the RSOM struct
%
% EXAMPLES
%
%   sM = rsom_lininit(D, [10 10]);
%   sM = rsom_lininit(D, [10 10], 'lattice', 'hexa');
%
% SEE ALSO
%
%   rsom_randinit    Initialize a RSOM randomly
%   rsom_batchtrain  Train a RSOM

% Contributed to SOM Toolbox vs2, December 7th, 2012 by Alexander Schulz
% Copyright (c) Alexander Schulz
% http://www.cis.hut.fi/projects/somtoolbox/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Prase input
p = inputParser;
p.addRequired('Distances', @isnumeric);
p.addRequired('msize', @isnumeric);

p.addParamValue('lattice', 'rect', @ischar);
p.addParamValue('shape', 'sheet', @ischar);
p.addParamValue('name', '', @ischar);

p.parse(Distances, msize, varargin{:});

lattice = p.Results.lattice;
shape   = p.Results.shape;
name    = p.Results.name;
clear p;
mdim     = length(msize);
nData    = size(Distances,1);

%% Init and normalize
% create the topology struct
sTopol = som_topol_struct('msize', msize, ...
                          'lattice', lattice, ...
                          'shape', shape);


Coords = som_unit_coords(msize,'rect','sheet');
% normalize the neuron positions to mean 0
Coords = bsxfun(@minus, Coords, mean(Coords,1));

% compute mds embedding of the distance matrix
[V, e]=cmdscale(Distances.^(1/2));
% keep only the first mdim dimensions
V = V(:, 1:mdim);
e = e(1:mdim);
scale = std(V);

% scale Coords to the scale of the data
%Coords = bsxfun(@rdivide, Coords, std(Coords));
%Coords = bsxfun(@times, Coords, scale);
Coords = bsxfun(@rdivide, Coords, std(Coords)+1e-10);

% project neurons onto the data
% e is n * variance
standardDev = sqrt(e/nData);
V = bsxfun(@rdivide, V, nData*standardDev'+1e-10);
%V = bsxfun(@rdivide, V, sqrt(sum(V.^2, 2))+1e-10);

cCodebook = Coords * V';
cCodebook = cCodebook + 1/nData;

%% create the resulting struct
sTrain = som_train_struct('algorithm','lininit');
           
sMap = struct('type', 'rsom_map', ...
              'cCodebook', cCodebook, ...
              'topol', sTopol, ...
              'labels', cell(1), ...
              'neigh', 'gaussian', ...
              'trainhist', cell(1), ...
              'name', name);

sTrain = som_set(sTrain,'time',datestr(now,0));
sMap.trainhist = sTrain;

