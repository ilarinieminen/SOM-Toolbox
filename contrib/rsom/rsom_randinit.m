function sMap = rsom_randinit(nData, msize, varargin)

%RSOM_RANDINIT Initialize a relational Self-Organizing Map with random
%values.
%
% sMap = rsom_randinit(nData, msize, [argID, value, ...])
%
%  sMap = rsom_randinit(nData, [10 10])
%  sMap = rsom_randinit(nData, [20 3], 'lattice', 'hexa');
%
% Input and output arguments: 
%   nData      (scalar) number of training data
%   mize       (vector) number neurons per dimension
%   [argID,    (string) See below.
%    value]    (varies) 
%
%   sMap       (struct) initialized RSOM map struct
%
% Here are the valid argument IDs and corresponding values.
%   'lattice'     (string) map lattice: 'hexa' or 'rect'
%   'shape'       (string) map shape: 'sheet', 'cyl' or 'toroid'
%   'name'        (string) name of the map
%
% For more help, try 'type rsom_randinit' or check out online documentation.
% See also RSOM_LININIT, RSOM_BATCHTRAIN.

%%%%%%%%%%%%% DETAILED DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% rsom_randinit
%
% PURPOSE
%
% Initializes a Relational Self-Organizing Map with random values. 
%
% SYNTAX
%
%  sM = rsom_randinit(nData, msize);
%  sM = rsom_randinit(...,'argID',value,...);
%
% DESCRIPTION
%
% Computes a random initialization of the grid.
% 
% REFERENCES
%
%   Barbara Hammer, Alexander Hasenfuss: Topographic Mapping of Large
%   Dissimilarity Data Sets. Neural Computation 22(9): 2229-2284 (2010)
%
% REQUIRED INPUT ARGUMENTS
%
%   nData      (scalar) number of training data
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
%   nData = size(DissimMat, 1);
%   sM = rsom_randinit(nData, [10 10]);
%   sM = rsom_randinit(nData, [10 10], 'lattice', 'hexa');
%
% SEE ALSO
%
%   rsom_lininit     Initialize a RSOM linearly
%   rsom_batchtrain  Train a RSOM

% Contributed to SOM Toolbox vs2, December 7th, 2012 by Alexander Schulz
% Copyright (c) Alexander Schulz
% http://www.cis.hut.fi/projects/somtoolbox/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Prase input
p = inputParser;
p.addRequired('nData', @(x) ~(x-floor(x)) & (length(x) == 1));
p.addRequired('msize', @isnumeric);

p.addParamValue('lattice', 'rect', @ischar);
p.addParamValue('shape', 'sheet', @ischar);
p.addParamValue('name', '', @ischar);

p.parse(nData, msize, varargin{:});

lattice = p.Results.lattice;
shape   = p.Results.shape;
name    = p.Results.name;
clear p;
nNeurons = prod(msize);

%% Init and normalize
cNeurons = rand(nNeurons,nData);
cNeurons = cNeurons ./ (sum(cNeurons,2) * ones(1,nData));


%% create structs
sTopol = som_topol_struct('msize', msize, ...
                          'lattice', lattice, ...
                          'shape', shape);

sTrain = som_train_struct('algorithm','randinit');
           
sMap = struct('type', 'rsom_map', ...
              'cCodebook', cNeurons, ...
              'topol', sTopol, ...
              'labels', cell(1), ...
              'neigh', 'gaussian', ...
              'trainhist', cell(1), ...
              'name', name);

sTrain = som_set(sTrain,'time',datestr(now,0));
sMap.trainhist = sTrain;