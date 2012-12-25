function [sMap, distNeuronsToData, sTrain] = rsom_batchtrain(sMap, Distances, varargin)

%RSOM_BATCHTRAIN Use a batch algorithm to train the Self-Organizing Map on
%relational data.
%
% [sM,distToData,sT] = som_batchtrain(sM, D, [argID, value, ...])
%
%  sM = som_batchtrain(sM, D);
%  [sM, distToData] = som_batchtrain(sM, D, 'epochs', 50, 'lambda0', 2.5);
%
% Input and output arguments: 
%   sM         (struct) rsom map struct, the trained and updated map is 
%                       returned
%   D          (matrix) dissimilarity data for training, size nData x nData
%   [argID,    (string) See below.
%    value]    (varies) 
%
%   distToData (matrix) distance  from the neurons to data, 
%                       size nNeurons x nData
%   sT         (struct) learning parameters used during the training
%
% Here are the valid argument IDs and corresponding values.
%   'lambda0'     (scalar) initial decay constant
%   'lambdaFin'   (scalar) final decay value
%   'epochs'      (scalar) number of epochs
%
% For more help, try 'type rsom_batchtrain' or check out online documentation.
% See also RSOM_RANDINIT, RSOM_LININIT.

%%%%%%%%%%%%% DETAILED DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% rsom_batchtrain
%
% PURPOSE
%
% Trains a Relational Self-Organizing Map on dissimilarity data using the
% batch algorithm. 
%
% SYNTAX
%
%  sM = rsom_batchtrain(sM,D);
%  sM = rsom_batchtrain(...,'argID',value,...);
%
% DESCRIPTION
%
% Trains the given and already initialized RSOM with the batch algorithm.
% Training data are expected to be given in the form of dissimilarity data
% D. If D is generated from euclidean data, i.e.
% (D)_ij = ||x_i - x_j||^2, training the relational SOM is equivalent to
% the standard SOM and the neuron position in the input space can be
% acquired by sMap.cCodebook * x, where x are training data. 
%
% REFERENCES
%
%   Barbara Hammer, Alexander Hasenfuss: Topographic Mapping of Large
%   Dissimilarity Data Sets. Neural Computation 22(9): 2229-2284 (2010)
%
% REQUIRED INPUT ARGUMENTS
%
%   sM         (struct) initialized rsom map struct
%   D          (matrix) dissimilarity data for training, size nData x nData
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
%   'lambda0'     (scalar) initial decay constant
%   'lambdaFin'   (scalar) final decay value
%   'epochs'      (scalar) number of training epochs
%
% OUTPUT ARGUMENTS
% 
%   sM         (struct) the trained rsom map struct. The current training
%                       is added to the training history (sM.trainhist).
%   distToData (matrix) 
%   sT         (struct) train struct; information of the accomplished 
%                       training
%
% EXAMPLES
%
%  mData = size(D,1);
%  sM = rsom_randinit(nData, [20 3]); % create a 20 x 3 grid
%  sM = rsom_batchtrain(sM,D);
%
% SEE ALSO
%
%   rsom_randinit   Initialize a RSOM randomly

% Contributed to SOM Toolbox vs2, December 7th, 2012 by Alexander Schulz
% Copyright (c) Alexander Hasenfuss
% http://www.cis.hut.fi/projects/somtoolbox/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Prase input
p = inputParser;
p.addRequired('sMap', @isstruct);
p.addRequired('Distances', @isnumeric);

% Get size
m = size( Distances, 1 );
GridDistances = som_unit_dists(sMap.topol).^2;
n = size( GridDistances, 1 );

p.addParamValue('lambda0', n/2, @isnumeric);
p.addParamValue('lambdaFin', 0.01, @isnumeric);
p.addParamValue('epochs', 100, @(x) ~(x-floor(x)) & (length(x) == 1));

p.parse(sMap, Distances, varargin{:});
lambda0   = p.Results.lambda0;
lambdaFin = p.Results.lambdaFin;
epochs    = p.Results.epochs;

%% Initialize
lambda = lambda0 * (lambdaFin/lambda0).^((0:(epochs-1))/(epochs-1));

QuantizationErrors = zeros(1,epochs);

cNeurons = sMap.cCodebook;

sTrain = som_train_struct('algorithm','batch');
sTrain = som_set(sTrain, 'neigh', sMap.neigh, 'trainlen', epochs, ...
    'radius_ini', sqrt(lambda0/2), 'radius_fin', sqrt(lambdaFin/2));

%% Perform optimization
for i=1:epochs,

  % Determine quadratic distances between neurons and data
   distNeuronsToData = determine_distance_relational_neurons_to_data( Distances, cNeurons );
   
  % Determine winner w_I*(x_j) for every datapoint (w* := [w_I*(x_1), ..., w_I*(x_m)]) 
   [err, w_star] = min( distNeuronsToData, [], 1 );

  % Determine k~_ij (distances to winner neurons)
   KK = GridDistances(:,w_star);
  
  % Neighbourhood Function 
   HH = exp(-(KK)/lambda(i));  

  % Update rule
   cNeurons = HH ./ (sum(HH,2) * ones(1,m));

  % Determine quantization error
   QuantizationErrors(i) = determine_qerror( sqrt( distNeuronsToData ) );
   
end

% Determine distance between neurons and data
distNeuronsToData = sqrt( determine_distance_relational_neurons_to_data( Distances, cNeurons ) );
sMap.cCodebook = cNeurons;

sTrain = som_set(sTrain,'time',datestr(now,0));
tl = length(sMap.trainhist);
sMap.trainhist(tl+1) = sTrain;


function distNeuronsToData = determine_distance_relational_neurons_to_data( Distances, cNeurons )

%
% Determines dissimilarities between relational neurons and data on which they are defined
%
% Distances ... (m x m)-Matrix of dissimilarities between data points
% cNeurons ... Coefficient (n x m)-Matrix of relational neurons
%
% Use squared distances to compare with other methods on euclidean distances !
%
% A.Hasenfuss (c) 2007-2008 (Optimized Version)
%

[n,m] = size( cNeurons );

% Determine distances between neurons and data
  temp1 = cNeurons * Distances;
 
  temp2 = zeros(n,1);
  for i = 1:n
    temp2(i,1) = 0.5 * (temp1(i,:) * cNeurons(i,:)');
  end
 
  distNeuronsToData = temp1 - temp2*ones(1,m);
  
  
  
function QuantizationError = determine_qerror( distNeuronsToData )

%
% Determines the quantization error for median neurons
%
% A.Hasenfuss (c) 2006
%

[WinnerDist, WinnerIdx] = min( distNeuronsToData, [], 1  );
QuantizationError = sum( WinnerDist.^2 )/2;
  