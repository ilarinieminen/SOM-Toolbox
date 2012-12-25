function [cNeurons, distNeuronsToData, QuantizationErrors] = rneural_gas( Distances, n, epochs, lambda0 )

%RNEURAL_GAS quantizes the data space using the Relational Batch Neural Gas (RNG)
%algorithm.
%
% [cNeurons, distToData, QuantErr] = rneural_gas(Dists, nNeurons, nEpochs, [lambda0])
%
%  cNeurons = rneural_gas(Dists, 15, 20)
%  [cNeurons, distToData] = rneural_gas(Dists, 1, 20, 10)
%
% Input and output arguments ([]'s are optional):
%   Dists           (matrix) pairwise dissimilarities
%   nNeurons        (scalar) number of neurons
%   nEpochs         (scalar) number of training epochs
%   [lambda0]       (scalar) initial decay constant
%
%   cNeurons        (matrix) coefficients for neuron representation
%   distToData      (matrix) distances between neurons and data points
%   QuantErr        (vector) quantization error for each epoch
%
% For more help, try 'type rneural_gas' or check out online documentation.
% See also NEURAL_GAS, KMEANS_CLUSTERS.

%%%%%%%%%%%%% DETAILED DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% rneural_gas
%
% PURPOSE
%
% Trains a set of neurons usding the Relational Batch Neural Gas algorithm.
%
% SYNTAX
%
%   cNeurons = rneural_gas(Dists, nNeurons, nEpochs);
%   [cNeurons, distToData, QuantErr] = rneural_gas(Dists, nNeurons, nEpochs, [lambda0]);
%
% DESCRIPTION
%
% Relational Batch Neural Gas (RNG) performs a vector quantization of a 
% data set which is represented by dissimilarities. If euclidean data are
% available, applying RNEURAL_GAS on the Distance matrix Dist where 
% (Dist)_ij = ||x_i - x_j||^2 is equivalent to standard Batch Neural Gas.
% 
% REFERENCES
%
%   Barbara Hammer, Alexander Hasenfuss: Topographic Mapping of Large
%   Dissimilarity Data Sets. Neural Computation 22(9): 2229-2284 (2010)
%
% REQUIRED INPUT ARGUMENTS
% 
%   Dists           (matrix) (m,m) matrix of pairwise dissimilarities; 
%                            for euclidean data use squared distances
%   nNeurons        (scalar) number of neurons
%   nEpochs         (scalar) number of training epochs
%
% OPTIONAL INPUT ARGUMENTS 
%
%   [lambda0]       (scalar) initial decay constant (default=nNeurons/2)
%
% OUTPUT ARGUMENTS
%
%   cNeurons        (matrix) (nNeurons x m) coefficients for neuron 
%                            representation
%   distToData      (matrix) (nNeurons x m) matrix with distances between 
%                            neurons and data points
%   QuantErr        (vector) (1 x epochs) vector with quantization error 
%                            for each epoch
% 
% EXAMPLES
%
%   nNeurons = 15;
%   nEpochs  = 20;
%   cNeurons = rneural_gas(Dists, nNeurons, nEpochs);
%
% SEE ALSO
%
%   neural_gas        Use the Neural Gas algorithm for training.
%   kmeans_clusters   Use the k-means algorithm for vector quantization.

% Contributed to SOM Toolbox vs2, December 3rd, 2012 by Alexander Schulz
% Copyright (c) Alexander Hasenfuss
% http://www.cis.hut.fi/projects/somtoolbox/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Check parameters
%error(nargchk(3, 4, nargin));
narginchk(3,4);

% Get sizes
m = size( Distances, 1 );

% Set parameters
if nargin<4 || isempty(lambda0) || isnan(lambda0), lambda0 = n/2; end

% lambda
lambda = lambda0 * (0.01/lambda0).^((0:(epochs-1))/epochs);

%% Initialize
HH = zeros(n,m);
Ranking = zeros(n,m);
QuantizationErrors = zeros(1,epochs);

% Init and normalize
cNeurons = rand(n,m);
cNeurons = cNeurons ./ (sum(cNeurons,2) * ones(1,m));

% Do it...
for i=1:epochs

  % Determine dissimilarities between neurons and data
   distNeuronsToData = determine_distance_relational_neurons_to_data( Distances, cNeurons );        
  
  % Determine ranking (k_ij)
   [err, idx] = sort( distNeuronsToData, 1 );
   for k = 1:m
     Ranking(idx(:,k),k) = (1:n)';
   end

  % Calculate h(rank)
   Hvector = exp(-(0:(n-1))/lambda(i));  

  % Calculate h_lamdba(k_ij) matrix and normalize
   for l=1:m
     HH(:,l) = Hvector(Ranking(:,l));
   end

  % Update rule
   cNeurons = HH ./ (sum(HH,2) * ones(1,m));

  % Determine quantization error
   QuantizationErrors(i) = determine_qerror_linear( distNeuronsToData );
   
end

% Determine distance between final neurons and data
distNeuronsToData = determine_distance_relational_neurons_to_data( Distances, cNeurons );
   

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
  
  
function QuantizationError = determine_qerror_linear( distNeuronsToData )

%
% Determines the linear quantization error for median neurons
%
% A.Hasenfuss (c) 2007-2008
%

[WinnerDist, WinnerIdx] = min( distNeuronsToData, [], 1 );
QuantizationError = 0.5 * sum( WinnerDist );
