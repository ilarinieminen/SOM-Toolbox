
function [net, data] = gtminit2(net, options, data, init, samp_type, varargin)
%GTMINIT Initialise the weights and latent sample in a GTM given reference vectors.
%
%	Description
%	NET = GTMINIT(NET, OPTIONS, DATA, INIT, SAMPTYPE) takes a GTM NET and
%	generates a sample of latent data points and sets the centres (and
%	widths if appropriate) of NET.RBFNET.
%
% INIT defines initial locations of reference vectors (Gaussian means).
%
%	If the SAMPTYPE is 'REGULAR', then regular grids of latent data
%	points and RBF centres are created.  The dimension of the latent data
%	space must be 1 or 2.  For one-dimensional latent space, the
%	LSAMPSIZE parameter gives the number of latent points and the
%	RBFSAMPSIZE parameter gives the number of RBF centres.  For a two-
%	dimensional latent space, these parameters must be vectors of length
%	2 with the number of points in each of the x and y directions to
%	create a rectangular grid.  The widths of the RBF basis functions are
%	set by a call to RBFSETFW passing OPTIONS(7) as the scaling
%	parameter.
%
%	If the SAMPTYPE is 'UNIFORM' or 'GAUSSIAN' then the latent data is
%	found by sampling from a uniform or Gaussian distribution
%	correspondingly.  The RBF basis function parameters are set by a call
%	to RBFSETBF with the DATA parameter as dataset and the OPTIONS
%	vector.
%
%	Finally, the output layer weights of the RBF are initialised by
%	mapping the mean of the latent variable to the mean of the target
%	variable, and the L-dimensional latent variale variance to the
%	variance of the targets along the first L principal components.
%
%	See also
%	GTMINIT, GTM, GTMEM, PCA, RBFSETBF, RBFSETFW
%

%	Copyright (c) Ian T Nabney (1996-2001)
% Modified by Tommi Vatanen (2012)

% Check for consistency
errstring = consist(net, 'gtm', data);
if ~isempty(errstring)
  error(errstring);
end

% Check type of sample
stypes = {'regular', 'uniform', 'gaussian', 'special'};
if (strcmp(samp_type, stypes)) == 0
  error('Undefined sample type.')
end

if net.dim_latent > size(data, 2)
  error('Latent space dimension must not be greater than data dimension')
end
nlatent = net.gmmnet.ncentres;
nhidden = net.rbfnet.nhidden;

% Create latent data sample and set RBF centres

switch samp_type
case 'regular'
   if nargin ~= 7
      error('Regular type must specify latent and RBF shapes');
   end
   l_samp_size = varargin{1};
   rbf_samp_size = varargin{2};
   if round(l_samp_size) ~= l_samp_size
      error('Latent sample specification must contain integers')
   end
   % Check existence and size of rbf specification
   if any(size(rbf_samp_size) ~= [1 net.dim_latent]) || ...
         prod(rbf_samp_size) ~= nhidden
      error('Incorrect specification of RBF centres')
   end
   % Check dimension and type of latent data specification
   if any(size(l_samp_size) ~= [1 net.dim_latent]) || ...
         prod(l_samp_size) ~= nlatent
      error('Incorrect dimension of latent sample spec.')
   end
   if net.dim_latent == 1
      net.X = [-1:2/(l_samp_size-1):1]';
      net.rbfnet.c = [-1:2/(rbf_samp_size-1):1]';
      net.rbfnet = rbfsetfw(net.rbfnet, options(7));
   elseif net.dim_latent == 2
      net.X = gtm_rctg(l_samp_size);
      net.rbfnet.c = gtm_rctg(rbf_samp_size);
      net.rbfnet = rbfsetfw(net.rbfnet, options(7));
   else
      error('For regular sample, input dimension must be 1 or 2.')
   end
   
case {'uniform', 'gaussian'}
   if strcmp(samp_type, 'uniform')
      net.X = 2 * (rand(nlatent, net.dim_latent) - 0.5);
   else
      % Sample from N(0, 0.25) distribution to ensure most latent 
      % data is inside square
      net.X = randn(nlatent, net.dim_latent)/2;
   end   
   net.rbfnet = rbfsetbf(net.rbfnet, options, net.X);
   
%   case {'special'}
%     net.X = [-0.5 sqrt(3)/2;0.5 sqrt(3)/2;0 0];
%     net.rbfnet.c = [-0.5 sqrt(3)/2;0.5 sqrt(3)/2;0 0];
%     net.rbfnet = rbfsetfw(net.rbfnet, options(7));
otherwise
   % Shouldn't get here
   error('Invalid sample type');
end

% use given unit vectors as the initial reference vectors
[~, Phi] = rbffwd(net.rbfnet, net.X);
% initparam are the reference vectors
net.rbfnet.w2 = Phi \ init; 
net.rbfnet.b2 = zeros(1,size(data,2));

net.gmmnet.centres = rbffwd(net.rbfnet, net.X);
d = dist2(net.gmmnet.centres, net.gmmnet.centres) + ...
  diag(ones(net.gmmnet.ncentres, 1)*realmax);
sigma = mean(min(d))/2;
    
net.gmmnet.covars = sigma*ones(1, net.gmmnet.ncentres);


% Sub-function to create the sample data in 2d
function sample = gtm_rctg(samp_size)

xDim = samp_size(1);
yDim = samp_size(2);
% Produce a grid with the right number of rows and columns
[X, Y] = meshgrid([0:1:(xDim-1)], [(yDim-1):-1:0]);

% Change grid representation
sample = [X(:), Y(:)];

% Shift grid to correct position and scale it
maxXY= max(sample);
sample(:,1) = 2*(sample(:,1) - maxXY(1)/2)./maxXY(1);
sample(:,2) = 2*(sample(:,2) - maxXY(2)/2)./maxXY(2);

return;