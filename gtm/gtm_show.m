function h=gtm_show(net, varargin)

% GTM_SHOW Basic GTM visualization
%
% h = gtm_show(net, ['argID', value, ...])
% 
%  gtm_show(net);
%  gtm_show(net,'mags');
%  gtm_show(net,'mags','data',sD,'groups',labels);
%
% Input and output arguments ([]'s are optional):
%  net         (struct) map struct
%  [argID,     (string) Additional parameters are given as argID, value
%    value]    (varies) pairs. See below for list of valid IDs and values.
%
%  h           (handle) handle to the plot figure
%
% Here are the valid argument IDs and corresponding values. M is
% the number of map units
%  'mags'      Show magnification factors
%  'data'      Data to be plotted 
%  'groups'    Groups or labels of the plotted data
%                (plot each group with distinctive marker)
%  'mapping'   'mean' (default) or 'mode'
%              Plot expected values of hidden variables (mean)
%              or best matching units (mode)
%
% SEE ALSO
%  
%  gtm_make         Create, initialize and train a GTM.
%
% Copyright (c) 1999-2012 by the SOM toolbox programming team.
% 
% Version 2.1 by Tommi Vatanen

errstring = consist(net, 'gtm');
if ~isempty(errstring)
  error(errstring);
end

comp_planes = 0;
mags = 0;
data = 0;
groups = 0;
mapping = 'mean';

% varargin
i=1; 
while i<=length(varargin), 
  argok = 1; 
  if ischar(varargin{i}), 
    switch varargin{i}, 
      % argument IDs
      case 'comp', i=i+1; comp_planes = varargin{i};
      case 'mags', mags = 1;
      case 'data', i=i+1; data = varargin{i};
        errstring = consist(net, 'gtm', data);
        if ~isempty(errstring)
          error(errstring);
        end
      case 'groups', i=i+1; groups = varargin{i};
      case 'mapping', i=i+1; mapping = varargin{i};
    end
  else
    argok = 0; 
  end
  if ~argok, 
    disp(['(gtm_show) Ignoring invalid argument #' num2str(i+1)]); 
  end
  i = i+1; 
end

% component planes
if comp_planes~=0
  for i = comp_planes
    figure;
    
  end
end

h = figure;

% magnification factors
if mags
  mags = gtmmag(net, net.X);
  % Reshape into grid form
  mags = reshape(mags, fliplr(net.msize));
  imagesc(net.X(:, 1), net.X(:,2), mags);
  colormap(1-gray);
  hold on;
end

if data
  % map data to latent space
  switch mapping
    case 'mean'
      mappedData = gtmlmean(net, data);
    case 'mode'
      mappedData = gtmlmode(net, data);      
    otherwise
      error('Unknown mapping, can be mean or mode.')
  end
  
  % does user provide groups (labels)
  if groups
    gscatter(mappedData(:,1), mappedData(:,2), groups);
  else
    scatter(mappedData(:,1), mappedData(:,2));
  end
end

xbroaden = 1/(1*net.msize(1));
ybroaden = 1/(1*net.msize(2));
axis([-1-xbroaden 1+xbroaden -1-ybroaden 1+ybroaden]);
set(gca, 'xtick', [], 'ytick', []);
tmp = get(gcf, 'pos');
tmp(3) = tmp(4)*net.msize(1)/net.msize(2);
set(gcf, 'pos', tmp)