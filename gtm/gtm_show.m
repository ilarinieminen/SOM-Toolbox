function h=gtm_show(net, varargin)

% GTM_SHOW Basic GTM visualizations
%
% h = som_show(sMap, ['argID', value, ...])
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
%  h           (struct) struct with the following fields:
%   .plane     (vector) handles to the axes objecets (subplots)
%   .colorbar  (vector) handles to the colorbars. Colorbars for empty
%                       grids & RGB color planes do not exist: the
%                       value for them in the vector is -1.
%   .label     (vector) handles to the axis labels
%
% Here are the valid argument IDs and corresponding values. M is
% the number of map units
%  'comp'               Which component planes to draw, title is
%                       the name of the component (from sMap.comp_names) 
%              (vector) a vector of component indices
%              (string) 'all' (or '' or []) for all components
%  'compi'              as 'comp' but uses interpolated shading
%  'umat'               Show u-matrix calculated using specified 
%                       components 
%              (vector) a vector of component indeces
%              (string) 'all' (or '' or []) to use all components
%              (cell)   of form {v, str} uses v as the vector, and put
%                       str as title instead of the default 'U-matrix'
%  'umati'              as 'umat' but uses interpolated shading of colors 
%  'empty'     (string) Make an empty plane using given string as title
%  'color'              Set arbitrary unit colors explicitly  
%              (matrix) size Mx1 or Mx3, Mx1 matrix uses indexed
%                       coloring;  Mx3 matrix (RGB triples as rows)
%                       defines fixed unit colors
%              (cell)   of from {color, str}. 'color' is the Mx1
%                       or Mx3 RGB triple matrix and 'str' is title 
%                       string
%  'colori'             as 'color' but uses interpolated shading of colors 
%  'norm'      (string) 'n' or 'd': Whether to show normalized 'n' or 
%                       denormalized 'd' data values on the
%                       colorbar. By default denormalized values are used.
%  'bar'       (string) Colorbar direction: 'horiz', 'vert' (default)
%                       or 'none'
%  'size'               size of the units
%              (scalar) same size for each unit, default is 1
%              (vector) size Mx1, individual size for each unit
%  'edge'      (string) Unit edges on component planes 'on'
%                       (default) or 'off'
%  'footnote'  (string) Footnote string, sMap.name by default
%  'colormap'  (matrix) user defined colormap 
%  'subplots'  (vector) size 1 x 2, the number of subplots in y- and
%                       and x-directions (as in SUBPLOT command)
%
% See also SOM_SHOW

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