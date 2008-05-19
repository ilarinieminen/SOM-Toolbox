function p=vis_patch(lattice)

% VIS_PATCH Defines the basic patches (hexa and rect) used in SOM_CPLANE
%
%  p = vis_patch(lattice)
%
%  Input and output arguments: 
%   lattice   (string) 'rect', 'hexa' or 'hexagon'
%
%   p         (matrix) size Lx2, defines the vertices of the patch
%
% This function produces vertex coordinates for a patch presenting
% a map unit in hexagonal or rectangular lattice with its centre in (0,0). 
%
% For more help, try 'type vis_patch' or check out online documentation.
% See also SOM_CPLANE, PATCH.

%%%% DETAILED DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% vis_patch
%
% SYNTAX
%  
%  p = vis_patch(lattice)
% 
% DESCRIPTION
% 
% Forms a map unit patch for SOM_CPLANE function. Mainly a subroutine
% of SOM_CPLANE, although can be used for on its own as well.
%
% REQUIRED INPUT ARGUMENTS
%
% lattice (string)
% 
%    'hexa'  produces vertex coordiantes for a hexagoanl patch which 
%            has its center on (0,0), unit width, and a height of
%            1.3334 units. This is not a regular hexagon but such that
%            the node which has topological coordinates (a,b) has its
%            center in the visualization at coordinates (a,b) for odd
%            a and at (a,b+.5) for even a. The non-regular look of the
%            patch is taken care simply by changing the axis ratio.
%
%    'rect'  produces vertex coordinates for a uniform rectangular patch.
%            having its center on (0,0) and unit sidelength. Used as a
%            subfunction in SOM_CPLANE.
%
%    'hexagon' produces vertex coordinates for a regular hexagonal patch.
%            It may be used in som_cplane if, for some reason, truly
%            regular hexagons are needed instead of the default unit
%            markers which are not uniform, but have integer
%            y-coordinates in the lattice.
%
% OUTPUT ARGUMENTS
%
% p (matrix) The 2-dimensional vertex coordinates: 
%
%   case 'rect'        case 'hexa'                case 'hexagon'
%    p=[[-.5 -.5];...     p=[[0     0.6667];...    p=[[0     0.5774];...
%       [-.5  .5];...        [0.5   0.3333];...       [0.5   0.2887];...
%       [ .5  .5];...        [0.5  -0.3333];...       [0.5  -0.2887];...
%       [ .5 -.5]];          [0    -0.6667];...       [0    -0.5774];...
%                            [-0.5 -0.3333];...       [-0.5 -0.2887];...
%                            [-0.5  0.3333]];         [-0.5  0.2887]]; 
%
% EXAMPLES
%
%  som_cplane(vis_patch('rect'),[6 5],'none');
%  % this produces the same result as som_cplane('rect',[6 5], 'none') 
%  
%  som_cplane(vis_patch('hexa'), vis_unit_coord('hexa',[6 5]), 'none');
%  % produces in principle the same result as 
%  % som_cplane(vis_patch('hexa'),[6 5],'none'), 
%  % _but_ in this case the axis are not rescaled and the non-regular 
%  % shape of hexagons can be seen.
%
%  som_cplane(vis_patch('hexagon'), som_unit_coords([6 5],'hexa'), 'none');
%  % produces a truly regular hexa lattice 
%
% SEE ALSO 
%
%  vis_unit_coord   The default 'hexa' and 'rect' coordinates in visualizations
%  som_unit_coords  Locations of units on the SOM grid.

% Copyright (c) 1999-2000 by the SOM toolbox programming team.
% http://www.cis.hut.fi/projects/somtoolbox/             

% Version 2.0beta Johan 041099

if ~ischar(lattice)
  error('Input argument should be a string')
else
  switch lattice
  case 'rect'
    p=[[-.5 -.5]; ...
      [-.5 .5];...
      [.5 .5];...
      [.5 -.5]];
  case 'hexagon'
    p=[[0 0.5774];...
      [0.5 0.2887];...
      [0.5 -0.2887];...
      [0 -0.5774];...
      [-0.5 -0.2887];...
      [-0.5 0.2887]];
  case 'hexa'
    p=[[0 0.6667];...
      [0.5 0.3333];...
      [0.5 -0.3333];...
      [0 -0.6667];...
      [-0.5 -0.3333];...
      [-0.5 0.3333]];
  otherwise
    error('Unknown lattice');
  end
end


