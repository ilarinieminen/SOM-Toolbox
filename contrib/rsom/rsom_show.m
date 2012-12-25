function rsom_show(sMap, D, varargin)

%RSOM_SHOW visualizes information about a given RSOM, currently only the 
%U-matrix.
%
% rsom_show(sMap, D);
%
% Input and output arguments: 
%   sMap     (struct) RSOM map struct
%   D        (matrix) Dissimilarity matrix used for training of the RSOM
%
% For more help, try 'type rsom_lininit' or check out online documentation.
% See also RSOM_UMAT.

%%%%%%%%%%%%% DETAILED DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% rsom_show
%
% PURPOSE
%
% Visualizes information concerning the given RSOM, currently only the 
% U-matrix.
%
% SYNTAX
%
%  sM = rsom_show(sMap, D);
%
% DESCRIPTION
%
% Visualizes the U-matrix defined by the SOM.
%
% REFERENCES
%
%   Barbara Hammer, Alexander Hasenfuss: Topographic Mapping of Large
%   Dissimilarity Data Sets. Neural Computation 22(9): 2229-2284 (2010)
%
% REQUIRED INPUT ARGUMENTS
%
%   sMap     (struct) RSOM map struct
%   D        (matrix) Dissimilarity matrix used for training of the RSOM
%
% EXAMPLES
%
%   sM = rsom_show(sMap, D);
%
% SEE ALSO
%
%   rsom_umat    Calculates the U-matrix

% Contributed to SOM Toolbox vs2, December 7th, 2012 by Alexander Schulz
% Copyright (c) Alexander Schulz
% http://www.cis.hut.fi/projects/somtoolbox/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h_axes = subplot(1,1,1);

lattice = sMap.topol.lattice;
msize   = sMap.topol.msize;

u = rsom_umat(sMap, D);
u = u(:);

tmp_h=som_cplane([lattice 'U'],msize,u);
set(tmp_h,'EdgeColor','none');
set(h_axes,'Tag','Uplane');
h_label=xlabel('U-matrix');
set(h_label,'Visible','on','verticalalignment','top');
set(gca,'plotboxaspectratio',[msize(2) msize(1) msize(1)]);


set(h_label,'interpreter','none'); 
colormap(flipud(gray(64).^2));
colorbar;