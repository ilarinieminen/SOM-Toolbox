function [nargin_,varargout]=vis_planeGetArgs(varargin)

% VIS_PLANEGETARGS Subfunction for som_*plane: extracts topolopy 
%                  information from the first arguments.
%
% [nargin,varargout]=vis_planeGetArgs(varargin)
%
%  Input and output arguments: 
%   varargin   (varies) arguments given to som_*plane function
%   nargin_    (scalar) number of arguments that nargchk of som_*plane "should see"
%                       +number_of_varargins if varargin{1} is not map/topol struct
%                       +number_of_varargins+1 if varargin{2} is a map/topol struct
%   varargout  (varies) the arguments that som_*plane "should see"
%
% Basically, this function allows topology information to be given 
% in various ways: either as a map/topology struct, or as a argument pair:
% lattice, msize. The topology is always converted into the (lattice, msize)
% argument pair.
%  - if first input argument (varargin{1}) is a map or topol struct 
%    the function extracts lattice and msize fields to two first 
%    output variables after 'nargin_'. 
%  - otherwise it copies the input arguments to the output arguments 
%    after 'nargin_'. 
% If there are too many inputs (as compared to number of outputs), the 
% last ones are ignored. If too few, they are replaced by empty values 
% in outputs.
%
% Example of usage: 
%   function definition: h = som_cplane(varargin)
%   first code line:     [nargin,lattice,msize,color,size,pos]=vis_planeGetArgs(varargin);
%
% See also SOM_CPLANE, SOM_BARPLANE, SOM_PLOTPLANE, SOM_PIEPLANE.

% Copyright (c) 2000 by the SOM toolbox programming team.
% http://www.cis.hut.fi/projects/somtoolbox/             

% Version 2.0beta Johan 240300

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nout=nargout-1;

% Set first all varargins to contain empty (==default values in som_*plane)

for i=1:nout, varargout{i}=[]; end

nargin_ = nargin;
% Struct: might be map or topol

if isstruct(varargin{1}),

  % Get topol from topol field
  if isfield(varargin{1},'topol'), topol=varargin{1}.topol;
  else topol=varargin{1}; % assume that this is topol struct 
  end

  if ~isstruct(topol),
    % topol not a struct !?
    warning('Field ''topol'' is not a struct.');
    varargout{1}=varargin{1};
    varargoutC=2;
    nargin_ = nargin;
  elseif ~isfield(topol,'msize') | ~isfield(topol,'lattice'),
    % Field missing?!
    warning('Invalid topology struct.');
    varargout{1}=topol;
    varargoutC=2;
    nargin_ = nargin;
  else
    varargout{1}=topol.lattice;
    varargout{2}=topol.msize;
    % increment input arg. counter
    varargoutC=3;
    nargin_ = nargin+1;    
  end

elseif iscell(varargin{1}), 

  c = varargin{1}; 
  lattice = 'hexa'; shape = 'sheet'; msize = [1 1]; 
  for i=1:length(c), 
    if ischar(c{i}), 
      switch c{i}, 
      case {'hexa','hexaU','rect','rectU'}, lattice = c{i}; 
      case {'sheet','cyl','toroid'}, shape = c{i}; 
      end
    else
      msize = c{i}; 
    end 
  end
  varargout{1} = lattice;
  varargout{2} = msize;
  varargoutC=3;
  nargin_ = nargin+1;    

else

  % should be a lattice (string) 
  varargout{1}=varargin{1};
  varargoutC=2;
  nargin_=nargin;

end

for i=2:nargin, varargout{varargoutC+i-2}=varargin{i}; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
