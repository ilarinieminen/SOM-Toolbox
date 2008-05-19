function [handles,msg,lattice,msize,dim,normalization,comps]=vis_som_show_data(p,f)

% VIS_SOM_SHOW_DATA Checks and returns UserData and subplot handles stored by SOM_SHOW
%
% [handles,msg,lattice,msize,dim,normalization,comps] = vis_som_show_data(p, f)
%
%  Input and output arguments ([]'s are optional): 
%   [p]           (vector) subplot numbers 
%                 (string) 'all' to process all subplots, this is default
%                          'comp' to process only subplots which have
%                          component planes
%   [f]           (double) figure handle, default is current figure
%
%   handles       (vector) handles of requested subplots
%   msg           (string) error message or empty string (no error)
%   lattice       (string) map lattice: 'hexa' or 'rect'
%   msize         (vector) map grid size in figure
%   dim           (scalar) map data dimension in figure
%   normalization (struct) normalization struct used in the map in figure
%   comps         (vector) the component indexes in figure
%
% This function gets the handles of component planes and u-matrices in
% subplots p from figure f. SOM_SHOW writes the handles into the
% UserData field of the figure where their order won't be mixed
% up. This function reads the data according to the vector p. If the
% figure has been manipulated (original planes are missing) the function
% warns user or returns error string.
% 
% The main purpose for this is to be a subfuncion for SOM_SHOW_ADD,
% SOM_SHOW_CLEAR and SOM_RECOLORBAR functions, but it may be used on
% command line in the followong manner:
%
%  % plots text on the fifth plane
%  axes(vis_som_show_data(5)); hold on; text(1,3,'I am here');
%    
% See also SOM_SHOW, SOM_SHOW_ADD.

% Copyright (c) 1997-2000 by the SOM toolbox programming team.
% http://www.cis.hut.fi/projects/somtoolbox/             

% Version 2.0beta Johan 201099 juuso 160600

%% Check input args %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

error(nargchk(0, 2, nargin))  % check no. of input args 

%% Init %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

handles=[];                                     % initialize output
normalize=[];
comps=[];
dim=[];
msize=[];
lattice=[];
msg=[];      

cr=sprintf('\n');                               % carriage return            

if nargin < 2 | isempty(f)
  f=gcf;                                        % default figure
end

if nargin < 1 | isempty(p)                      % default p 
  p= 'all';
end

% Find component planes and u-matrices from the figure and get the 
% UserData field where the handles for the components are 
% in the original order. 
% If the fields are corrupted, return an error message.

h_real = [findobj(f, 'Tag', 'Cplane'); ...
    findobj(f, 'Tag', 'Uplane'); ...
    findobj(f,'Tag','CplaneI'); ...
    findobj(f,'Tag','UplaneI')];
eval( 'h_stored=getfield(get(f,''UserData''),''subplotorder'');' , ...
    'msg=[ msg cr '' Missing SOM_SHOW.subplotorder''];');  
eval( 'normalization=getfield(get(f,''UserData''),''comp_norm'');' , ...
    'msg=[msg cr '' Missing SOM_SHOW.comp_norm''];');
eval( 'comps=getfield(get(f,''UserData''),''comps'');' , ...
    'msg=[msg cr '' Missing SOM_SHOW.comps''];');    
eval( 'msize=getfield(get(f,''UserData''),''msize'');' , ...
    'msg=[msg cr '' Missing SOM_SHOW.msize''];');    
eval( 'dim=getfield(get(f,''UserData''),''dim'');' , ...
    'msg=[msg cr '' Missing SOM_SHOW.dim''];');    
eval( 'lattice=getfield(get(f,''UserData''),''lattice'');' , ...
    'msg=[msg cr '' Missing SOM_SHOW.lattice''];');    
if ~isempty(msg), 
  msg=['The figure does not contain SOM_SHOW visualization or is corrupted.'...
	cr msg cr cr ...
	'This command may be applied only to a SOM_SHOW visualization.']; 
  return; 
end

%% Check arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

index=ismember(h_stored, h_real);  % the original order for plot axes 

if ~prod(double(index))                    % missing planes?!
                                           % double added by kr 1.10.02
  l1= 'Some of the original planes seems to be missing.';
  l2= 'Subplot numbers now refer to the existing ones.';
  warning([l1 cr l2]);
end

if ~prod(double(ismember(h_real, h_stored))) % extra planes?! 
                                             % double added by kr 5.9.02
  warning('There seems to be new planes. Subplot numbers refer to the old ones.');
end

h_stored=h_stored(index);          % existing original plots in original order

if ischar(p)                       % check if p is 'all'
  switch(p)
   case 'all'                                   
    p=1:size(h_stored,1);          % all original subplots
   case 'comp'
    p=find(comps>0); 
   otherwise
    msg= 'String value for subplot number vector has to be ''all''!';
    return;
  end
end

if ~vis_valuetype(p,{ '1xn','nx1'}) % check the size
  msg= 'Subplot numbers (argument p in help text) have to be in a vector!';
  return
end

if min(p) < 1                      % check for invalid values
  msg= 'Subplot numbers (argument p in help text) must be at least 1!';
  return
end

%% p is too large

if max(p) > size(h_stored,1)
  l1= 'There are not so many existing subplots created by SOM_SHOW in the';
  l2= 'figure as you are trying to refer with subplot numbers.';
  l3= 'This is probably caused by having a too large subplot number.';
  l4= 'However, the reason may be invalid manipulation of';
  l5= 'this figure object or a program failure, too.';
  msg=([l1 cr l2 cr l3 cr cr l4 cr l5]);
  return;
end

%% Action and building output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

handles=h_stored(p);
comps=comps(p);

