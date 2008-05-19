function h=vis_footnote(txt)

% VIS_FOOTNOTE Adds a movable text to the current figure
%
%  h = vis_footnote(T)
%
%  Input and output arguments ([]'s are optional)
%   [T]  (string) text to be written
%        (scalar) font size to use in all strings 
%
%   h    (vector) handles to axis objects created by this function 
%
% This function sets a text to the current figure. If T is a string,
% it's written as it is to the same place. If T is a scalar, the font
% size of all text objects created by this function are changed to the
% pointsize T. If no input argument is given the function only returns
% the handles to all objects created by this function. The texts may
% be dragged to a new location at any time using mouse.  Note that the
% current axis will be the parent of the text object after dragging.
%
% String 'Info' is set to the Tag property field of the objects. 
% 
% EXAMPLES
%
% % add movable texts to the current figure and change their
% % fontsize to 20 points
% vis_footnote('Faa'); vis_footnote('Foo'); vis_footnote(20);
% 
% % delete all objects created by this function from the current figure
% delete(vis_footnote);
% 
% See also SOM_SHOW.

% Copyright (c) 1997-2000 by the SOM toolbox programming team.
% http://www.cis.hut.fi/projects/somtoolbox/             

% Version 2.0beta Johan 080698

%% Check arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

error(nargchk(0, 1, nargin))  % check no. of input args

%% Init %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the handles to the existing Info-axes objects

h_infotxt=findobj(gcf,'tag','Info','type','text');
h_infoax=findobj(gcf,'tag','Info','type','axes');

%% Action %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If no arguments are given, return the old axes handles

if nargin == 0 | isempty(txt),
  ;  
elseif ischar(txt)                    % text: set new text
  [t,h_]=movetext(txt);
  h_infoax=[h_; h_infoax];
elseif vis_valuetype(txt,{'1x1'})      % scalar: change font size  
  set(h_infotxt,'fontunits','points');
  set(h_infotxt,'fontsize',txt);
else
  error('Input argument should be a string or a scalar.');
end

%% Build output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargout>0     % output only if necessary
  h=h_infoax;
end

%%% SUBFUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

function [t,h]=movetext(txt)
% Moves the text. See also VIS_FOOTNOTEBUTTONDOWNFCN
%
%
initpos=[0.05 0.05 0.01 0.01];   

memaxes = gca;                   % Memorize the gca

%% Create new axis on the lower left corner.
%% This will be the parent for the text object

h = axes('position',initpos,'units','normalized');
set(h,'visible','off');          % hide axis

t = text(0,0,txt);               % write text 
set(t,'tag','Info');             % set tag
set(h,'tag','Info');             % set tag

set(t,'verticalalignment','bottom');  % set text alignment
set(t,'horizontalalignment','left');

%% Set ButtonDownFcn

set(t,'buttondownfcn','vis_footnoteButtonDownFcn') 

axes(memaxes);                   % Reset original gca


