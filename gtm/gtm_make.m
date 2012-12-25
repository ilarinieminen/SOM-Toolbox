function net = gtm_make(D, varargin)

%GTM_MAKE Create, initialize and train Generative Topographic Mapping.
%
% gMap = gtm_make(D, [[argID,] value, ...])
%
%  gMap = gtm_make(D);
%  gMap = gtm_make(D, 'munits', 20);
%  gMap = gtm_make(D, 'msize', [4 6], 'rbfgrid', [3 3]);
%
%  Input and output arguments ([]'s are optional): 
%   D        (matrix) training data, size dlen x dim
%            (struct) data struct
%   [argID,  (string) See below. The values which are unambiguous can 
%    value]  (varies) be given without the preceeding argID.
%
%   gMap     (struct) GTM struct
%
% Here are the valid argument IDs and corresponding values. The values 
% which are unambiguous (marked with '*') can be given without the
% preceeding argID.
%   'init'       *(string) initialization: 'lininit' or 'sominit' (default)
%   'algorithm'  *(string) 'batch' or 'sequential' (default) training
%   'blen'        (scalar) size of one batch in sequential training,
%                          default = 10
%   'rbfgrid'     (vector) size 1 x 2, size of the RBF grid in the nonlinear
%                          mapping
%   'regul'       (scalar) regularization parameter alpha, 
%                          default = 0, i.e., no regularization
%   'munits'      (scalar) the preferred number of map units
%   'msize'       (vector) map grid size
%   'mapsize'    *(string) do you want a 'small', 'normal' or 'big' map
%                          Any explicit settings of munits or msize override this.
%   'topol'      *(struct) topology struct
%   'comp_names'  (string array / cellstr) component names, size dim x 1
%   'tracking'    (scalar) how much to report, default = 1
%   'stopcrit'    (scalar) threshold to stop EM algorithm, default = 1e-2
%   'maxiters'    (scalar) max iterations in EM algorithm, default = 1000
%
% For more help, try 'type gtm_make' or check out online documentation.
% See also SOM_MAKE, SOM_MAP_STRUCT, SOM_TOPOL_STRUCT, SOM_LININIT
%
% EXAMPLES
%
%  To simply train a map with default parameters: 
%
%   gMap = gtm_make(D); 
%  
%  With the optional arguments, the initialization and training can be
%  influenced. To change map size, use 'msize', 'munits' or 'mapsize'
%  arguments:  
%
%   gMap = gtm_make(D,'mapsize','big'); or gMap=gtm_make(D,'big');
%   gMap = gtm_make(D,'munits', 100);
%   gMap = gtm_make(D,'msize', [20 10]); 
%
%  The 'tracking' argument can be used to control the amout of reporting
%  during training. The argument is used in this function, and it is
%  passed to the training functions. To make the function work silently
%  set it to 0.
%
%   gMap = gtm_make(D,'tracking',0); 
%
% SEE ALSO
%  
%  som_make         Create, initialize and train a SOM.
%
% Copyright (c) 1999-2012 by the SOM toolbox programming team.
% 
%
% Version 2.1 by Tommi Vatanen

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% check arguments

% D
if isstruct(D) 
  data_name = D.name; 
  comp_names = D.comp_names;
  comp_norm = D.comp_norm;
  D = D.data;
else 
  data_name = inputname(1);
  sDummy = som_data_struct(D(1,:)); 
  comp_names = sDummy.comp_names;
  comp_norm = sDummy.comp_norm;
end
[dlen dim] = size(D);
missingData = any(any(isnan(D)));

% defaults

emopt = foptions;   % options for EM algorithm
emopt(1) = 0;       % verbosity level
emopt(3) = 1e-2;    % stopping criterion
emopt(13) = 10;     % blen in sequential training
emopt(14) = 100;    % max iters

initopt = foptions; % options for gtm initialization
initopt(7) = 1;     % width factor of RBFs

mapsize = '';
alpha = 0;
sM = som_map_struct(dim, 'hexa'); 
sTopol = sM.topol;
munits = prod(sTopol.msize); % should be zero
tracking = 1;
rbfgrid = 0;
initalg = 'sominit';
algorithm = 'sequential';

% varargin
i=1; 
while i<=length(varargin), 
  argok = 1; 
  if ischar(varargin{i}), 
    switch varargin{i}, 
      % argument IDs
      case 'rbfgrid',    i=i+1; rbfgrid = varargin{i};
        numRbfCenters = prod(rbfgrid);
      case 'algorithm',   i=i+1; algorithm = varargin{i};
      case 'blen',        i=i+1; emopt(13) = varargin{i};
      case 'munits',     i=i+1; munits = varargin{i};
      case 'regul',      i=i+1; alpha = varargin{i};
      case 'msize',      i=i+1; sTopol.msize = varargin{i};
        munits = prod(sTopol.msize);
      case 'mapsize',    i=i+1; mapsize = varargin{i};
      case {'topol','som_topol','sTopol'},
        i=i+1; sTopol = varargin{i}; munits = prod(sTopol.msize);
      case 'tracking',   i=i+1; tracking = varargin{i};
      case 'init',       i=i+1; initalg = varargin{i};
      case 'stopcrit',   i=i+1; emopt(3) = varargin{i};
      case 'maxiters',   i=i+1; emopt(14) = varargin{i};
        % unambiguous values
      case {'small','normal','big'}, mapsize = varargin{i};
      case {'sominit','lininit'}, initalg = varargin{i};
      case {'sequential','batch'}, algorithm = varargin{i};
      otherwise argok=0;
    end
  elseif isstruct(varargin{i}) && isfield(varargin{i},'type'), 
    switch varargin{i}(1).type, 
     case 'som_topol', sTopol = varargin{i}; 
     otherwise argok=0; 
    end
  else
    argok = 0; 
  end
  if ~argok, 
    disp(['(gtm_make) Ignoring invalid argument #' num2str(i+1)]); 
  end
  i = i+1; 
end

emopt(1) = tracking-1;       % verbosity level

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make the map struct

%% map size
if isempty(sTopol.msize) || ~prod(sTopol.msize), 
  if tracking>0, fprintf(1,'Determining map size...\n'); end
  if ~munits,     
    sTemp = som_topol_struct('dlen',dlen);
    munits = prod(sTemp.msize);
    switch mapsize,
     case 'small', munits = max(9,ceil(munits/4));
     case 'big',   munits = munits*4;
     otherwise % nil
    end
  end
  sTemp = som_topol_struct('data',D,'munits',munits);
  sTopol.msize = sTemp.msize;
  munits = prod(sTemp.msize);
  if tracking>0, 
    fprintf(1,' map size [%d, %d]\n',sTopol.msize(1), sTopol.msize(2));   
  end
end
latentShapeGTM = fliplr(sTopol.msize);

% Heuristic selection of RBF grid. Experimental selection is encouraged.
if rbfgrid == 0
  rbfgrid = [ceil(log(max(sTopol.msize))) ceil(log(max(sTopol.msize)))];
  numRbfCenters = prod(rbfgrid);
  if tracking>0, 
    fprintf(' rbf grid [%d, %d]\n',rbfgrid(1), rbfgrid(2));   
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initialization

net = gtm(2, munits, dim, numRbfCenters,'gaussian', alpha(1));

if tracking>0, fprintf(1,'Initialization...\n'); end

switch initalg, 
  case 'lininit'
    if missingData
      sMapInit = som_map_struct(dim,sTopol);
      sMapInit = som_lininit(D, sMapInit);
      [net, ~] = gtminit2(net, initopt, D, sMapInit.codebook, 'regular' ...
        ,latentShapeGTM,rbfgrid);  
    else
      net = gtminit(net, initopt, D, 'regular', latentShapeGTM, rbfgrid);
    end
  case 'sominit'
    sMapInit = som_map_struct(dim,sTopol); 
    sMapInit = som_lininit(D, sMapInit); 
    sTrain = som_train_struct(sMapInit,'dlen',dlen,'algorithm','imp','phase','rough');
    if missingData
      if tracking>0, fprintf(1,'initializing with Imputation SOM\n'); end
      sMapInit = som_impbatch(sMapInit,D,sTrain,'tracking',tracking);
    else
      sMapInit = som_batchtrain(sMapInit,D,sTrain,'tracking',tracking);
    end
    [net, ~] = gtminit2(net, initopt, D, sMapInit.codebook, 'regular', ...
      latentShapeGTM, rbfgrid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% training

switch algorithm
  case 'sequential'
    if tracking>0, fprintf(1,'Training using sequential EM algorithm...\n'); end
    [net,~,errlog] = gtmemseq(net, D, emopt);
  case 'batch'
    if missingData
      if tracking>0, fprintf(1,'Training using batch EM algorithm for sparse data...\n'); end
      [net,~,errlog] = gtmem2(net, D, emopt);
    else
      if tracking>0, fprintf(1,'Training using batch EM algorithm...\n'); end        
      [net,~,errlog] = gtmem(net, D, emopt);
    end
end

net.msize = latentShapeGTM;
net.errlog = errlog;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%