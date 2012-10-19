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

%%%%%%%%%%%%% DETAILED DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% gtm_make
%
% PURPOSE
%
% Creates, initializes and trains a GTM using default parameters.
%
% SYNTAX
%
%  gMap = gtm_make(D);
%  gMap = gtm_make(...,'argID',value,...);
%  gMap = gtm_make(...,value,...);
%
% DESCRIPTION
%
% Creates, initializes and trains a GTM with default parameters. Uses functions
% SOM_TOPOL_STRUCT, SOM_TRAIN_STRUCT, SOM_DATA_STRUCT and SOM_MAP_STRUCT to come
% up with the default values.
%
% First, the number of map units is determined. Unless they are
% explicitly defined, function SOM_TOPOL_STRUCT is used to determine this.
% It uses a heuristic formula of 'munits = 5*dlen^0.54321'. The 'mapsize'
% argument influences the final number of map units: a 'big' map has 
% x4 the default number of map units and a 'small' map has x0.25 the
% default number of map units. 
%
% After the number of map units has been determined, the map size is 
% determined. Basically, the two biggest eigenvalues of the training
% data are calculated and the ratio between sidelengths of the map grid
% is set to this ratio. The actual sidelengths are then set so that 
% their product is as close to the desired number of map units as
% possible.
%
% Then the GTM is initialized. First, linear initialization along two
% greatest eigenvectors is tried, but if this can't be done (the
% eigenvectors cannot be calculated), random initialization is used
% instead.  After initialization, the SOM is trained in two phases:
% first rough training and then fine-tuning. If the 'tracking'
% argument is greater than zero, the average quantization error and
% topographic error of the final map are calculated.
%
% REQUIRED INPUT ARGUMENTS
%
%  D           The data to use in the training.
%     (struct) A data struct. If a struct is given, '.comp_names' field as 
%              well as '.comp_norm' field is copied to the map struct.
%     (matrix) A data matrix, size dlen x dim. The data matrix may
%              contain unknown values, indicated by NaNs. 
%  
% OPTIONAL INPUT ARGUMENTS 
%
%  argID (string) Argument identifier string (see below).
%  value (varies) Value for the argument (see below).
%
% Here are the valid argument IDs and corresponding values. The values 
% which are unambiguous (marked with '*') can be given without the
% preceeding argID.
%   'init'       *(string) initialization: 'randinit' or 'lininit' (default)
%   'munits'      (scalar) the preferred number of map units
%   'msize'       (vector) map grid size
%   'mapsize'    *(string) do you want a 'small', 'normal' or 'big' map
%                          Any explicit settings of munits or msize override this.
%   'topol'      *(struct) topology struct
%   'comp_names'  (string array / cellstr) component names, size dim x 1
%   'tracking'    (scalar) how much to report, default = 1
%
% OUTPUT ARGUMENTS
% 
%  gMap (struct) the trained map struct
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
%  som_map_struct   Create a map struct.
%  som_topol_struct Default values for SOM topology.
%  som_train_struct Default values for SOM training parameters.
%  som_randinint    Random initialization algorithm.

% Copyright (c) 1999-2000 by the SOM toolbox programming team.
% http://

% Version 2.0beta juuso 111199

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
      [net, ~] = gtminit_imputation(net, initopt, D, 'units', ...
        sMapInit.codebook,'regular',latentShapeGTM,rbfgrid);  
    else
      net = gtminit(net, initopt, D, 'regular', latentShapeGTM, rbfgrid);
    end
  case 'sominit'
    sMapInit = som_map_struct(dim,sTopol); 
    sMapInit = som_lininit(D, sMapInit); 
    sTrain = som_train_struct(sMapInit,'dlen',dlen,'algorithm','batch','phase','rough');
    if missingData
      if tracking>0, fprintf(1,'initializing with imputation SOM\n'); end
      sMapInit = som_impbatch(sMapInit,D,sTrain,'tracking',tracking);
    else
      sMapInit = som_batchtrain(sMapInit,D,sTrain,'tracking',tracking);
    end
    [net, ~] = gtminit_imputation(net, initopt, D, 'units', ...
      sMapInit.codebook,'regular',latentShapeGTM,rbfgrid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% training

switch algorithm
  case 'sequential'
    if tracking>0, fprintf(1,'Training using sequential EM algorithm...\n'); end
    [net,~,errlog] = gtmem_new(net, D, emopt);
  case 'batch'
    if missingData
      if tracking>0, fprintf(1,'Training using batch EM algorithm for sparse data...\n'); end
      [net,~,errlog] = gtmem_imputation(net, D, emopt);
    else
      if tracking>0, fprintf(1,'Training using batch EM algorithm...\n'); end        
      [net,~,errlog] = gtmem(net, D, emopt);
    end
end

net.msize = latentShapeGTM;
net.errlog = errlog;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%