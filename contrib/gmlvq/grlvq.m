function [sCodebook, settings] = grlvq(sCodebook, sData, varargin)

%GRLVQ trains a codebook using the Generalized Relevance LVQ (GRLVQ) 
%algorithm.
%
% sM = grlvq(sM, sD, [argID, value, ...])
%
%  sM = grlvq(sM, sD)
%  sM = grlvq(sM, sD, 'learningRatePrototypes', [0.01 0.001], 'relevanceStart', 10);
%
% Required Input Arguments:
%   sM                     (struct) map struct containing the labels as the 
%                                   first column of .labels
%   sD                     (struct) data struct containing the labels as 
%                                   the first column of .labels
%
% Optional Input Arguments:
%   PrototypesPerClass     (vector) number of prototypes per class used;
%                                   a number or a vector with a number 
%                                   for each class (default=1)
%   initialRelevances      (vector) initial Relevances 
%   regularization         (scalar) regularization for relevances
%   testSet                (matrix) a test set to compute the test error;
%                                   the last column should be the labels 
%   comparable             (scalar) a flag which resets the random 
%                                   generator to produce comparable results
%                                   if set to 1
%   optimization           (string) choice for the optimization technique: 
%                                   sgd or fminlbfgs (default=fminlbfgs)
%  Parameter for the stochastic gradient descent sgd:
%   nb_epochs              (scalar) the number of epochs (default=100)
%   learningRatePrototypes (vector) learning rate for the prototypes; could 
%                                   be the start and end value or a vector
%                                   of length nb_epochs
%   learningRateRelevances (vector) learning rate for the relevances; could 
%                                   be the start and end value or a vector 
%                                   of length nb_epochs
%   relevanceStart         (scalar) epoch to start the relevance learning
%                                   (default=1)
%  Parameter for the build-in function fminlbfgs:
%   threshstop             (scalar) the training error for early stopping
%                                   (default=0)
%   useEarlyStopping       (scalar) use early stopping based on threshstop
%                                   (default=1)
%   Display                (string) the optimization output 'iter' or 'off'
%                                   (default='off')
%   GradObj                (string) turn the usage of gradient information 
%                                   on or off (default='on')
%   HessUpdate             (string) the update can be 'lbfgs', 'bfgs' or 
%                                   'steepdesc' (default='lbfgs')
%   TolFun                 (scalar) the tolerance (default=1e-6)
%   MaxIter                (scalar) the maximal number of iterations 
%                                   (default=2500)
%   MaxFunEvals                     the maximal number of function
%                                   evaluations (default=1000000)
%   TolX                   (scalar) tolerance on the minimum 
%                                   (default=1e-10)
%   DiffMinChange          (scalar) minimal change (default=1e-10)
%
% Output Arguments:
%   sM                     (struct) map struct containing the optimized 
%                                   prototypes and the relevances lambda 
%   settings               (struct) information on the settings of the
%                                   algorithm
%
% NOTE: does not take the vector mask into account but trains a relevance
%       for each dimension.
% 
% For more help, try 'type grlvq', or check out the online documentation.
% See also GMLVQ, GRLVQ_CORE, LVQ3, LVQ1.

%%%%%%%%%%%%% DETAILED DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% grlvq
%
% PURPOSE
%
% Trains a codebook and a global metric with the Generalized Relevance LVQ
% (GRLVQ) algorithm. 
%
% SYNTAX
%
%   sM = grlvq(sM, sD)
%   sM = grlvq(sM, sD, 'learningRatePrototypes', [0.01 0.001], 'relevanceStart', 10);
%   [sM, settings] = grlvq(sM, sD, argID, value, ...);
%
% DESCRIPTION
%
% GRLVQ constitutes an enhancement of the LVQ algorithm by minimizing an
% explicit error function and optimizing the metric of the space where the
% computation is performed. Distances are calculated by 
%   d(x,y) = (x - y) * diag(lambda) * (x - y)'.
% The vector lambda is optimized by the algorithm.
%
% Two optimization techniques are implemented: 
%     - stochastic gradient descent (sgd)
%     - limited memory Quasi Newton Broyden-Fletcher-Goldfarb-Shanno
%       (L-BFGS)
%
% Classification is performed similarly to standard LVQ, with the only
% difference that for the computation of distances the above specified 
% formula is applied.
%
% This function provides an interface in the style of the som-toolbox and
% is basically a wrapper for the method grlvq_core. It deals with structs
% rather then directly with data matrices.
%
% REFERENCES
%
%   Barbara Hammer, Thomas Villmann: Generalized relevance learning vector 
%   quantization. Neural Networks 15(8-9): 1059-1068 (2002).
%
%   P. Schneider, K. Bunte, B. Hammer and M. Biehl: Regularization in 
%   Matrix Relevance Learning, IEEE Transactions on Neural Networks, vol. 
%   21, nb. 5, pp. 831-840, 2010.
%
% REQUIRED INPUT ARGUMENTS
%   sM                     (struct) map struct containing the initialized
%                                   prototype vectors as rows of the
%                                   .codebooks matrix and the labels as the
%                                   first column of .labels
%   sD                     (struct) data struct containing the data vectors
%                                   as rows of .data and the labels as the
%                                   first column of .labels
%
% OPTIONAL INPUT ARGUMENTS
%   PrototypesPerClass     (vector) number of prototypes per class used;
%                                   a number or a vector with a number 
%                                   for each class (default=1)
%   initialRelevances      (vector) initial Relevances
%                                   (default: vector of ones)
%   regularization         (scalar) regularization for relevances
%                                   (detault=0)
%   comparable             (scalar) a flag which resets the random 
%                                   generator to produce comparable results
%                                   if set to 1 (default=0)
%   optimization           (string) choice for the optimization technique: 
%                                   sgd (stochastic gradient descent) or 
%                                   fminlbfgs (Limited memory Quasi Newton 
%                                   Broyden-Fletcher-Goldfarb-Shanno)  
%                                   (default=fminlbfgs)
%  Parameter for the stochastic gradient descent sgd:
%   nb_epochs              (scalar) the number of epochs (default=100)
%   learningRatePrototypes (vector) learning rate for the prototypes; could 
%                                   be the start and end value or a vector
%                                   of length nb_epochs
%   learningRateRelevances (vector) learning rate for the relevances; could 
%                                   be the start and end value or a vector 
%                                   of length nb_epochs
%   relevanceStart         (scalar) epoch to start the relevance learning 
%                                   (default=1)
%  Parameter for the build-in function fminlbfgs:
%   threshstop             (scalar) the training error for early stopping
%                                   (default=0)
%   useEarlyStopping       (scalar) use early stopping based on threshstop
%                                   (default=1)
%   Display                (string) the optimization output 'iter' or 'off'
%                                   (default='off')
%   GradObj                (string) turn the usage of gradient information 
%                                   on or off (default='on')
%   HessUpdate             (string) the update can be 'lbfgs', 'bfgs' or 
%                                   'steepdesc' (default=lbfgs)
%   TolFun                 (scalar) the tolerance (default=1e-6)
%   MaxIter                (scalar) the maximal number of iterations 
%                                   (default=2500)
%   MaxFunEvals                     the maximal number of function
%                                   evaluations (default=1000000)
%   TolX                   (scalar) tolerance on the minimum 
%                                   (default=1e-10)
%   DiffMinChange          (scalar) minimal change (default=1e-10)
%
% OUTPUT ARGUMENTS
%   sM                     (struct) map struct containing the optimized 
%                                   prototypes as well as the relevances 
%                                   for each dimensionlambda 
%   settings               (struct) information on the settings of the
%                                   algorithm
%
% EXAMPLES
% 
%   lab   = unique(sD.labels(:,1));           % different classes
%   nProt = length(lab)*5;                    % 5 prototypes for each    
%   sM = som_randinit(sD,'msize',[nProt 1]);  % initial prototypes
%   sM.labels = [lab;lab;lab;lab;lab];        % and their classes
%   sM = grlvq(sM,sD);                        % use GRLVQ to adjust the
%                                             % prototypes and the metric
%
% SEE ALSO
%
%   gmlvq         Use the GMLVQ algorithm for training.
%   grlvq_core    Access the GRLVQ functionality without using structs.
%   lvq3          Use LVQ3 algorithm for training.
%   lvq1          Use LVQ1 algorithm for training.

% Contributed to SOM Toolbox vs2, December 3rd, 2012 by Alexander Schulz
% Copyright (c) Alexander Schulz
% http://www.cis.hut.fi/projects/somtoolbox/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



cods      = sCodebook.codebook;
codLabels = sCodebook.labels;
data      = sData.data;
labels    = sData.labels;

% convert labels from cell to a vector
labelsV    = som_label2num(labels);
codLabelsV = som_label2num(codLabels);

% calculate the number of prototypes per class
n_prot_per_class = zeros(1, max(codLabelsV));
for i=1:max(codLabelsV)
    n_prot_per_class(i) = sum(codLabelsV == i);
end


% call the main function
[model, settings]  = grlvq_core(data, labelsV, 'initialPrototypes', ... 
    [cods, codLabelsV], 'PrototypesPerClass', n_prot_per_class, varargin{:});



% write the results in the output struct
sCodebook.codebook = model.w;
sCodebook.lambda   = model.lambda;

if strcmp(settings.optimization, 'sgd')
    trainlen   = settings.nb_epochs;
    alpha_ini  = settings.learningRatePrototypes(1);
    alpha_type = 'power';
else
    trainlen   = NaN;
    alpha_ini  = NaN;
    alpha_type = '';
end



sTrain = som_set('som_train','algorithm','GRLVQ',...
		 'data_name',sData.name,...
		 'neigh','',...
		 'mask',ones(size(model.w,2),1),...
		 'radius_ini',NaN,...
		 'radius_fin',NaN,...
		 'alpha_ini',alpha_ini,... 
		 'alpha_type',alpha_type,...
		 'trainlen',trainlen,...
		 'time',datestr(now,0));
sCodebook.trainhist(end+1) = sTrain;


