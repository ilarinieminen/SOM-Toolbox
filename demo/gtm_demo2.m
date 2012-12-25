% GTM_DEMO2 Basic properties and behaviour of the GTM functionality in 
% SOM toolbox.
%
% Contributed to SOM Toolbox 2.1, Dec 16, 2012 by Tommi Vatanen
%
% This demonstation uses ISOLET data from UCI machine learning repository
% (http://archive.ics.uci.edu/ml/datasets/ISOLET). 
%
% Copyright (c) Tommi Vatanen 2012


% load isoletData and isoletLabels
load isoletdata

% Normalize the data, use SOM toolbox functionality
sD = som_normalize(isoletData, 'var');

%% Cross-validate regularization parameter alpha
cv = cvpartition(size(sD,1), 'kfold', 10);
alphas = logspace(-1,-5,20);
logprobs = zeros(10,numel(alphas));

for i = 1:numel(alphas)
  for j = 1:1
    net = gtm_make(sD, 'rbfgrid', [6 6], 'regul', alphas(i));
    logprobs(j,i) = -sum(log(gtmprob(net, sD)));
  end
end

%% plot cross-validation results
semilogx(alphas,mean(logprobs),'bo-');

% the plot shows that alpha = 2-3 is en elbow point, thus suitable as
% regularization parameter

alpha = 18e-5

%% Use wrapper function gtm_make to train a GTM
net = gtm_make(sD, 'rbfgrid', [5 5],'regul', alpha);

gtm_show(net, 'mags', 'data', sD, 'groups', wineLabels);