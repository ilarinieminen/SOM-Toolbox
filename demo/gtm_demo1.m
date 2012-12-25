% GTM_DEMO1 Basic properties and behaviour of the GTM functionality in 
% SOM toolbox.
%
% Contributed to SOM Toolbox 2.1, Oct 19th, 2012 by Tommi Vatanen
%
% This demonstation uses wine data from UCI machine learning repository
% (http://archive.ics.uci.edu/ml/datasets/Wine). It is a simple data with
% wines from three different wine regions. Data contains 13 variables for
% each wine (wineData) and class labels (wineLabels) that is the region the
% wine is from.
%
% Copyright (c) Tommi Vatanen 2012


% load wineData and wineLabels
load winedata

% Normalize the data, use SOM toolbox functionality
sD = som_normalize(wineData, 'var');

% Cross-validate regularization parameter alpha
cv = cvpartition(size(sD,1), 'kfold', 10);
alphas = logspace(2,-2,10);
logprobs = zeros(10,numel(alphas));

for i = 1:numel(alphas)
  for j = 1:1
    net = gtm_make(sD(cv.training(j),:), 'rbfgrid', [6 6], 'regul', alphas(i));
    logprobs(j,i) = -sum(log(gtmprob(net, sD(cv.test(j),:))));
  end
end

%% plot cross-validation results
semilogx(alphas,mean(logprobs),'bo-');

% the plot shows that alpha = 2-3 is en elbow point, thus suitable as
% regularization parameter% GTM_DEMO1 Basic properties and behaviour of the GTM functionality in 
% SOM toolbox.
%
% Contributed to SOM Toolbox 2.1, Oct 19th, 2012 by Tommi Vatanen
%
% This demonstation uses wine data from UCI machine learning repository
% (http://archive.ics.uci.edu/ml/datasets/Wine). It is a simple data with
% wines from three different wine regions. Data contains 13 variables for
% each wine (wineData) and class labels (wineLabels) that is the region the
% wine is from.
%
% Copyright (c) Tommi Vatanen 2012


% load wineData and wineLabels
load winedata

% Normalize the data, use SOM toolbox functionality
sD = som_normalize(wineData, 'var');

% Cross-validate regularization parameter alpha
cv = cvpartition(size(sD,1), 'kfold', 10);
alphas = logspace(2,-2,10);
logprobs = zeros(10,numel(alphas));

for i = 1:numel(alphas)
  for j = 1:1
    net = gtm_make(sD(cv.training(j),:), 'rbfgrid', [5 5], 'regul', alphas(i));
    logprobs(j,i) = -sum(log(gtmprob(net, sD(cv.test(j),:))));
  end
end

%% plot cross-validation results
semilogx(alphas,mean(logprobs),'bo-');

% the plot shows that alpha = XXX is en elbow point, thus suitable as
% regularization parameter

alpha = XXX

%% Use wrapper function gtm_make to train a GTM
net = gtm_make(sD, 'rbfgrid', [5 5],'regul', alpha);
gtm_show(net, 'mags', 'data', sD, 'groups', wineLabels);

