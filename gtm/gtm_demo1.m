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


% load wineData and wineLabels
load winedata

% Normalize the data, use SOM toolbox functionality
sD = som_normalize(wineData, 'var');

% Use wrapper function gtm_make to train a GTM
net = gtm_make(sD);

gtm_show(net, 'mags', 'data', sD, 'groups', wineLabels)

% set(gcf,'units','centimeters');
% set(gcf,'DefaultAxesFontSize',8)
% set(gcf,'DefaultTextFontSize',8)
% set(gcf,'pos',[13.6 5 8 12])
