function [data_rec] = gtmmapimputate(net, data)

%GTMEXPIMPUTE Impute missing data using GTM MAP estimates
%
%	Description
%	 DATA_REC = GTMEXPIMPUTE(NET, DATA) takes a GTM structure NET, and
%	imputes the missing values in DATA according to MAP estimates of hidden
% variables given observed data.
%
%	See also
%	GTMEXPIMPUTE, GTM, GTMEM, GTMLMEAN, GMLMODE, GMMPROB

% Copyright (c) Tommi Vatanen (2012)

% Check for consistency
errstring = consist(net, 'gtm', data);
if ~isempty(errstring)
  error(errstring);
end

data_rec = data;
missing = isnan(data_rec);

% 
net.gmmnet.centres = rbffwd(net.rbfnet, net.X);
R = gmmpost(net.gmmnet, data);

R_tmp = zeros(size(R));
R_max = repmat(max(R,[],2), [1 size(R,2)]);
R_tmp(R==R_max) = 1;

% pick up map reference vectors
rec = R_tmp*net.gmmnet.centres;
data_rec(missing) = rec(missing);
