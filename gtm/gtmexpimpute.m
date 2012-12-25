function [data_rec] = gtmexpimpute(net, data)
%GTMEXPIMPUTE Impute missing data using GTM expected values
%
%	Description
%	 DATA_REC = GTMEXPIMPUTE(NET, DATA) takes a GTM structure NET, and
%	imputes the missing values in DATA according to expectation of hidden
% variables given observed data.
%
%	See also
%	GTM, GTMEM, GTMLMEAN, GMLMODE, GMMPROB

% Copyright (c) Tommi Vatanen (2012)


% Check for consistency
errstring = consist(net, 'gtm', data);
if ~isempty(errstring)
  error(errstring);
end

data_rec = data;
missing = isnan(data_rec);
 
net.gmmnet.centres = rbffwd(net.rbfnet, net.X);
R = gmmpost(net.gmmnet, data);

rec = R*net.gmmnet.centres;
data_rec(missing) = rec(missing);