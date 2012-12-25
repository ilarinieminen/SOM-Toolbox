function variance = gtmvariance(x, c, ND, var_old, R)

% GMMVARIANCE Calculates variance for GTM with missing values.
%
% See also
% GTMEM2, GTM

% Copyright (c) Tommi Vatanen 2012


[ndata, dimx] = size(x);
[ncentres, dimc] = size(c);
if dimx ~= dimc
	error('Data dimension does not match dimension of centres')
end

if any(any(isnan(x)))
  n2 = zeros(ndata,ncentres);
  for i = 1:ndata
    ind = ~isnan(x(i,:));
    n2(i,:) = sum((repmat(x(i,ind), [ncentres 1]) - c(:,ind)).^2,2)' ...
      + sum(isnan(x(i,:)))*var_old;
  end
else
  n2 = (ones(ncentres, 1) * sum((x.^2)', 1))' + ...
    ones(ndata, 1) * sum((c.^2)',1) - ...
    2.*(x*(c'));
end

% Rounding errors occasionally cause negative entries in n2
if any(any(n2<0))
  n2(n2<0) = 0;
end

variance = (sum(sum(n2.*R))/ND);