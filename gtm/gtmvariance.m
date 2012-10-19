function variance = gtmvariance(x, c, ND, var_old, R)
%GMMVARIANCE Calculates variance for gtm with missing values.

% x == t, c == Phi*W

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
%   for i = 1:ndata
%     for j = 1:ncentres
%       for k = 1:dimx
%         if ~isnan(x(i,k))
%           n2(i,j) = n2(i,j) + (x(i,k)-c(j,k))^2;
%         end
%       end
%     end
%   end
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