function [net, options, errlog] = gtmemseq(net, t, options)
%GTMEM	Sequential EM algorithm for Generative Topographic Mapping.
%
%	Description
%	[NET, OPTIONS, ERRLOG] = GTMEM(NET, T, OPTIONS) uses sequential Expectation
%	Maximization algorithm to estimate the parameters of a GTM defined by
%	a data structure NET. The matrix T represents the data whose
%	expectation is maximized, with each row corresponding to a vector.
%	It is assumed that the latent data NET.X has been set following a
%	call to GTMINIT, for example.    The optional parameters have the
%	following interpretations.
%
%	OPTIONS(1) is set to 1 to display error values; also logs error
%	values in the return argument ERRLOG. If OPTIONS(1) is set to 0, then
%	only warning messages are displayed.  If OPTIONS(1) is -1, then
%	nothing is displayed.
%
%	OPTIONS(3) is a measure of the absolute precision required of the
%	error function at the solution. If the change in log likelihood
%	between two steps of the EM algorithm is less than this value, then
%	the function terminates.
%
% OPTIONS(13) is the size of the batch in sequential training; default 10.
%
%	OPTIONS(14) is the maximum number of iterations; default 100.
%
%	The optional return value OPTIONS contains the final error value
%	(i.e. data log likelihood) in OPTIONS(8).
%
%	See also
%	GTM, GTMINIT
%
% [1] Tommi Vatanen. Missing value imputation using subspace methods with 
% applications on survey data. Master's thesis, Aalto University, 
% Espoo, Finland, 2012.

%	Copyright (c) Ian T Nabney (1996-2001)
% Modified by Tommi Vatanen (2012)

% Check that inputs are consistent
errstring = consist(net, 'gtm', t);
if ~isempty(errstring)
  error(errstring);
end

% Sort out the options
if (options(14))
  niters = options(14);
else
  niters = 100;
end

if (options(13))
  blen = options(13);
else
  blen = 10;
end

display = options(1);
store = 0;
if (nargout > 2)
  store = 1;	% Store the error values to return them
  errlog = zeros(1, niters);
end
test = 0;
if options(3) > 0.0
  test = 1;	% Test log likelihood for termination
end

% Calculate various quantities that remain constant during training
[ndata, tdim] = size(t);
ND = ndata*tdim;
[net.gmmnet.centres, Phi] = rbffwd(net.rbfnet, net.X);
Phi = [Phi ones(size(net.X, 1), 1)];
PhiT = Phi';
[K, Mplus1] = size(Phi);

A = zeros(Mplus1, Mplus1);
cholDcmp = zeros(Mplus1, Mplus1);
% Use a sparse representation for the weight regularizing matrix.
if (net.rbfnet.alpha > 0)
  Alpha = net.rbfnet.alpha*speye(Mplus1);
  Alpha(Mplus1, Mplus1) = 0;
end 

%variance = net.gmmnet.covars(1,1);
[R, act, d] = gtmpost(net, t);
Rt = R'*t;
G = spdiags(sum(R)', 0, K, K);

for n = 1:niters

  if (display || store || test)
    [~, act, d] = gtmpost(net, t);
    prob = act*(net.gmmnet.priors)';
    % Error value is negative log likelihood of data
    e = - sum(log(max(prob,realmin)));
    if store
      errlog(n) = e;
    end
    if display > 0
      fprintf(1, 'Cycle %4d  Error %11.6f\n', n, e);
    end
    if test
      if (n > 1 && abs(eold - e) < options(3))
        options(8) = e;
        errlog = errlog(errlog~=0);
        return;
      else
        eold = e;
      end
    end
  end
  
  % Sequential training - learning is done a block of data (inds) at a time
  % rather than in a single sweep to save memory consumption. (see [1])
  R_old = R;  
  i0 = 0;     
  while i0+1<=ndata
    inds = (i0+1):min(ndata,i0+blen); 
    i0 = i0+blen;      
    
    % Calculate responsibilities
    [R_new, act_new] = gtmpost(net, t(inds,:));
    % Calculate error value if needed
    R(inds,:) = R_new;
    act(inds,:) = act_new;
    
    % update estimate of G and R*X
    % there might be a faster way to do this(?)
    if length(inds) > 1
      G = G + spdiags(sum(R(inds,:) - R_old(inds,:))', 0, K, K);
    else
      G = G + spdiags((R(inds,:) - R_old(inds,:))', 0, K, K);
    end
    
    Rt = Rt + (R(inds,:) - R_old(inds,:))'*t(inds,:);
    
    % Calculate matrix be inverted (Phi'*G*Phi + alpha*I in the papers).
    % Sparse representation of G normally executes faster and saves memory
    if (net.rbfnet.alpha > 0)
      A = full(PhiT*G*Phi + (Alpha.*net.gmmnet.covars(1)));
    else
      A = full(PhiT*G*Phi);
    end
    % A is a symmetric matrix likely to be positive definite, so try
    % fast Cholesky decomposition to calculate W, otherwise use SVD.
    % (PhiT*(R*t)) is computed right-to-left, as R
    % and t are normally (much) larger than PhiT.
    [cholDcmp singular] = chol(A);
    if (singular)
      if (display)
        fprintf(1, ...
          'gtmem: Warning -- M-Step matrix singular, using pinv.\n');
      end
      W = pinv(A)*(PhiT*(Rt));
    else
      W = cholDcmp \ (cholDcmp' \ (PhiT*(Rt)));
    end
    % Put new weights into network to calculate responsibilities
    net.rbfnet.w2 = W(1:net.rbfnet.nhidden, :);
    net.rbfnet.b2 = W(net.rbfnet.nhidden+1, :);
    % Calculate new distances
    
    d(inds,:) = dist2(t(inds,:), Phi*W);
    
    % Calculate new value for beta
    variance = sum(sum(d.*R))/ND;
    net.gmmnet.covars = ones(1, net.gmmnet.ncentres)*variance;
  end
end

options(8) = -sum(log(gtmprob(net, t)));
%if (display >= 0)
%  disp(maxitmess);
%end
