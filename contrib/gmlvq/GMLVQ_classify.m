function estimatedLabels = GMLVQ_classify(Data, model)

%GMLVQ_classify.m - classifies the given data with the given model
%
% estimatedLabels = GMLVQ_classify(Data, model)
%
% Input and output arguments: 
%  trainSet          (struct)
%                    (matrix) matrix with training samples in its rows
%  model    : GMLVQ model with prototypes w their labels c_w and the matrix omega
% 
%  estimatedLabels       : the estimated labels
%
%  example for usage:
%  trainSet = [1,2,3;4,5,6;7,8,9];
%  trainLab = [1;1;2];
%  GMLVQ_model=gmlvq_core(trainSet,trainLab); % minimal parameters required
%  estimatedTrainLabels = GMLVQ_classify(trainSet, GMLVQ_model);
%  trainError = mean( trainLab ~= estimatedTrainLabels );
%  
% Kerstin Bunte
% kerstin.bunte@googlemail.com
% Mon Nov 05 09:05:52 CEST 2012
%
% Conditions of GNU General Public License, version 2 and BSD License apply.
% See file 'license-gpl2.txt' and 'BSD_license.txt' enclosed in this package.
% Programs are not for use in critical applications!
%

% Contributed to SOM Toolbox vs2, December 3rd, 2012 by Alexander Schulz
% Copyright (c) Kerstin Bunte
% http://www.cis.hut.fi/projects/somtoolbox/

if isstruct(Data)
    Data = Data.data;
end

if ~isfield(model, 'w')
    model.c_w = som_label2num(model.labels);
    dist = computeDistance(Data, model.codebook, model);
else
    dist = computeDistance(Data, model.w, model);
end

[~, index] = min(dist,[],2);

estimatedLabels = model.c_w(index);



function distance = computeDistance(X, W, model)
nb_samples = size(X,1);
distance = zeros(nb_samples,length(model.c_w));


% tic;
% for i = 1:size(W,1)
%     delta = X - ones(P,1) * W(i,:);
%     delta = bsxfun(@minus, X, W(i,:));
%     delta(isnan(delta)) = 0;
    % Hadamard product: to skip unnecessary calculation between two different examples
%     distance(1:nb_samples,i) = sum( ((delta*model.omega'*model.omega).*delta) ,2 );   
%     wi = W(i,:);
%     distance(1:nb_samples,i) = sum( (( (X-wi(ones(nb_samples,1),:)) *model.omega'*model.omega).*(X-wi(ones(nb_samples,1),:))) ,2 );
% end
% disp(toc);
if isfield(model,'psis')
    if length(model.psis)~=size(W,1)
        classes = unique(model.c_w);
        for i = 1:size(W,1)
            matrixIdx = classes==model.c_w(i);
%             delta = X-ones(nb_samples,1)*W(i,:);            
%             distance(1:nb_samples,i) = sum( ((delta* model.psis{matrixIdx}'*model.psis{matrixIdx}).*delta) ,2 );
            distance(1:nb_samples,i) = sum((bsxfun(@minus, X, W(i,:))*model.psis{matrixIdx}').^2, 2);
        end
    else
        for i = 1:size(W,1)
%             delta = X-ones(nb_samples,1)*W(i,:);
%             distance(1:nb_samples,i) = sum( ((delta* model.psis{i}'*model.psis{i}).*delta) ,2 );
            distance(1:nb_samples,i) = sum((bsxfun(@minus, X, W(i,:))*model.psis{i}').^2, 2);
        end
    end
else
    if isfield(model,'lambda')
        for i = 1:size(W,1)
%             delta = X-ones(nb_samples,1)*W(i,:);
%             distance(1:nb_samples,i) = sum( (( delta *model.lambda).*delta) ,2 );
            delta = bsxfun(@minus, X, W(i,:));
            distance(1:nb_samples,i) = sum(bsxfun(@times,delta.^2,model.lambda), 2);
%             distance(1:nb_samples,i) = sum( (( delta *model.lambda).*delta) ,2 );
        end
    end
    if isfield(model,'omega')
        % tic;
        for i = 1:size(W,1)
%             delta = X-ones(nb_samples,1)*W(i,:);
%             distance(1:nb_samples,i) = sum( (( delta *model.omega'*model.omega).*delta) ,2 );
            distance(1:nb_samples,i) = sum((bsxfun(@minus, X, W(i,:))*model.omega').^2, 2);
        end
        % disp(toc);
    end
end
