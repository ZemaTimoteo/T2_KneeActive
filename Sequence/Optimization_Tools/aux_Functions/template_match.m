function [ ind_X ] = template_match( S,X )
%[ ind_X ] = template_match( S,X )
%   size(X) = [nTEs num_param]; % imagem em 1-by-mxn
%   size(S) = [nTEs Nsq]; % dicionario em nTEsxNsq-by-1

% Modifications AF 15.05.2020
%%

%normalize the matrices before the calculation of their inner product
%X = normc(X);
%S = normc(abs(S)); 

%normalize the matrices
for i=1:size(X(:,:),2)
    X(:,i)=X(:,i)./norm(X(:,i));
end
S = abs(S);
for i=1:size(S(:,:),2)
 
    S(:,i)=S(:,i)./norm(S(:,i));
end

%size(inner_product) = [num_param Nsq];
inner_product = X'*S;

%find the index with the highest inner product for each pixel
[~,ind_X] = max(abs(inner_product));

end

