function T = getTensorTranspose(n,m)
% Returns a matrix T such that vec(X')=T*vec(X) for any nxm matrix X.
%
% Example:
% n = 3;
% m = 4;
% X = zeros(n,m);
% X(:)=1:numel(X);
% T = getTensorTranspose(n,m);
% vec(X')
% T*vec(X)

[i,j] = meshgrid(1:n,1:m);
T = sparse(sub2ind([m n],j,i),sub2ind([n m],i,j),1);

