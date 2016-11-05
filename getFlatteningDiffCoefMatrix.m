function [T,areas] = getFlatteningDiffCoefMatrix(V,F)
%Input: V,F vertices and tris of a matrix
%output: T, a matrix, which once applied to vector W of vertex positions, 
% returns the linear trasnformations of each triangle taking it from V to W
% coordinates. Namely, T*W where W is a (x_1,y_1,x_2,y_2,...) vector, will
% return % a vector of the form (a_1,b_1,c_1,d_1,a_2,b_2,c_2... etc.),
% s.t. each quadruplet matches a matrix, A_i=[a_i,b_i;c_i d_i], 
% which is the linear transforatmion induced by W, of the i'th triangle

% [T,areas]=computeMeshTranformationCoeffsMex(F,V);


% return matrix such that:
% vec(Ai)=Acoeff*V(:) 

% init
tol = 1e2*eps;
dim = size(V,2);
n_vert = size(V,1);
n_faces = size(F,1);
n_simplex = size(F,2);
assert(dim==3 & n_simplex==3, 'Can only be applied to 3d *surfaces*');
B = eye(n_simplex) - 1/(n_simplex); % centering matrix
diffDim = dim-1; % dimension of jacobian (square)

% prepare indexing helpers
[I_mask, J_mask] = ndgrid(1:diffDim^2,1:n_simplex);
J_add_mask = reshape(bsxfun(@plus,zeros(diffDim,n_simplex),permute((0:diffDim-1)*n_vert,[3 1 2])),[],diffDim^2)';

% calculate jacobian coefficients
I = cell(n_faces,1);
J = cell(n_faces,1);
S = cell(n_faces,1);
for ii = 1:n_faces
    % get current triangle
    currF = F(ii,:);
    currV = V(currF,:);
    %     % trasform to plane
    % RFlat = compute2DEmbedding(currV')';
    % currVFlat = currV*RFlat;
    currVFlat = embedTriangle(currV);

    % calculate differential
    currAcoeff = (B*currVFlat)\B;
    %currAcoeff(abs(currAcoeff)<tol) = 0;
    % calculate indices in full tranformation matrix
    I{ii} = I_mask + (ii-1)*diffDim^2;
    J{ii} = currF(J_mask) + J_add_mask;
    S{ii} = repmat(currAcoeff, [diffDim,1]);
end

% gather
I = cat(1,I{:});
J = cat(1,J{:});
S = cat(1,S{:});
T = sparse(I,J,S,diffDim^2*n_faces,diffDim*n_vert);
areas=ones(length(F),1);
return
V=V';
F=F';
% prepare stuff
nVert = size(V,2);
nFaces = size(F,2);
B = eye(3) - 1/3; % centering matrix

% calculate flattenning coefficients
I = cell(nFaces,1);
J = cell(nFaces,1);
S = cell(nFaces,1);
counter = 1;
allFlatX=nan(length(F),3);
allFlatY=nan(length(F),3);
for ii = 1:nFaces
    % get current triangle
    currF = F(:,ii);
    currV = V(:,currF);
    % trasform to plane
    RFlat = find_2d_embedding(currV);
    currVFlat = RFlat*currV;
    allFlatX(ii,:)=currVFlat(1,:);
    allFlatY(ii,:)=currVFlat(2,:);
    % calculate differential
    currT = B/(currVFlat*B);
    currT(abs(currT)<=1e2*eps) = 0;
    % calculate indices in full tranformation matrix
    I{ii} = [counter counter counter; counter+1 counter+1 counter+1];
    J{ii} = [currF'; currF'];
    S{ii} = currT';
    counter = counter + 2;
end
areas=polyarea(allFlatX',allFlatY')';
% gather
I = cat(1,I{:});
J = cat(1,J{:});
S = cat(1,S{:});
% T = sparse(I,J,S,2*nFaces,nVert);
IX=I*2-1;
JX=J*2-1;
IY=I*2;
JY=J*2;
T = sparse([IX; IY],[JX ;JY],[S;S],2*2*nFaces,2*nVert);


% v = randn(size(V,2),2);
% 
% 
% Aold = permute(reshape(T*v,2,[],2),[1 3 2]);
% Anew = reshape(Tnew*v(:),2,2,[]);
% 
% for ii = 1:size(Aold,3)
%     assert(norm(Anew(:,:,ii)'*Anew(:,:,ii)-Aold(:,:,ii)'*Aold(:,:,ii),'fro')/norm(Anew(:,:,ii),'fro')^2<1e-10)
% end




end
function flatV = embedTriangle(V)

v1 = V(2,:)-V(1,:);
v2 = V(3,:)-V(1,:);

norm_v1 = norm(v1);
norm_v2 = norm(v2);
cos_theta = dot(v1,v2)/(norm_v1*norm_v2);
sin_theta = sqrt(1-cos_theta^2);


flatV = [0 0;...
    norm_v1 0;...
    norm_v2*cos_theta norm_v2*sin_theta];
end