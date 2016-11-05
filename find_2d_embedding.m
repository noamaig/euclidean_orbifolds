
function A = find_2d_embedding(V)
%given a 3X3 matrix whose columns are the vertices of a 3d triangle,
%returns the transformation A, s.t A*V gives embedding in 2D
v1=V(:,1)-V(:,2);
v2=V(:,1)-V(:,3);
%v3=cross(v1,v2);
A=orth([v1 v2])';
%make sure the embedding gives positive orinetation
if det(A*[v1 v2])<0
    A=A( [2,1],:);
end
% if nargin>1
%     if det([A;normal])<0
%         A=A(2:1,:);
%     end
% end

end

