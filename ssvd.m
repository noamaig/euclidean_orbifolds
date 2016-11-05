
function [ U,E,V,flip] = ssvd( A )
%ssvd implementation - return U,E,V so that U*E*V'=A, U,V are rotations, E
%is diagonal with E(end,end)=sign(det(A)), and lastly flip =sign(det(A))

%TODO: this can be improved to use only one determinant computation, of A
%itself as the det of U and V is irrelevant

%regular svd
[U,E,V]=svd(A);
%compute determinant of V
dV=det(V);
%get source and target dimensions
SD=size(A,2);
TD=size(A,1);
flip=0;
if dV<0 %if V is flipped, simply fix the determinant of both. This still keeps A=UEV'
    V(:,SD)=-V(:,SD);
    U(:,SD)=-U(:,SD);
end
%now det(V)>0, let's see if det(U)<0 - if so det(A)<0 
dU=det(U);
if dU<0
    if TD==SD %if we are of full dimension, only then can there be a flip
        flip=1;
    end
    U(:,TD)=-U(:,TD);%now fix det(U)>0
    E(TD,SD)=-E(TD,SD);%and correct E so that we keep A=UEV'
end
end

