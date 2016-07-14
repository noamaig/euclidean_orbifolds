function cols  = distortion_colormap(N)
if nargin==0
    N=1000;
end
CBcols = [1 1 1];%[0.85 0.85 0.85];%[0.9 0.9 0.9];
t=(1:N)/N;t=t';
cols = (1-sqrt(t))*CBcols + sqrt(t)*[0.9 0.7 0.7];
cols(cols>1) =1;


CBcols = [0.9 0.9 0.9];
t=(1:64) /64;t=t';
cols = t*[1 0 0]+(1-t)*CBcols ;
cols(cols>1) =1;
%caxis([1 10]);

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% my color scheme for distortion
CBcols = [0.85 0.85 0.85];%[0.9 0.9 0.9];
t=(1:N) /N;t=t';
cols = (1-t)*CBcols + t*[1 0.0 0.0];
cols(cols>1) =1;
cols=jet(N)*0.4+0.6;

t=(1:N) /N;t=t';
cols1 = (1-t)*[0.7 0.78 1] + t*[1 1 1];
cols2 = (1-t)*[1 1 1] + t*[1 0.5 0.5];
cols=[cols1;cols2];
cols(cols>1) =1;



end

