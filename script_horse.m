addpath(genpath(['toolbox_graph']));
addpath(genpath(['yalmip']));
addpath(genpath(['MeshLibrary']));
inds=[1
       1000
       100
       2000];
   singularities=[2 2 2];
M_orig=TriangleMesh('off','horse_5k.off');
[V,T,pathPairs,vcol]=flatten(M_orig,inds,singularities);
%%
t=Tiler(V,T,pathPairs);
t.tile(3);
figure(317);

clf;
% t.trans{end+1}=[1 0;0 1; 0 0];
for i=1:length(t.trans)
    A=t.trans{i};
    curV=V*A([1 2],:)'+repmat(A(3,:),length(V),1);
visualizeNew( curV,T,pathPairs,zeros(length(T),1),vcol);
end

% axis([-3 3 -3 3]);