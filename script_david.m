rng(1);
addpath(genpath(['toolbox_graph']));
addpath(genpath(['yalmip']));
addpath(genpath(['MeshLibrary']));
oldinds=[5809
    14550
    10027];
inds=[428 24764 5476];
% inds=[428 24764 5476 1000];
singularities=[4 4 ];
% inds=[428 24764 5476 7000];
% singularities=[2 2 2];
% M_orig=TriangleMesh('off','david_50k.off');
[V,T]=read_off('david_50k.off');
load_name='david';
V=V/mean(std(V'));
% V=bsxfun(@minus,V,mean(V));
M_orig=[];
M_orig.V=V;
M_orig.F=T;
flattener=Flattener(M_orig,inds,singularities);
flattener.flatten();

flattener.visualize();
return
% flattener.saveX3D('test',[],[max(flattener.frobenius),0]);
%%
t=Tiler(flattener.flat_V,flattener.flat_T,flattener.M_cut.pathPairs);
t.tile(2);
orgV=flattener.flat_V;

figure(317);

clf;
hold on
% t.trans{end+1}=[1 0;0 1; 0 0];
vcol=flattener.computeVCol();
for i=1:length(t.trans)
    if ~ismember(i,[1 4 5 8])
        continue
    end
    A=t.trans{i};
    curV=orgV*A([1 2],:)'+repmat(A(3,:),length(flattener.flat_V),1);
%     visualizeNew( curV,flattener.flat_T,flattener.M_cut.pathPairs,zeros(length(flattener.flat_T),1),vcol);
    flattener.flat_V=curV;
    flattener.visualize();
%     pause;
    
end

flattener.flat_V=orgV;
%%
hold off
lim=5.5;
axis([-lim+2 lim -lim lim-2]);
axis off;

figure(318);
clf
% V3=flattener.flat_V;
% V3=V3/max(max(V3)-min(V3));
% visualizeNew(V3,flattener.flat_T,flattener.M_cut.pathPairs,zeros(length(flattener.flat_T),1),flattener.vcol);
flattener.visualize
axis off
return
%% mkdir
dirname=mkdir_for_results(load_name);
%%
fname=[dirname '/' load_name];
flattener.save([fname]);
S1=flattener.getX3D([load_name],'vcol',[2,0]);
S=x3dScene(S1);
f=fopen([fname '.x3d'],'w');
fprintf(f,S);
fclose(f);
