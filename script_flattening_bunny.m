rng(1);
addpath(genpath(['toolbox_graph']));
addpath(genpath(['yalmip']));
addpath(genpath(['MeshLibrary']));
load_name='bunny_watertight';
f_ours=loadDataToFlattener([load_name '.mat']);
f_ours.flatten();
f_ours.computeDistortion();

%% visualize

figure(200);
clf

f_ours.visualize();
colormap(jet);
e_max=prctile(f_ours.frobenius,98);%    1.6474

caxis([0 e_max]);
axis off;


%% mkdir
dirname=mkdir_for_results(load_name);
%%
fname=[dirname '/' load_name];
f_ours.save([fname]);
S1=f_ours.getX3D([load_name],'vcol',[1,0]);
S=x3dScene(S1);
f=fopen([fname '.x3d'],'w');
fprintf(f,S);
fclose(f);
