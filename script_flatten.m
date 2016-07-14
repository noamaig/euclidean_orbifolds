rng(1);
addpath(genpath(['toolbox_graph']));
addpath(genpath(['yalmip']));
addpath(genpath(['MeshLibrary']));
load_name='eros';
f_ours=loadDataToFlattener([load_name '.mat']);
f_ours.flatten();


%% visualize

figure(200);
clf

f_ours.visualize('dirichlet');
colormap(jet);
caxis([0 15]);



%% mkdir
dirname=mkdir_for_results(load_name);
%%
fname=[dirname '/' load_name];
f_ours.save([fname '_ours']);
S1=f_ours.getX3D(['our_' load_name],[],[15,0]);
S=x3dScene(S1);
f=fopen([fname '_ours.x3d'],'w');
fprintf(f,S);
fclose(f);
