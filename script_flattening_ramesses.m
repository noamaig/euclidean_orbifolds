rng(1);
addpath(genpath(['toolbox_graph']));
addpath(genpath(['yalmip']));
addpath(genpath(['MeshLibrary']));
load_name='ramesses';
f_ours=loadDataToFlattener([load_name '.mat']);
f_ours.flatten();%'freesquare',false);


%% visualize

figure(200);
clf

f_ours.visualize('dirichlet');
colormap(jet);
caxis([0 2]);
return


%% mkdir
dirname=mkdir_for_results(load_name);
%%
fname=[dirname '/' load_name];
f_ours.save([fname]);
S1=f_ours.getX3D([load_name],'frobenius',[2,0]);
S=x3dScene(S1);
f=fopen([fname '.x3d'],'w');
fprintf(f,S);
fclose(f);
