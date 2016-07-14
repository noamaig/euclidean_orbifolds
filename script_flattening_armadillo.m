rng(1);
addpath(genpath(['toolbox_graph']));
addpath(genpath(['yalmip']));
addpath(genpath(['MeshLibrary']));
load_name='Armadillo';
f_ours=loadDataToFlattener([load_name '.mat']);
f_ours.flatten();%'freesquare',false);

vorig=f_ours.flat_V;
%% mkdir
dirname=mkdir_for_results(load_name);

%% visualize
f_ours.flat_V=vorig;
f_ours.computeDistortion(true);
figure(200);
clf
dist_lim=1.5;
f_ours.visualize('confdist',2,[1 dist_lim]);
%%
fname=[dirname '/' load_name];
f_ours.save([fname '_before']);
S1=f_ours.getX3D([load_name '_before'],'confdist',[dist_lim,1]);
S=x3dScene(S1);
f=fopen([fname '_before.x3d'],'w');
fprintf(f,S);
fclose(f);

%%
figure(201);

clf
f_ours.correctGlobalAffine();
f_ours.visualize('confdist',2,[1 dist_lim])
%%
fname=[dirname '/' load_name];
f_ours.save([fname '_after']);
S1=f_ours.getX3D([load_name '_after'],'confdist',[dist_lim,1]);
S=x3dScene(S1);
f=fopen([fname '_after.x3d'],'w');
fprintf(f,S);
fclose(f);
return
