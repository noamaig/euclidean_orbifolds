rng(1);
addpath(genpath(['toolbox_graph']));
addpath(genpath(['yalmip']));
addpath(genpath(['MeshLibrary']));
load_name='max';
flattener=loadDataToFlattener([load_name '.mat']);
flattener.flatten();


%% visualize

figure(200);
clf
max_e=1.5;%0.96

flattener.visualize('ncol',2);
d=flattener.smax./flattener.smin;
caxis([1 max_e]);
figure();
hist(d,100);


%% mkdir
dirname=mkdir_for_results(load_name);
%%
fname=[dirname '/' load_name];
flattener.save([fname]);
S1=flattener.getX3D([load_name],'vcol',[max_e,1]);
S=x3dScene(S1);
f=fopen([fname '.x3d'],'w');
fprintf(f,S);
fclose(f);
