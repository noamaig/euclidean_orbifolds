rng(1);
addpath(genpath(['toolbox_graph']));
addpath(genpath(['yalmip']));
addpath(genpath(['MeshLibrary']));
load_name='bimba_cvd';
flattener=loadDataToFlattener([load_name '.mat']);
% flattener.flatten();
% load results\15-05-21_1339_bimba_cvd\bimba_cvd.mat
[V2,T2] = subdivision(flattener.M_orig.V, flattener.M_orig.F, 'middle');
% [V2,T2] = subdivision(V2, T2, 'middle');
M_orig=[];
M_orig.V=V2;
M_orig.F=T2;
% flattener=Flattener(M_orig,    flattener.inds([3  2 1])',[6 2]);
flattener=Flattener(M_orig,    flattener.inds, flattener.singularities);
flattener.flatten();
%% visualize

figure(200);
clf
max_e=1.5;%0.96
% flattener.correctGlobalAffine
flattener.visualize('confdist',2,[1 1.5]);
d=flattener.smax./flattener.smin;
caxis([1 max_e]);
figure();
flattener.drawHist(3);

%% mkdir
dirname=mkdir_for_results(load_name);
%%
fname=[dirname '/' load_name];
flattener.save([fname]);
S1=flattener.getX3D([load_name],'confdist',[max_e,1]);
S=x3dScene(S1);
f=fopen([fname '.x3d'],'w');
fprintf(f,S);
fclose(f);
