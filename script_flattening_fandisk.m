rng(1);
addpath(genpath(['toolbox_graph']));
addpath(genpath(['yalmip']));
addpath(genpath(['MeshLibrary']));
load_name='fandisk';
flattener=loadDataToFlattener([load_name '.mat']);
flattener.singularities=[3 3];
flattener.flatten();
% load bimba_cvd.mat

%% visualize

%%
[V2,T2] = subdivision(flattener.M_orig.V, flattener.M_orig.F, 'middle');
[V3,T3] = subdivision(V2, T2, 'middle');
[V4,T4] = subdivision(V3, T3, 'middle');
% [V4,T4] = subdivision(V3, T3, 'middle');
% [V4,T4] = subdivision(V4, T4, 'middle');
% M_orig=[];
% M_orig.V=V2;
% M_orig.F=T2;
% f2=Flattener(M_orig,flattener.inds,flattener.singularities);
% f2.flatten();
% 
% M_orig=[];
% M_orig.V=V3;
% M_orig.F=T3;
% f3=Flattener(M_orig,flattener.inds,flattener.singularities);
% f3.flatten();


M_orig=[];
M_orig.V=V4;
M_orig.F=T4;
f4=Flattener(M_orig,flattener.inds,flattener.singularities);
f4.flatten();
%%
figure(200);
% flattener.visualize('confdist',2,[1 1.5]);
% figure(201);
% f2.visualize('confdist',2,[1 1.5]);
% figure(202);
% f3.visualize('confdist',2,[1 1.5]);
figure(203);
f4.visualize('confdist',2,[1 1.5]);



%% mkdir
dirname=mkdir_for_results(load_name);
%%
fname=[dirname '/' load_name '_before'];
flattener.save([fname]);
S1=flattener.getX3D([load_name '_before'],'confdist',[1.5,1]);
S=x3dScene(S1);
f=fopen([fname  '.x3d'],'w');
fprintf(f,S);
fclose(f);

fname=[dirname '/' load_name '_after'];
f4.save([fname]);
S1=f4.getX3D([load_name '_after'],'confdist',[1.5,1]);
S=x3dScene(S1);
f=fopen([fname  '.x3d'],'w');
fprintf(f,S);
fclose(f);
