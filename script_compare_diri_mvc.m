rng(1);
addpath(genpath(['toolbox_graph']));
addpath(genpath(['yalmip']));
addpath(genpath(['MeshLibrary']));
load_name='camel';
% load_name='gargoyle';
f_diri=loadDataToFlattener([load_name '.mat']);
f_diri.flatten();


f_mvc=loadDataToFlattener([load_name '.mat']);
f_mvc.setCutMesh(f_diri.M_cut);
f_mvc.flatten(false,false);

%%
%compute tutte's distortion first as it fails if mesh is not delaunay also
f_mvc.computeDistortion();

f_diri.computeDistortion();
clf
subplot(1,2,1);
f_diri.visualize();

subplot(1,2,2);
f_mvc.visualize();
return
%% mkdir
dirname=mkdir_for_results(load_name);
%%
fname=[dirname '/' load_name];
f_ours.save([fname '_ours']);
S1=f_ours.getX3D(['our_' load_name],'frobenius',[maxE,0]);
S=x3dScene(S1);
f=fopen([fname '_ours.x3d'],'w');
fprintf(f,S);
fclose(f);

f_square.save([fname '_square']);
S2=f_square.getMeshX3D(['square_' load_name],'frobenius',[maxE,0]);
S=x3dScene(S2);
f=fopen([fname '_square.x3d'],'w');
fprintf(f,S);
fclose(f);

f_circle.save([fname '_circle']);
S2=f_circle.getMeshX3D(['circle_' load_name],'frobenius',[maxE,0]);
S=x3dScene(S2);
f=fopen([fname '_circle.x3d'],'w');
fprintf(f,S);
fclose(f);