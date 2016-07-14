rng(1);
addpath(genpath(['toolbox_graph']));
addpath(genpath(['yalmip']));
addpath(genpath(['MeshLibrary']));
load_name='arma';
% load_name='gargoyle';
f_ours=loadDataToFlattener([load_name '.mat']);
f_ours.singularities=[4 4 ];
f_ours.flatten();


f_square=loadDataToFlattener([load_name '.mat']);
f_square.setCutMesh(f_ours.M_cut);
f_square.flatten('square');

f_circle=loadDataToFlattener([load_name '.mat']);
f_circle.setCutMesh(f_ours.M_cut);
f_circle.flatten('circle');
%%
%compute tutte's distortion first as it fails if mesh is not delaunay also
f_square.computeDistortion();

f_ours.computeDistortion();

energy_ours=sum(f_ours.frobenius.^2.*f_ours.areas)/sum(f_ours.areas);
energy_square=sum(f_square.frobenius.^2.*f_square.areas)/sum(f_square.areas);
energy_circle=sum(f_circle.frobenius.^2.*f_circle.areas)/sum(f_circle.areas);
%% visualize
max_prct=99.8;
min_prct=0;
maxE=min([prctile(f_ours.frobenius,max_prct),prctile(f_square.frobenius,max_prct),prctile(f_square.frobenius,max_prct)]);
minE=max([prctile(f_ours.frobenius,min_prct),prctile(f_square.frobenius,min_prct),prctile(f_square.frobenius,min_prct)]);


figure(200);
clf
subplot(1,3,3);
f_circle.visualize('dirichlet',3);
title(sprintf('%f',energy_circle));
axis off;
colormap(jet);
caxis([minE maxE]);

axis equal;

% ax=axis;
% ax=ax*1.22;
% axis(ax);

subplot(1,3,2);
f_square.visualize('dirichlet',3);
title(sprintf('%f',energy_square));
axis off;
colormap(jet);
caxis([minE maxE]);
% axis(ax);

subplot(1,3,1);
f_ours.visualize('dirichlet',3);
title(sprintf('%f',energy_ours));
axis off
caxis([minE maxE]);colormap(jet);

% axis(ax);
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