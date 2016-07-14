rng(1);
addpath(genpath(['toolbox_graph']));
addpath(genpath(['yalmip']));
addpath(genpath(['MeshLibrary']));
fname='eros';
f_ours=loadDataToFlattener([fname '.mat']);
f_ours.flatten();


f_tutte=loadDataToFlattener([fname '.mat']);
f_tutte.flatten('disk');
%%
f_ours.computeDistortion();
f_tutte.computeDistortion();
energy_ours=sum(f_ours.frobenius.*f_ours.areas)/sum(f_ours.areas);
energy_tutte=sum(f_tutte.frobenius.*f_tutte.areas)/sum(f_tutte.areas);
maxE=prctile(f_tutte.frobenius,95);
%%
figure(200);
clf
subplot(1,2,2);
f_ours.visualize('dirichlet');
title(sprintf('%f',energy_ours));
axis off

ax=axis;
ax=ax*1.1;
axis(ax);
subplot(1,2,1);
f_tutte.visualize('dirichlet');
title(sprintf('%f',energy_tutte));
axis off;

colormap(jet);
caxis([0 maxE]);
axis(ax);
%%
S1=f_ours.getX3D(['our_' fname],[],[maxE,0]);
S2=f_tutte.getX3D(['tutte_' fname],[],[maxE,0]);
S=d3xScene(S1,S2);
f=fopen([fname '.x3d'],'w');
            fprintf(f,S);
            fclose(f);
