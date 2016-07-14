addpath(genpath(['toolbox_graph']));
addpath(genpath(['yalmip']));
addpath(genpath(['MeshLibrary']));
inds_s=[1
    1000
    100
    2000];

M_orig_s=TriangleMesh('off','horse_5k.off');
inds_t=[10
    1500
    300
    2400];
singularities=[2 2 2];
M_orig_t=TriangleMesh('off','horse_5k.off');

[fs,ft]=flatten_compatible(M_orig_s,M_orig_t,inds_s,inds_t,singularities);
% [V_t,T_t,pathPairs_t,vcol_t,M_cut_t]=flatten(M_orig_t,inds_t,singularities);

%%
close all
figure(1);
visualizeNew( V_s,T_s,pathPairs_s,zeros(length(T_s),1),vcol_s);
figure(2);
visualizeNew( V_t,T_t,pathPairs_t,zeros(length(T_t),1),vcol_t);

map=FlatteningsToMap(fs,ft);

% gS={};
% T={};
% pathEnds=[];
% for i=1:length(pathPairs_s)
%     p1=V_s(pathPairs_s{i}(1,1),:);
%     p2=V_s(pathPairs_s{i}(end,1),:);
%     q1=V_s(pathPairs_s{i}(1,2),:);
%     q2=V_s(pathPairs_s{i}(end,2),:);
%     [A,t]=computeSimilarity(p1,p2,q1,q2);
%     gS{i}=A;
%     T{i}=t';
%     pathEnds=[pathEnds;reshape(pathPairs_s{i}([1 end],:),4,1) reshape(pathPairs_t{i}([1 end],:),4,1)];
% end
% pathEnds=unique(pathEnds,'rows');
% source=[];
% source.TR=triangulation(T_s,V_s);
% source.Y=V_s;
% source.pathPairs=pathPairs_s;
% target=[];
% target.TR=triangulation(T_t,V_t);
% target.Y=V_t;
% target.pathPairs=pathPairs_t;
% lifter1to2=RobustLifter(source,target,gS,T,pathEnds);
% lifter1to2.lift();
% lifter2to1=RobustLifter(target,source,gS,T,pathEnds(:,[2 1]));
% lifter2to1.lift();
% map=UncutSurfMap(M_orig_s,M_orig_t,M_cut_s,M_cut_t,lifter1to2,lifter2to1);

%%

disp('============= Generating textured models ===========');
UVGenerator=VisibleUVGenerator;
UVGenerator.SetSource(M_orig_s);
UVGenerator.SetTarget(M_orig_t);
UVGenerator.SetBC(map.barCoords2to1);
UVGenerator.FlipNormals(1);
UVGenerator.SetThreshold(0.5)
pos=[1 0 0;
    -1 0 0;
    0 1 0;
    0 -1 0;
    0 0 1;
    0 0 -1;
    ];
dirname='.';
v_scalars1=[];
v_scalars2=[];
v_scalars1.org_x=M_orig_s.V(1,:);
v_scalars1.org_y=M_orig_s.V(2,:);
v_scalars1.org_z=M_orig_s.V(3,:);
mapped=map.barCoords1to2*M_orig_t.V';
v_scalars1.mapped_x=mapped(:,1);
v_scalars1.mapped_y=mapped(:,2);
v_scalars1.mapped_z=mapped(:,3);

v_scalars2.org_x=M_orig_t.V(1,:);
v_scalars2.org_y=M_orig_t.V(2,:);
v_scalars2.org_z=M_orig_t.V(3,:);
mapped=map.barCoords2to1*M_orig_s.V';
v_scalars2.mapped_x=mapped(:,1);
v_scalars2.mapped_y=mapped(:,2);
v_scalars2.mapped_z=mapped(:,3);


for i=1:size(pos,1)
    
    [UVSource,UVTarget]=UVGenerator.GenerateUVs(pos(i,:));
    UVSource(isnan(UVSource))=100;
    UVTarget(isnan(UVTarget))=100;
    
    output_vtk_surf([dirname sprintf('/textured1_%d_%d_%d.vtk',pos(i,:))],M_orig_s.V',M_orig_s.F',v_scalars1,[],UVSource);
    output_vtk_surf([dirname sprintf('/textured2_%d_%d_%d.vtk',pos(i,:))],M_orig_t.V',M_orig_t.F',v_scalars2,[],UVTarget);
end