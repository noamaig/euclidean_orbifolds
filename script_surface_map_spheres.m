%% =====================================================================
%  ===                      Setting up stuff                         ===
%  =====================================================================

%Read the mesh
[V1,T1]=read_off('max_map.off');V1=V1';T1=T1';
[V2,T2]=read_off('julius_map.off');V2=V2';T2=T2';

%Some color-related variables - no need to concern yourself with these :)
cone_colors=[1 0.8 0;0.7 0 1; 0 0.5 0.8;0 0 0.5];

%The cone positions and the choice of orbifold structure, defining the
%desired embedding

%this uses the type IV orbifold which can satisfy 4 given correspondences
%but pays for that by being "only" quasiconformal and not conformal.
% orbifold_type=OrbifoldType.Parallelogram;
% % 
% cones1=[  25372
%     12611
%     5868
%     12618];
% cones2=[  16813
%     17530
%     17692
%     16986
%     ];

%this uses only 3 points and hence yields a conformal map between the surfaces
orbifold_type=OrbifoldType.Square;
% 
cones1=[  25372
    12611
    
    12618];
cones2=[  16813
    17530
    
    16986
    ];



%% =======================================================================
%  =======       The actual algorithm! cutting and flattening      =======
%  =======================================================================

BC=map_spheres(V1,T1,V2,T2,cones1,cones2,orbifold_type);
V_mapped=BC*V2;


%% =======================================================================
% =======       And now a nice visualization....
%=========================================================================

tr=triangulation(T1,V1);
VN=tr.faceNormal();
VN=(VN+1)/2;
VN(all(VN==0,2),:)=1;
VN=normr(VN);


fig=figure(1);
    clf

set(gcf, 'Position', get(0,'Screensize'));


M=struct('cdata',[],'colormap',[]);
step=0.02;
count=0;
ax=[];
for c=[0:step:1-step 1:-step:0]
    clf
    V=V_mapped*c+(1-c)*V1;
    patch('faces',T1,'vertices',V,'facecolor','flat','edgecolor','k','edgealpha',1,'FaceVertexCDATA',VN);
    campos([-39.5547   51.5112 -130.9298]);
    camup([0.1032    0.9349    0.3395]);
    axis([ -10.4249    8.5982   -6.1319    5.0566   -7.3091    6.6766])
    campos([-39.5547   51.5112 -130.9298]);
    camup([0.1032    0.9349    0.3395]);
    % campos([58.6269   45.4568 -126.0013]);
    % camup([ -0.1328    0.9493    0.2850]);
    axis off
    camlight
    hold on
    fig.PaperPositionMode = 'auto';
    for i=1:length(cones2)
        VV=V(cones1(i),:);
        scatter3(VV(:,1),VV(:,2),VV(:,3),60,cone_colors(i,:),'fill');
    end
    drawnow
    pause(0.001);
    count=count+1;

    

end
