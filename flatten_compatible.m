function [fs,ft]=flatten_compatible(M_orig_s,M_orig_t,inds_s,inds_t,singularities)
% addpath(genpath(['toolbox_graph']));
% addpath(genpath(['yalmip']));
% addpath(genpath(['MeshLibrary']));
DIRICHLET=1;
% rng(2);
tid=tic;


tid=tic;
M_to_cut_s=TriangleMesh('VF',M_orig_s.V,M_orig_s.F);
M_to_cut_t=TriangleMesh('VF',M_orig_t.V,M_orig_t.F);
% M_to_cut=TriangleMesh('VF',M_orig.V+diag([100 1 1])*rand(size(M_orig.V))*3200,M_orig.F);
% M_to_cut=TriangleMesh('VF',diag([1 5 1])*M_orig.V+rand(size(M_orig.V))*200,M_orig.F);

cutter=CompatibleMeshCutter();
cutter.setMesh1(M_to_cut_s,inds_s);
cutter.setMesh2(M_to_cut_t,inds_t);
fixedPairs=[ones(1,length(inds_s)-1)*length(inds_s);1:length(inds_s)-1]';
if length(inds_s)==3
    root=length(inds_s);
    fixedPairs=[ones(1,length(inds_s)-1)*root;1:length(inds_s)-1]';
else
    root=1;
    fixedPairs=[1 3;3 4;4 2];
end
cutter.setFixedPairs(fixedPairs);
cutter.Cut();
[M_s,M_t]=cutter.GatherResults;

fs=Flattener(M_orig_s,inds_s,singularities,M_s);
ft=Flattener(M_orig_t,inds_t,singularities,M_t);
fs.flatten(false,false);
ft.flatten(false,false);
% function [V,T,pathPairs,vcol,M_cut]=helper(M,M_orig,inds,singularities)
% M.V=M_orig.V(:,M.New2Old);
% 
% M_cut=M;
% M_cut.V=M.V';
% M_cut.T=M.F';
% pathPairs=convertFromCuttingToMesh(M);
% pathPairsOrdered={};
% for i=1:length(pathPairs)
%     ind=find(inds==M.New2Old(pathPairs{i}(end,1)));
%     pathPairsOrdered{ind}=pathPairs{i};
% end
% M_cut.pathPairs=pathPairsOrdered;
% startP=M_cut.pathPairs{1}(end,1);
% assert(startP==M_cut.pathPairs{1}(end,2));
% % fprintf('cutting: %f seconds\n',toc(tid));
% tid=tic;
% TR=triangulation(M_cut.T,M_cut.V);
% pathPairs=M_cut.pathPairs;
% pathEnds=[];
% for i=1:length(M_cut.pathPairs)
%     pathEnds=[pathEnds M_cut.pathPairs{i}([1 end],:)];
% end
% pathEnds=unique(pathEnds);
% all_binds = TR.freeBoundary();
% assert(all(all_binds(:,2)==all_binds([2:end,1],1)));
% all_binds=all_binds(:,1);
% ind=find(all_binds==startP);
% all_binds=all_binds([ind:end,1:ind-1]);
% p=find(ismember(all_binds,pathEnds));
% p=all_binds(p);
% cons=PosConstraints(length(M_cut.V));
% 
% % cons.addConstraint(p,0.25,0);
% % cons.addConstraint(p(1),1,1);
% 
% 
% P={};
% angs={};
% theta=2*pi*(1:length(p))/length(p);
% 
% % coords=[
% %         0 1;
% %         0 2;
% %         1 2;
% %         1 1;
% %         1 0; 0 0];
% %     coords=coords(end:-1:1,:);
% coords=[cos(theta)' sin(theta)'];
% if length(inds)==4
%     tcoords=[0 0;0 0;1 0;0 0;1 1];
% else
%     tcoords=coords;
% end
% 
% for i=1:length(p)
%     
%     if ismember(i,[1 3 5])
%         cons.addConstraint(p(i),1,tcoords(i,:)');
%     end
%     P{p(i)}=coords(i,:);
%     ind=find(M.New2Old(p(i))==inds);
%     if ind<=length(singularities)
%         angs{p(i)}=singularities(ind);
%     end
% end
% for i=1:length(M_cut.pathPairs)
%     path1=M_cut.pathPairs{i}(:,1);
%     path2=M_cut.pathPairs{i}(:,2);
%     sign=-1;
%     if path1(end)==path2(end)
%         path1=path1(end:-1:1);
%         path2=path2(end:-1:1);
%         sign=1;
%     end
%     %     p1=P{path1(1)};
%     %     p2=P{path1(end)};
%     %     q1=P{path2(1)};
%     %     q2=P{path2(end)};
%     
%     %     %adding constraints on 2nd point
%     %     v=(p2-p1)/length(path1);
%     %     ind=path1(2);
%     %     cons.addConstraint(ind,1,(p1+v)');
%     %     ind=path1(end-1);
%     %     cons.addConstraint(ind,1,(p2-v)');
%     %%%%%%%%%%
%     
%     %     [ R,t ] = computeSimilarity( p1,p2,q1,q2 );
%     %     R=round(R);
%     %     t=round(t);
%     ang=angs{path1(1)};
%     if ~isempty(ang)
%         ang=sign*ang;
%         R=[cos(2*pi/ang) -sin(2*pi/ang);sin(2*pi/ang) cos(2*pi/ang)];
%         cons.addTransConstraints(path1,path2,R)
%     end
% end
% fprintf('constraint generation: %f seconds\n',toc(tid));
% tid=tic;
% 
% [V2A,areas] = getFlatteningDiffCoefMatrix(M_cut.V',M_cut.T'); % calculate map between 2d vertices to differentials
% 
% fprintf('flattening coef matrix: %f seconds\n',toc(tid));
% tid=tic;
% %
% 
% %% Dirichlet Laplacian...
% if DIRICHLET
%     RealV2A=sparse(size(V2A,1)*2,size(V2A,2)*2);
%     RealV2A(1:2:end,1:2:end)=V2A;
%     RealV2A(2:2:end,2:2:end)=V2A;
%     % X=V(T,1)';
%     % Y=V(T,2)';
%     % areas=ones(length(T),1);%polyarea(X,Y);
%     areas4=repmat(reshape(areas,1,length(areas)),4,1);
%     areas4=areas4(:);
%     % areas4=spdiags(areas4',1:length(areas4),length(areas4),length(areas4));
%     areas4= accumarray([(1:length(areas4))' (1:length(areas4))'],areas4,[],[],[],true);
%     L=RealV2A'*areas4*RealV2A;
% else
%     L = compute_combinatorial_laplacian( triangulation2adjacency(TR.ConnectivityList) );
%     RealL=sparse(size(L,1)*2,size(L,2)*2);
%     RealL(1:2:end,1:2:end)=L;
%     RealL(2:2:end,2:2:end)=L;
%     L=RealL;
% end
% x=computeFlattening(cons.A,cons.b,L);
% fprintf('compute: %f seconds\n',toc(tid));
% tid=tic;
% 
% X=x(1:2:end);
% Y=x(2:2:end);
% 
% V=[X Y];
% 
% 
% As=V2A*V;
% Astack = permute(reshape(As,2,[],2),[1 3 2]);
% %%
% dets=zeros(size(Astack,3),1);
% dists=zeros(size(Astack,3),1);
% for i=1:size(Astack,3)
%     A=Astack(:,:,i);
%     dets(i)=det(A);
%     dists(i)=norm(A,'fro');
% end
% flipped=dets<-1e-8;
% % close all;
% figure(501);
% clf;
% TR=triangulation(M_orig.F',M_orig.V');
% vn = vertexNormal(TR);
% vcol=(vn(M.New2Old,:)+1)/2;
% vcol=rgb2hsv(vcol);
% vcol(:,2)=0.4;
% vcol(:,3)=1;
% vcol=hsv2rgb(vcol);
% visualizeNew( V,M_cut.T,M_cut.pathPairs,flipped,vcol);
% fprintf('drawing: %f seconds\n',toc(tid));
% assert(~any(flipped));
% output_vtk_surf('test.vtk',M_cut.V,M_cut.T,[],[],V);
% T=M_cut.T;
% end
end