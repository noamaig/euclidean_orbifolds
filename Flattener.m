classdef Flattener < handle
    % Object that generates the embedding.
    
    properties
        M_orig; %struct of original mesh
        inds; %indices of the cones
        singularities; %cone angles
        flat_V; %will hold flat vertices after solcing
        flat_T; %the triangles of the mesh after cutting (if needed)
        M_cut; %the cut mesh
        flipped; %boolean - true for flipped (det<0) tris
        V2A; %matrix 
        As; %the affine transforatmion of each triangle
        dets; %determinant of linaer trans of each triangle
        frobenius; %frobenius norm of linear trans of each tri
        smin; %smaller singular value of linear trans of each tri
        smax; %larger singular value of linear trans of each tri
        areas; %area of each tri
        %colors for rendering
        cut_colors={[0 0.4 1],[0 0.7 0.3],[ 0.8 0.2 0.0],[0.8,0.6,0]};
        LM_colors={[0.6 0.2 0.1],[1 0.8 0],[1 0 1],[0 0.8 0.4]};
        
        orgBoundary=[]; %keeps the parts of the mesh which are not a result 
        % of cuttin (i.e., originally were boundary).
        isDisc=false; %flag - is this a sphere or disc?
        verbose=true; %should output some info?
    end
    
    methods
        function obj=Flattener(M_orig,inds,singularities,M_cut)
            %M_orig - struct with M_orig.F = (3,|T|) matrix of faces,
            %M_orig.V =(3,|V|) matrix of vertices.
            %inds - the indices of the cones
            %singularities - the cone singularities
            %M_cut - if specified, use this for the cut mesh
            assert(max(inds)<=length(M_orig.V)&&min(inds)>=1,'cone indices should be positive integers <= |V|');
            obj.M_orig=M_orig;
            %indices of the cones
            obj.inds=inds;
            %the singularities to use (specified according to the index of the cones on the ends of the cut) 
            obj.singularities=singularities;
            if nargin>3
                %an already-cut mesh, if already computed
                obj.setCutMesh(M_cut);
            end
            %computing the boundary
            TR=triangulation(obj.M_orig.F',obj.M_orig.V');
            t=(TR.freeBoundary());
            %if there was a boundary, remember this was a disk
            if ~isempty(t)
                obj.orgBoundary=t(:,1);
                if ~all(ismember(inds,obj.orgBoundary))
                    error('in case of a disk orbifold, all cones should be boundary vertices!');
                end
                ind=find(obj.inds(1)==obj.orgBoundary);
                obj.orgBoundary=[obj.orgBoundary(ind:end)' obj.orgBoundary(1:ind-1)']';
            end
            obj.isDisc=~isempty(obj.orgBoundary);
            
        end
        function B=computeAffineMinimizer(obj)
            %return 2x2 matrix B s.t. the embedding \Phi becomes as conformal
            %as possible in the least-squares sense when applying B*\Phi
            obj.computeAs();
            As = obj.V2A*obj.flat_V(:);
            a=As(1:4:end);
            b=As(3:4:end);
            c=As(2:4:end);
            d=As(4:4:end);
            correctAs=zeros(size(As));
            correctAs([1:4:end 2:4:end 3:4:end 4:4:end])=[a b c d];
            B=minimize_global_lscm(correctAs,obj.areas);
        end
        function correctGlobalAffine(obj)
            %modify the embedding so that it minimizes the LSCM energy
            %(applicable only to orbifold IV)
            assert(length(obj.inds)==4,'cannot modify embedding of oribfolds which are not type IV');
            B=obj.computeAffineMinimizer();
            obj.flat_V=obj.flat_V*B';
            obj.fixToAxis();
        end
        function setCutMesh(obj,M_cut)
            %set the cut mesh (see CutMesh.m for the actual structure)
            M_cut.V=obj.M_orig.V(:,M_cut.cutIndsToUncutInds);
            
            if size(M_cut.V,2)~=3
                M_cut.V=M_cut.V';
            end
            if size(M_cut.T,2)~=3
                M_cut.T=M_cut.T';
            end
            obj.M_cut=CutMesh(M_cut.V,M_cut.T,M_cut.pathPairs,M_cut.cutIndsToUncutInds,M_cut.uncutIndsToCutInds);
        end
        
        function cut(obj)
            %cut the mesh before flattening
            if obj.verbose
                fprintf('*** cutting: ');
                tid=tic;
            end
            if length(obj.inds)==3
                root=length(obj.inds);
                fixedPairs=[ones(1,length(obj.inds)-1)*root;1:length(obj.inds)-1]';
            else
                root=1;
                fixedPairs=[1 3;3 4;4 2];
            end
            tree=sparse(fixedPairs(:,1),fixedPairs(:,2),1,length(obj.inds),length(obj.inds));
            cutter=TreeCutter(obj.M_orig.V',obj.M_orig.F',tree,obj.inds,root);
            cutter.cutTree();
            obj.setCutMesh(cutter);
            if obj.verbose
                toc(tid)
            end
        end
        
        function setLinearConstraints(obj,cons,v1,v2,inds)
            %add constraints so that all vertices in inds are sequentially
            %placed on a line between v1 and v2, with spacing proportional
            %to edge lengths
            %measure the distance between consecutive vertices
            d=sqrt(sum(obj.M_cut.V(inds(1:end-1),:)-obj.M_cut.V(inds(2:end),:),2).^2);
            d=[0;cumsum(d)/sum(d)];
            %add the constraint
            for i=1:length(inds)-1
                cons.addConstraint(inds(i),1,v1*(1-d(i))+v2*d(i));
            end
        end
        function flatten(obj,convexBoundary,DIRICHLET)
            %flatten the mesh to one of the orbifolds. In the end the
            %position of each vertex is stored in the property `flat_V`. 
            %INPUT:
            %convexBoundary - omitted or false for free boundary, of one of 
            % the 4 sphere orbifolds as specified in the constructor; true 
            % for the classic Tutte embedding with fixed boundary into a disk; 
            % 'square' for "classic" Tutte on a square with prefixed boundary 
            % map; 'freesquare' for the square disk orbifold; 'freetri' for 
            % the triangle disk orbifold.
            % Dirichlet - true or omitted for cotan weights, false for mvc weights.
            if nargin<2
                convexBoundary=false;
            end
            if nargin<3
                DIRICHLET=true;
            end
           
            if isempty(obj.M_cut)
                if ~obj.isDisc
                    
                    
                    obj.cut();
                else
                    if strcmp(convexBoundary,'square') || strcmp(convexBoundary,'freesquare')
                        obj.orgBoundaryToPaths();
                    elseif strcmp(convexBoundary,'freetri')
                        obj.orgBoundaryToPaths();
                    else
                        error;
                    end
                end
            end
            if obj.verbose
                fprintf('*** flattening: ');
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%
            % Boundary conditions
            %%%%%%%%%%%%%%%%%%%%%%%
            startP=obj.M_cut.uncutIndsToCutInds{obj.inds(1)};
            assert(length(startP)==1);
            
            tid=tic;
            TR=triangulation(obj.M_cut.T,obj.M_cut.V);
            
            
            
            
            cons=PosConstraints(length(obj.M_cut.V));
            
            
            if convexBoundary
                if strcmp(convexBoundary,'square')
                    %set boundary vertices to fixed positions on square.
                    pathEnds=[];
                    for i=1:length(obj.M_cut.pathPairs)
                        pathEnds=[pathEnds obj.M_cut.pathPairs{i}([1 end],:)];
                    end
                    pathEnds=unique(pathEnds);
                    all_binds = TR.freeBoundary();
                    assert(all(all_binds(:,2)==all_binds([2:end,1],1)));
                    all_binds=all_binds(:,1);
                    all_binds=all_binds(end:-1:1);
                    ind=find(all_binds==startP);
                    all_binds=all_binds([ind:end,1:ind-1]);
                    p=find(ismember(all_binds,pathEnds));
                    obj.setLinearConstraints(cons,[-1 1]',[1 1]',all_binds(p(1):p(2)));
                    obj.setLinearConstraints(cons,[1 1]',[1 -1]',all_binds(p(2):p(3)));
                    obj.setLinearConstraints(cons,[1 -1]',[-1 -1]',all_binds(p(3):p(4)));
                    obj.setLinearConstraints(cons,[-1 -1]',[-1 1]',[all_binds(p(4):length(all_binds)); all_binds(1)]);
                elseif strcmp(convexBoundary,'freesquare')
                   %set boundary vertices to lie on infinite lines
                   %supporting the edges of a square
                    pathEnds=[];
                    for i=1:length(obj.M_cut.pathPairs)
                        pathEnds=[pathEnds obj.M_cut.pathPairs{i}([1 end],:)];
                    end
                    pathEnds=unique(pathEnds);
                    all_binds = TR.freeBoundary();
                    all_binds=all_binds(end:-1:1,1);
                    ind=find(all_binds==startP);
                    all_binds=all_binds([ind:end,1:ind-1]);
                    p=find(ismember(all_binds,pathEnds));
                    %                     p=all_binds(p);
                    cons.addConstraint(all_binds(p(1)),1,[-1 1]');
                    cons.addConstraint(all_binds(p(2)),1,[1 1]');
                    cons.addConstraint(all_binds(p(3)),1,[1 -1]');
                    cons.addConstraint(all_binds(p(4)),1,[-1 -1]');
                    for i=p(1)+1:p(2)-1
                        cons.addLineConstraint(all_binds(i),[0 1],1);
                    end
                    
                    for i=p(2)+1:p(3)-1
                        cons.addLineConstraint(all_binds(i),[1 0],1);
                    end
                    
                    
                    for i=p(3)+1:p(4)-1
                        cons.addLineConstraint(all_binds(i),[0 -1],1);
                    end
                    for i=p(4)+1:length(all_binds)
                        cons.addLineConstraint(all_binds(i),[-1 0],1);
                    end
                    
                    
                elseif strcmp(convexBoundary,'freetri')
                    %set boundary vertices to lie on infinite lines
                   %supporting the edges of a triangle
                    pathEnds=[];
                    for i=1:length(obj.M_cut.pathPairs)
                        pathEnds=[pathEnds obj.M_cut.pathPairs{i}([1 end],:)];
                    end
                    pathEnds=unique(pathEnds);
                    all_binds = TR.freeBoundary();
                    all_binds=all_binds(end:-1:1,1);
                    ind=find(all_binds==startP);
                    all_binds=all_binds([ind:end,1:ind-1]);
                    p=find(ismember(all_binds,pathEnds));
                    %                     p=all_binds(p);
                    cons.addConstraint(all_binds(p(1)),1,[0 0]'); %90
                    cons.addConstraint(all_binds(p(2)),1,[0 1]'); %45
                    cons.addConstraint(all_binds(p(3)),1,[1 0]');  %45
                    
                    for i=p(1)+1:p(2)-1
                        cons.addLineConstraint(all_binds(i),[1 0],0);
                    end
                    
                    for i=p(2)+1:p(3)-1
                        cons.addLineConstraint(all_binds(i),normr([1 1]),sqrt(2)/2);
                    end
                    
                    
                    
                    for i=p(3)+1:length(all_binds)
                        cons.addLineConstraint(all_binds(i),[0 1],0);
                    end
                    
                else
                    %set boundary vertices to lie on fixed positions in the
                    %unit square
                    all_binds = TR.freeBoundary();
                    d=sqrt(sum(obj.M_cut.V(all_binds(:,1),:)-obj.M_cut.V(all_binds(:,2),:),2).^2);
                    %want area of 4 (as the square), so pi r^2=4 -->
                    %r=2/sqrt(pi)
                    R=2/sqrt(pi);
                    theta=2*pi*cumsum(d)/sum(d);
                    for i=1:length(all_binds)
                        ind=all_binds(i,1);
                        t=theta(i);
                        cons.addConstraint(ind,1,R*[cos(t) sin(t)]');
                    end
                end
            else
                %a sphere orbifold, boundary is associated by rotations
                %find the indices of the cones on the cut mesh
                pathEnds=[];
                for i=1:length(obj.M_cut.pathPairs)
                    pathEnds=[pathEnds obj.M_cut.pathPairs{i}([1 end],:)];
                end
                pathEnds=unique(pathEnds);
                all_binds = TR.freeBoundary();
                assert(all(all_binds(:,2)==all_binds([2:end,1],1)));
                %holds the ordered list of boundary vertices of cut mesh
                all_binds=all_binds(:,1);
                %find indices of boundary vertices which are cones
                ind=find(all_binds==startP);
                %cyclic shift all_binds so that all_binds(1) is a cone
                all_binds=all_binds([ind:end,1:ind-1]);
                %find the 3 cones
                p=find(ismember(all_binds,pathEnds));
                p=all_binds(p);
                
                
                angs={};
                %positions of cones in the plane
                theta=2*pi*(1:length(p))/length(p)+pi/4;
                coords=[cos(theta)' sin(theta)']*sqrt(2);
                if length(obj.inds)==4 %type IV orbifold
                    tcoords=[0 -0.5;0 0.5];
                elseif all(obj.singularities==[4 4]) %type I orbifold
                    tcoords=[-1 -1;1 1];
                else %type II or III
                    tcoords=coords;
                end
                
                for i=1:length(p) %iterate over cones
                    %find which cone is it according to their order in p.
                    ind=find(obj.M_cut.cutIndsToUncutInds(p(i))==obj.inds);
                    if ind<=2 %for cones 1 or 2
                        %set positional constraint on cone's position in
                        %plane
                        cons.addConstraint(p(i),1,tcoords(ind,:)');
                        %set the angle according to specified singularity
                        angs{p(i)}=obj.singularities(ind);
                    else
                        %other cones will not get an angle
                        angs{p(i)}=[];
                    end
                    if length(obj.inds)==4&& i==2 %if orbifold type IV and it is second cone in flattening, add constraint
                        cons.addConstraint(p(i),1,[1 -0.5]');
                    end
                end
                %add rotational constraints
                for i=1:length(obj.M_cut.pathPairs)
                    %getting the corresponding indices of the two paths on
                    %the two sides of the cut
                    path1=obj.M_cut.pathPairs{i}(:,1);
                    path2=obj.M_cut.pathPairs{i}(:,2);
                    
                    %making sure that if paths share the same vertex on
                    %both sides, it's the first vertex.
                    sign=-1;
                    if path1(end)==path2(end)
                        path1=path1(end:-1:1);
                        path2=path2(end:-1:1);
                        sign=1;
                    end
                    %get the angle associated with this path
                    ang=angs{path1(1)};
                    
                    %if there is no associated angle, just force
                    %translational constraints
                    if isempty(ang)
                        ang=1;
                    end
                    
                    if ~isempty(ang)
                        %correct the angle in case we reversed ordering
                        %previously 
                        ang=sign*ang;
                        %rotation matrix according to ang
                        R=[cos(2*pi/ang) -sin(2*pi/ang);sin(2*pi/ang) cos(2*pi/ang)];
                        %add the constrain that the two paths are identical
                        %up to a rotation
                        cons.addTransConstraints(path1,path2,R)
                    end
                end
            end
            if obj.verbose
                fprintf('constraint generation: %f seconds\n',toc(tid));
                tid=tic;
            end
            
            %%%%%%%%%%
            % Energy (laplacain)
            %%%%%%%%%%
            if DIRICHLET
                L=cotmatrix(obj.M_cut.V,obj.M_cut.T);
                m=min(min(tril(L,-1)));
                if m<0
                    warning('Mesh is not Delaunay!!');
                    clamp=1e-2;
                    fprintf('clamping negative weights to %f\n',clamp);
                    L(L<0)=clamp;
                    % now fix the laplacian
                    inds=sub2ind(size(L),1:length(L),1:length(L));
                    L(inds)=0;
                    L(inds)=-sum(L);
                end
            else
                L=mean_value_laplacian(obj.M_cut.V,obj.M_cut.T);
            end
            %duplicating laplacian to work on each coordinate.
            RealL=sparse(size(L,1)*2,size(L,2)*2);
            RealL(1:2:end,1:2:end)=L;
            RealL(2:2:end,2:2:end)=L;
            L=RealL;
            if obj.verbose
                fprintf('compute: %f seconds\n',toc(tid));
                
                tidc=tic;
            end
            %compute the flattening by solving the boundary conditions
            %while satisfying the convex combination property with L
            x=computeFlattening(cons.A,cons.b,L);
            if obj.verbose
                fprintf('lin solve: %f seconds\n',toc(tidc));
            end
            %create the images of each vertex
            X=x(1:2:end);
            Y=x(2:2:end);
            
            obj.flat_V=[X Y];
            obj.flat_T=obj.M_cut.T;
            
            
            
        end
        function orgBoundaryToPaths(obj)
            %in case the mesh is a disk, create a "fake" boundary of correspondences 
            %(to fit %our internal structure)
            pathPairs={};
            for i=1:length(obj.inds)-1
                s=find(obj.orgBoundary==obj.inds(i));
                e=find(obj.orgBoundary==obj.inds(i+1));
                p=obj.orgBoundary(s:e);
                pathPairs{i}=[p p];
            end
             s=find(obj.orgBoundary==obj.inds(end));
            p=obj.orgBoundary([s:length(obj.orgBoundary) 1]);
            pathPairs{end+1}=[p p];
            
            obj.M_cut=CutMesh(obj.M_orig.V',obj.M_orig.F',pathPairs,1:length(obj.M_orig.V),num2cell(1:length(obj.M_orig.V)));
        end
        
        function computeV2A(obj)
            
            % calculate map between 2d vertices to matrix mapping each
            % triangle induced by the placement of 2d vertices.
            [obj.V2A,obj.areas] = getFlatteningDiffCoefMatrix(obj.M_cut.V,obj.M_cut.T); 
            
        end
        
        
        
        function computeAs(obj)
            % calculate map between 2d vertices to differentials
            if isempty(obj.V2A)
                obj.computeV2A();
            end
            
            %obj.As = permute(reshape(obj.V2A*obj.flat_V,2,[],2),[1 3 2]);
            obj.As = reshape(obj.V2A*obj.flat_V(:),2,2,[]);
            
        end
        function fixToAxis(obj)
            %make sure the type IV orbifold is axis-aligned.
            p=obj.flat_V(obj.M_cut.uncutIndsToCutInds{obj.inds(3)},:);
            p1=p(1,:);
            p2=p(2,:);
            [ M,t ] = computeSimilarity( p1,p2,[1 0],[-1 0] );
            obj.flat_V=bsxfun(@plus,obj.flat_V*M',t');
        end
        function fixToGrid(obj,scale)
            %make sure the type IV orbifold lands exactly on an integer
            %grid.
            obj.fixToAxis();
            obj.flat_V=obj.flat_V*scale;
            ind1=obj.M_cut.uncutIndsToCutInds{obj.inds(3)}(1);
            ind2=obj.M_cut.uncutIndsToCutInds{obj.inds(3)}(2);
            
            X=obj.flat_V;
            if X(ind1,1)>X(ind2,1)
                ind3=ind1;
                ind1=ind2;
                ind2=ind3;
            end
            X(:,1)=X(:,1)-X(ind1,1);
            m=X(ind2,1);
            X=2*round(m)*X/m;
            
            
            ind3=obj.M_cut.uncutIndsToCutInds{obj.inds(4)}(1);
            ind4=obj.M_cut.uncutIndsToCutInds{obj.inds(4)}(2);
            p1=X(ind3,:);
            p2=X(ind4,:);
            p3=X(ind2,:);
            q1=round(p1);
            q2=round(p2);
            q3=p3;
            [ A,t ] = affineTransFrom3Points( p1,p2,p3,q1,q2,q3 );
            X=bsxfun(@plus,X*A',t');
            obj.flat_V=X;
            
        end
    end
    
end


