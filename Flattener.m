classdef Flattener < handle
    %FLATTENER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        M_orig;
        inds;
        singularities;
        flat_V;
        flat_T;
        %         vcol;
        %         ncol;
        M_cut;
        flipped;
        V2A;
        As;
        dets;
        frobenius;
        smin;
        smax;
        areas;
        cut_colors={[0 0.4 1],[0 0.7 0.3],[ 0.8 0.2 0.0],[0.8,0.6,0]};
        LM_colors={[0.6 0.2 0.1],[1 0.8 0],[1 0 1],[0 0.8 0.4]};
        orgBoundary=[];
        isDisc=false;
    end
    
    methods
        function obj=Flattener(M_orig,inds,singularities,M_cut)
            %             assert(length(inds)==length(singularities)+1);
            obj.M_orig=M_orig;
            obj.inds=inds;
            obj.singularities=singularities;
            if nargin>3
                obj.setCutMesh(M_cut);
            end
            TR=triangulation(obj.M_orig.F',obj.M_orig.V');
           t=(TR.freeBoundary());
           if ~isempty(t)
           obj.orgBoundary=t(:,1);
           end
           obj.isDisc=~isempty(obj.orgBoundary);
            
        end
        function B=computeAffineMinimizer(obj,opt)
            obj.computeAs();
            
            As = obj.V2A*obj.flat_V(:);
            a=As(1:4:end);
            b=As(3:4:end);
            c=As(2:4:end);
            d=As(4:4:end);
            correctAs=zeros(size(As));
            correctAs([1:4:end 2:4:end 3:4:end 4:4:end])=[a b c d];
            B=minimize_global_lscm(correctAs,obj.areas,opt);
        end
        function correctGlobalAffine(obj,opt)
            if nargin==1
                opt='LSCM';
            end
            B=obj.computeAffineMinimizer(opt);
            obj.flat_V=obj.flat_V*B';
            obj.fixToAxis();
            obj.computeDistortion(true);
        end
        function setCutMesh(obj,M_cut)
            if isa(M_cut,'TreeCutter')
                M_cut.V=obj.M_orig.V(:,M_cut.new2old);
                obj.M_cut.V=M_cut.V;
                obj.M_cut.T=M_cut.T;
                if size(obj.M_cut.V,2)~=3
                    obj.M_cut.V=obj.M_cut.V';
                end
                if size(obj.M_cut.T,2)~=3
                    obj.M_cut.T=obj.M_cut.T';
                end
                obj.M_cut.pathPairs=M_cut.pathPairs;
                obj.M_cut.New2Old=M_cut.new2old;
                obj.M_cut.Old2New=M_cut.old2new;
            else
                M_cut.V=obj.M_orig.V(:,M_cut.New2Old);
                obj.M_cut.V=M_cut.V;
                obj.M_cut.T=M_cut.F;
                if size(obj.M_cut.V,2)~=3
                    obj.M_cut.V=obj.M_cut.V';
                end
                if size(obj.M_cut.T,2)~=3
                    obj.M_cut.T=obj.M_cut.T';
                end
                obj.M_cut.pathPairs=M_cut.pathPairs;
                obj.M_cut.New2Old=M_cut.New2Old;
                obj.M_cut.Old2New=M_cut.Old2New;
            end
        end
        
       
       
        
        function cut(obj)
            fprintf('*** cutting: ');
            tid=tic;
            if length(obj.inds)==3
                root=length(obj.inds);
                fixedPairs=[ones(1,length(obj.inds)-1)*root;1:length(obj.inds)-1]';
            else
                root=1;
                fixedPairs=[1 3;3 4;4 2];
            end
            % cutter=SingleMeshCutter(M_to_cut,inds);
            %cutter.setFixedPairs(fixedPairs);
            
            tree=sparse(fixedPairs(:,1),fixedPairs(:,2),1,length(obj.inds),length(obj.inds));
            %for david
%                         cutter=TreeCutter(obj.M_orig.V'+rand(size(obj.M_orig.V'))*diag([100 1 1])*3200,obj.M_orig.F',tree,obj.inds,root);
                        cutter=TreeCutter(obj.M_orig.V'+randn(size(obj.M_orig.V'))*3200,obj.M_orig.F',tree,obj.inds,root);
%             cutter=TreeCutter(obj.M_orig.V',obj.M_orig.F',tree,obj.inds,root);
            cutter.cutTree();
            % M=cutter.GatherResults;
            
            
            obj.setCutMesh(cutter);
            toc(tid)
        end
        function vprops=face2vertex_int(obj,property,areas)
            inds=reshape(obj.M_orig.F,1,3*length(obj.M_orig.F))';
            areas=kron(areas,[1; 1; 1]);
            props=kron(property,[1; 1; 1]);
            vareas=accumarray(inds,areas);
            vprops=accumarray(inds,areas.*props);
            vprops=vprops./vareas;
            
        end
        function setLinearConstraints(obj,cons,v1,v2,inds)
            d=sqrt(sum(obj.M_cut.V(inds(1:end-1),:)-obj.M_cut.V(inds(2:end),:),2).^2);
            d=[0;cumsum(d)/sum(d)];
            for i=1:length(inds)-1
                cons.addConstraint(inds(i),1,v1*(1-d(i))+v2*d(i));
            end
        end
        function flatten(obj,convexBoundary,DIRICHLET)
            %convexBoundary - 1. omitted or false for free boundary
            %2. 'square' for arrangment on a square
            %3. true for disc
            if nargin<2
                convexBoundary=false;
            end
            if nargin<3
                DIRICHLET=true;
            end
            %show the chosen indices
            fprintf('%d ',obj.inds');
            
            if isempty(obj.M_cut)
                if ~obj.isDisc
                
                     
                    obj.cut();
                else
                    if strcmp(convexBoundary,'square') || strcmp(convexBoundary,'freesquare')
                        obj.orgBoundaryToPaths(4);
                    elseif strcmp(convexBoundary,'freetri')
                        obj.orgBoundaryToPaths(3);
                    else
                        error;
                    end
                end
            end
            fprintf('*** flattening: ');
            
                startP=obj.M_cut.Old2New{obj.inds(1)};
            assert(length(startP)==1);
          
            tid=tic;
            TR=triangulation(obj.M_cut.T,obj.M_cut.V);
            
            
            
            
            cons=PosConstraints(length(obj.M_cut.V));
            
            
            if convexBoundary
                if strcmp(convexBoundary,'square')
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
%                     if obj.isDisc
%                         obj.orgBoundaryToPaths(4);
%                     end
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
%                     if obj.isDisc
%                         obj.orgBoundaryToPaths(4);
%                     end
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
                pathEnds=[];
                for i=1:length(obj.M_cut.pathPairs)
                    pathEnds=[pathEnds obj.M_cut.pathPairs{i}([1 end],:)];
                end
                pathEnds=unique(pathEnds);
                all_binds = TR.freeBoundary();
                assert(all(all_binds(:,2)==all_binds([2:end,1],1)));
                all_binds=all_binds(:,1);
                ind=find(all_binds==startP);
                all_binds=all_binds([ind:end,1:ind-1]);
                p=find(ismember(all_binds,pathEnds));
                p=all_binds(p);
                
                %                 P={};
                angs={};
                theta=2*pi*(1:length(p))/length(p)+pi/4;
                coords=[cos(theta)' sin(theta)']*sqrt(2);
                if length(obj.inds)==4
                    %                     tcoords=[0 -0.5;1 -0.5;1 0;0 0.5;1 1];
                    tcoords=[0 -0.5;0 0.5];
                elseif all(obj.singularities==[4 4])
                    tcoords=[-1 -1;1 1];
                else
                    tcoords=coords;
                end
%                 tcoords=[0 -0.5;0 0.5];warning('ifproblemremovethis');
                for i=1:length(p)
                    ind=find(obj.M_cut.New2Old(p(i))==obj.inds);
                    if ind<=2
                        cons.addConstraint(p(i),1,tcoords(ind,:)');
                        
                        angs{p(i)}=obj.singularities(ind);
                    else
                        angs{p(i)}=[];
                    end
                    if length(obj.inds)==4&& i==2
                        cons.addConstraint(p(i),1,[1 -0.5]');
                    end
                end
                for i=1:length(obj.M_cut.pathPairs)
                    path1=obj.M_cut.pathPairs{i}(:,1);
                    path2=obj.M_cut.pathPairs{i}(:,2);
                    sign=-1;
                    if path1(end)==path2(end)
                        path1=path1(end:-1:1);
                        path2=path2(end:-1:1);
                        sign=1;
                    end
                    
                    ang=angs{path1(1)};
                    if isempty(ang)
                        ang=1;
                    end
                    
                    if ~isempty(ang)
                        ang=sign*ang;
                        R=[cos(2*pi/ang) -sin(2*pi/ang);sin(2*pi/ang) cos(2*pi/ang)];
                        cons.addTransConstraints(path1,path2,R)
                    end
                end
            end
            fprintf('constraint generation: %f seconds\n',toc(tid));
            tid=tic;
            
            
            
            %% Dirichlet Laplacian...
            if DIRICHLET
                L=cotmatrix(obj.M_cut.V,obj.M_cut.T);
                m=min(min(tril(L,-1)));
                if m<0
                    warning('Mesh is not Delaunay!!');
                end
            else
                L=mean_value_laplacian(obj.M_cut.V,obj.M_cut.T);
            end
            
                RealL=sparse(size(L,1)*2,size(L,2)*2);
                RealL(1:2:end,1:2:end)=L;
                RealL(2:2:end,2:2:end)=L;
                L=RealL;
                        fprintf('compute: %f seconds\n',toc(tid));

            tidc=tic;
            x=computeFlattening(cons.A,cons.b,L);
            fprintf('lin solve: %f seconds\n',toc(tidc));
            tid=tic;
            
            X=x(1:2:end);
            Y=x(2:2:end);
            
            obj.flat_V=[X Y];
            
            %%
            
            %
            obj.flat_T=obj.M_cut.T;
            
            
            
        end
        function orgBoundaryToPaths(obj,N)
            pathPairs={};
            inds=obj.orgBoundary;
            d=sqrt(sum(obj.M_orig.V(:,inds(1:end-1))'-obj.M_orig.V(:,inds(2:end))',2).^2);
            d=[0;cumsum(d)/sum(d)];
           
            inds=[1];
            for i=1:N-1
                ind=find(d>i/N,1);
                inds(end+1)=ind;
            end
%             inds(end+1)=length(obj.orgBoundary);
            
            for i=1:length(inds)-1
                p=obj.orgBoundary(inds(i):inds(i+1));
                pathPairs{i}=[p p];
            end
             p=obj.orgBoundary([inds(i+1):length(obj.orgBoundary) 1]);
                pathPairs{end+1}=[p p];
            obj.inds=obj.orgBoundary(inds);
            M_cut=[];
            M_cut.V=obj.M_orig.V';
            M_cut.T=obj.M_orig.F';
            M_cut.pathPairs=pathPairs;
            M_cut.Old2New=num2cell(1:length(obj.M_orig.V));
            M_cut.New2Old=1:length(obj.M_orig.V);
            obj.M_cut=M_cut;
        end
        function [vcol,ncol]=computeVCol(obj)
            
            TR=triangulation(obj.M_orig.F',obj.M_orig.V');
            vn = TR.vertexNormal();
            %             fn = TR.faceNormal();
            %             ncol=fn(:,1);
            cutVN= vn(obj.M_cut.New2Old,:);
            v=(1-cutVN(:,1))/2;
            %             v=v/4+0.75;
            ncol=v;%cmap(round(1+v*(length(cmap)-1)),:);
            %             ncol=repmat((cutVN(:,1)+1)/2,1,3);
            vcol=(cutVN+1)/2;
            vcol=rgb2hsv(vcol);
            vcol(:,2)=0.2;
            vcol(:,3)=1;
            vcol=hsv2rgb(vcol);
            
            vn=obj.M_orig.V';
            vn=vn-min(vn(:));
            vn=vn/max(vn(:));
            
            %             ncol=fn(:,1);
            %             cutVN= vn(obj.M_cut.New2Old,:);
            v=obj.M_cut.V;
            %vn=vn-min(vn(:));
            %vn=vn/max(vn(:));
            v = bsxfun(@minus,v,min(v));
            v = bsxfun(@rdivide,v,max(v));
            vcol=v*0.5+0.5;
            %             vcol=rgb2hsv(vcol);
            %             vcol(:,2)= vcol(:,2)*0.5;
            %             vcol(:,3)=vcol(:,3)/2+0.5;
            %             vcol=hsv2rgb(vcol);
            %             vcol = normr(vcol);
        end
        function computeV2A(obj)
            fprintf('==== Computing V2A matrix ===\n');
            tid=tic;
            [obj.V2A,obj.areas] = getFlatteningDiffCoefMatrix(obj.M_cut.V,obj.M_cut.T); % calculate map between 2d vertices to differentials
            fprintf('done computing, %f seconds\n',toc(tid));
        end
        
        
        
        function computeAs(obj)
            if isempty(obj.V2A)
                obj.computeV2A();
            end
            fprintf('==== Computing differentials ===\n');
            tid=tic;
            %obj.As = permute(reshape(obj.V2A*obj.flat_V,2,[],2),[1 3 2]);
            obj.As = reshape(obj.V2A*obj.flat_V(:),2,2,[]);
            fprintf('done computing, %f seconds\n',toc(tid));
        end
        function computeDistortion(obj,force)
            if nargin<2
                force=false;
            end
            if ~force && ~isempty(obj.dets)
                
                return;
            end
            if isempty(obj.As)||force
                obj.computeAs();
            end
            fprintf('==== Computing distortion etc. ===\n');
            tid=tic;
            As=obj.As;
            
            
            a = squeeze(As(1,1,:))';
            b = squeeze(As(1,2,:))';
            c = squeeze(As(2,1,:))';
            d = squeeze(As(2,2,:))';%the entries of A
            
            obj.dets=(a.*d-b.*c)';
            obj.frobenius=sqrt(a.^2+b.^2+c.^2+d.^2)';
            
            
            alpha = [a+d;b-c]/2; %2XM
            beta = [a-d;b+c]/2;
            alpha=sqrt(sum(alpha.^2));
            beta=sqrt(sum(beta.^2));
            obj.smax=(alpha+beta)';
            obj.smin=abs((alpha-beta)');
            
            obj.flipped=obj.dets<-1e-5;
            assert(~any(obj.flipped));
            fprintf('done computing, %f seconds\n',toc(tid));
        end
        function drawHist(obj,cutoff)
            prop=obj.smax./obj.smin;
            %             cutoff=min(cutoff,max(prop));
            m=mean(prop);
            
            prop(prop>cutoff)=cutoff;
            hist(prop,300);
            hold on
            yl=ylim;
            plot([m,m],yl,'r--','LineWidth',2)
            
            hold off
            ax=gca;
            ax.XTick = [];%[m,cutoff];
            geq=[];
            if max(prop)==cutoff
                geq='\geq';
            end
            xl=xlim;
            dist=(xl(2)-xl(1))/100;
            fsize=22;
            text(m+dist,yl(2),sprintf('%0.3f',m),'VerticalAlignment','top','fontsize',fsize);
            text(1,0,'1','VerticalAlignment','top','fontsize',fsize,'HorizontalAlignment','center');
            text(cutoff,0,[geq num2str(cutoff)],'VerticalAlignment','top','HorizontalAlignment','center','fontsize',fsize);
            %             ax.XTickLabel = {,[geq num2str(cutoff)]};
            ax.YTick=[];%ax.YTick(end);
            line(xl,[0 0 ],'color','k');
            line([1 1],yl,'color','k');
            %ax.YTickLabel=ax.YTickLabel(end);
            axis off
        end
        function visualizeBasicPolygon(obj,patchcolor)
            pathEnds=[];
            for i=1:length(obj.M_cut.pathPairs)
                pathEnds=[pathEnds obj.M_cut.pathPairs{i}([1 end],:)];
            end
            pathEnds=unique(pathEnds);
            TR=triangulation(obj.M_cut.T,obj.flat_V);
            all_binds = TR.freeBoundary();
            assert(all(all_binds(:,2)==all_binds([2:end,1],1)));
            all_binds=all_binds(:,1);
            %                 ind=find(all_binds==startP);
            %                 all_binds=all_binds([ind:end,1:ind-1]);
            p=find(ismember(all_binds,pathEnds));
            p=all_binds(p);
            if nargin>1
            fill(obj.flat_V(p,1),obj.flat_V(p,2),patchcolor);
            end
            hold on;
            line_colors=obj.cut_colors;
            
            
            
            LINE_WIDTH=6;
            curY=obj.flat_V;
            for i=1:length(obj.M_cut.pathPairs)
                p=obj.M_cut.pathPairs{i};
                %scatter(curY(p(1,1),1),curY(p(1,1),2),80,'filled');
                c=line_colors{i};%hsv2rgb([i/length(path_pairs),1,1]);
                
                
                
                line(curY(p([1 end],1),1),curY(p([1 end],1),2),'linewidth',LINE_WIDTH,'Color',c);
                line(curY(p([1 end],2),1),curY(p([1 end],2),2),'linewidth',LINE_WIDTH,'Color',c);
                cwidth=0.05;
                
                
                
            end
            % caxis([1 c_lim])
            axis equal
            %saveas(h,sprintf('%d.png',i),'png');
            for i=1:length(obj.inds)
                curInds=obj.M_cut.Old2New{obj.inds(i)};
                for j=1:length(curInds);
                    p=curInds(j);
                    %scatter(curY(p(1,1),1),curY(p(1,1),2),80,'filled');
                    %     if size(curY,2)==2
                    
                    
                    
                    %                     drawn_landmarks(end+1:end+2)=p(end,:);
                    cp=curY(p,:);
                    %             rectangle('position',[cp-cwidth/2 cwidth cwidth],'curvature',1,'facecolor',landmarks_colors{cur_landmark_color},'edgecolor','none');
%                     obj.drawLandmark(cp,obj.LM_colors{i});
                    scatter(cp(1),cp(2),80,obj.LM_colors{i},'filled');
                end
                
                
                
            end
        end
        function  visualize( obj,coloring,dim,clims)
            hold on;
            if nargin<3||dim==2
                curY=obj.flat_V;
            else
                curY=obj.M_cut.V;
            end
%             curY=curY-min(curY(:));
%             curY=curY/max(curY(:));
            T=obj.flat_T;
            if nargin==1 || ~strcmp(coloring,'none')
            p=patch('Faces',T,'Vertices',curY,'EdgeColor','k','linewidth',0.1,'FaceColor','none');
            end
            if nargin==1 || strcmp(coloring,'ncol')
                [vcol,ncol]=obj.computeVCol();
                a=0.85;
                vals=(300:999)';
                %                 cmap=kron(vals/max(vals),[0.95 0.99 1])+kron(1-(vals)/max(vals),[0.25 0.2 0.25 ]);
                cmap=kron(vals/max(vals),[1 1 1])+kron(1-(vals)/max(vals),[0 0 0 ]);
                shadecol=cmap(round(ncol*(length(cmap)-1)+1),:);
                col=a*shadecol+(1-a)*vcol;
                ncol=ncol*0.666+0.333;
                col=vcol.*[ncol ncol ncol];
                set(p,'FaceVertexCData',col,'FaceColor','interp','edgecolor','none');%,'linewidth',0.1,'edgealpha',0.05);
                
                %                 colormap(cmap);
            elseif strcmp(coloring,'none')
                %do nothing
            elseif strcmp(coloring,'vcol')
                
                vcol=obj.computeVCol();
                set(p,'FaceVertexCData',vcol,'FaceColor','interp','edgecolor','k','linewidth',0.1);%,'edgealpha',0.1);
            elseif strcmp(coloring,'dirichlet')||strcmp(coloring,'dirichletlog')
                obj.computeDistortion();
                if strcmp(coloring,'dirichlet')
                    fro=obj.frobenius;
                else
                    fro=log(obj.frobenius);
                end
                set(p,'FaceColor','flat','FaceVertexCData',fro,'CDataMapping','scaled','edgecolor','none');
                colormap(jet);
            elseif strcmp(coloring,'confdist')
                obj.computeDistortion();
                TR=triangulation(obj.M_orig.F',obj.M_orig.V');
                
                fn = TR.faceNormal();
                prop=obj.smax./obj.smin;
                MAP_RES=1000;
                cmap=distortion_colormap(MAP_RES);
                prop=max(prop,clims(1));
                prop=min(prop,clims(2));
                prop=prop-min(prop);
                
                prop=prop/max(prop);
                prop=min(prop,1);
                fn=(fn(:,1)+1+2);
                fn=fn/max(fn);
                C1=cmap(round(prop*(length(cmap)-1)+1),:).*[fn fn fn];;
                set(p,'FaceColor','flat','FaceVertexCData',C1,'CDataMapping','scaled','edgecolor','none');
            else   %                 colormap(distortion_colormap());
                                set(p,'FaceColor',coloring,'edgecolor','k','linewidth',0.1);
            end
                                            set(p,'edgecolor','k','linewidth',0.2,'edgealpha',0.1);

            patch('Faces',T(obj.flipped==1,:),'Vertices',curY,'FaceColor','yellow','linewidth',10);
            
            %     patch(curY(T(bd.flipped==1,:),1),curY(T(bd.flipped==1,:),2),'y',parent{:})
            
            % set(gcf,'units','normalized','outerposition',[0 0 1 1])
            line_colors=obj.cut_colors;{[0 0.4 1],[0.7 0 0],[0 1 0]};
            
            
            
            LINE_WIDTH=3;
            
            for i=1:length(obj.M_cut.pathPairs)
                p=obj.M_cut.pathPairs{i};
                %scatter(curY(p(1,1),1),curY(p(1,1),2),80,'filled');
                c=line_colors{i};%hsv2rgb([i/length(path_pairs),1,1]);
                
                if size(curY,2)==2
                    
                    line(curY(p(:,1),1),curY(p(:,1),2),'linewidth',LINE_WIDTH,'Color',c);
                    line(curY(p(:,2),1),curY(p(:,2),2),'linewidth',LINE_WIDTH,'Color',c);
                    cwidth=0.05;
                    
                    
                else
                    
                    line(curY(p(:,1),1),curY(p(:,1),2),curY(p(:,1),3),'linewidth',3,'Color',c);
                    line(curY(p(:,2),1),curY(p(:,2),2),curY(p(:,1),3),'linewidth',3,'Color',c);
                end
            end
            % caxis([1 c_lim])
            axis equal
            %saveas(h,sprintf('%d.png',i),'png');
            for i=1:length(obj.inds)
                curInds=obj.M_cut.Old2New{obj.inds(i)};
                for j=1:length(curInds);
                    p=curInds(j);
                    %scatter(curY(p(1,1),1),curY(p(1,1),2),80,'filled');
                    %     if size(curY,2)==2
                    
                    
                    
                    %                     drawn_landmarks(end+1:end+2)=p(end,:);
                    cp=curY(p,:);
                    %             rectangle('position',[cp-cwidth/2 cwidth cwidth],'curvature',1,'facecolor',landmarks_colors{cur_landmark_color},'edgecolor','none');
                    obj.drawLandmark(cp,obj.LM_colors{i});
                end
                
                
                
            end
            hold off
            
        end
        function drawLandmark(obj,cp,col)
            if size(cp,2)==2
                cwidth=0.1;
                rectangle('position',[cp-cwidth/2 cwidth cwidth],'curvature',1,'facecolor',col,'edgecolor','none');
            else
                
                %     scatter3(cp(1),cp(2),cp(3),300,col,'filled');
                [x,y,z] = sphere;
                x=x/40;
                y=y/40;
                z=z/40;
                
                
                
                surf(x+cp(1),y+cp(2),z+cp(3),'edgecolor','none','facecolor',col) % centered at (3,-2,0)
                
            end
        end
        function fixToAxis(obj)
            p=obj.flat_V(obj.M_cut.Old2New{obj.inds(3)},:);
            p1=p(1,:);
            p2=p(2,:);
            [ M,t ] = computeSimilarity( p1,p2,[1 0],[-1 0] );
            obj.flat_V=bsxfun(@plus,obj.flat_V*M',t');
        end
        function fixToGrid(obj)
            obj.fixToAxis();
            obj.flat_V=obj.flat_V*10;
            ind1=obj.M_cut.Old2New{obj.inds(3)}(1);
            ind2=obj.M_cut.Old2New{obj.inds(3)}(2);
            
            X=obj.flat_V;
            if X(ind1,1)>X(ind2,1)
                ind3=ind1;
                ind1=ind2;
                ind2=ind3;
            end
            X(:,1)=X(:,1)-X(ind1,1);
            m=X(ind2,1);
            X=2*round(m)*X/m;
            
            
            ind3=obj.M_cut.Old2New{obj.inds(4)}(1);
            ind4=obj.M_cut.Old2New{obj.inds(4)}(2);
            p1=X(ind3,:);
            p2=X(ind4,:);
            p3=X(ind2,:);
            q1=round(p1);
            q2=round(p2);
            q3=p3;
            [ A,t ] = affineTransFrom3Points( p1,p2,p3,q1,q2,q3 );
            X=bsxfun(@plus,X*A',t');
            obj.flat_V=X;
            
% %             ind1=obj.M_cut.Old2New{obj.inds(1)};
% %             ind2=obj.M_cut.Old2New{obj.inds(2)};
% ind1=obj.M_cut.Old2New{obj.inds(3)}(1);
% ind2=obj.M_cut.Old2New{obj.inds(3)}(2);
%             obj.fixToAxis();
% %             obj.flat_V=obj.flat_V*10;
%             X=obj.flat_V(:,1);
%             if X(ind1,1)>X(ind2,1)
%                 ind3=ind1;
%                 ind1=ind2;
%                 ind2=ind3;
%             end
%             X=X-X(ind1);
%             m=X(ind2);
%             X=round(m)*X/m;
%             obj.flat_V(:,1)=X;
%             
%             
%             ind1=obj.M_cut.Old2New{obj.inds(1)};
%             ind2=obj.M_cut.Old2New{obj.inds(2)};
%             
%             X=obj.flat_V(:,2);
%             if X(ind1,1)>X(ind2,1)
%                 ind3=ind1;
%                 ind1=ind2;
%                 ind2=ind3;
%             end
%              X=X-X(ind1);
%             m=X(ind2);
%             X=round(m)*X/m;
%             obj.flat_V(:,2)=X;
%              obj.flat_V= obj.flat_V*2;
        end
        
        %         function UV=uncutUV(obj)
        %             UV=obj.flat_V(cellfun(@(C)C(1),flattener.M_cut.Old2New),:);
        %         end
        
    end
    
end


