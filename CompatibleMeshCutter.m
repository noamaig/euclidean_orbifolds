classdef CompatibleMeshCutter < handle
    %COMPATIBLEMESHSEGMENTER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        M1
        I1
        C1
        M2
        I2
        C2
        Method='Tree'
        Cutter1
        Cutter2
    end
    properties(SetObservable)
        Path1
        Path2
        PathInd1
        PathInd2
        FixedPairs
        Yparam1
        UseTutte=1;
        UseOriginalMetric=0;
    end
    
    methods
        function obj=CompatibleMeshCutter
            obj.Yparam1=1;
        end
        function setMesh1(obj,M1,I1)
            obj.M1=M1;
            obj.Cutter1=MeshCutter(M1);
            obj.I1=I1;
        end
        function setMesh2(obj,M2,I2)
            obj.M2=M2;
            obj.Cutter2=MeshCutter(M2);
            obj.I2=I2;
        end
        function setMethod(obj,Method)
            obj.Method=Method;
        end
        function setFixedPairs(obj,FixedPairs)
            obj.FixedPairs=FixedPairs;
        end
        function setYaronHeuristic1(obj,num)
            obj.Yparam1=num;
        end
        function setUseTutte(obj,flag)
            obj.UseTutte=flag;
        end
        function setUseOriginalMetric(obj,flag)
            obj.UseOriginalMetric=flag;
        end
        function Segs=Cut(obj)
            obj.CutMeshes;
            obj.FindSegments;
            [M1,M2]=obj.GatherResults;
            Segs={M1,M2};
        end
        function [M1,M2]=GatherResults(obj)
            Cutter1=obj.Cutter1;
            Cutter2=obj.Cutter2;
            M1=struct('V',Cutter1.CutM.V,'F',Cutter1.CutM.F,...
                'Old2New',{Cutter1.Old2New},'New2Old',{Cutter1.New2Old},...
                'I',obj.I1,'Paths',{obj.PathInd1},'Segments',obj.C1);
            M2=struct('V',Cutter2.CutM.V,'F',Cutter2.CutM.F,...
                'Old2New',{Cutter2.Old2New},'New2Old',{Cutter2.New2Old},...
                'I',obj.I2,'Paths',{obj.PathInd2},'Segments',obj.C2);
            pathPairs=convertFromCuttingToMesh(M1);
            M1.pathPairs=pathPairs;
            pathPairs=convertFromCuttingToMesh(M2);
            M2.pathPairs=pathPairs;
        end
    end
    methods(Access=private)
        function CutMeshes(obj)
            if ~obj.UseTutte && obj.UseOriginalMetric
                error('UseOriginalMetric can only be used with UseTutte');
            end
            assert(numel(obj.I1)==numel(obj.I2));
            
            FixedPairs=obj.FixedPairs;
            
            obj.Path1=cell(numel(obj.I1));
            obj.Path2=cell(numel(obj.I1));
            obj.PathInd1=cell(numel(obj.I1));
            obj.PathInd2=cell(numel(obj.I1));
            if isempty(FixedPairs)
                [i,j]=find(triu(ones(numel(obj.I1)),1));
                CandidatePairs=[i,j]; %list of yet to be connected pairs in the original mesh
            elseif ~iscell(FixedPairs)
                CandidatePairs=FixedPairs(1,:); %list of yet to be connected pairs in the original mesh
                FixedPairs(1,:)=[];
            else
                CandidatePairs=FixedPairs{1};
                FixedPairs=FixedPairs(2:end);
            end
            ConnectedPairs=zeros(0,2);
            % init
            I=CandidatePairs;
            In1=obj.I1;
            In2=obj.I2;
            N2O=1:numel(obj.I1);
            % TODO: To save time, don't look at pairs that we already know
            % are disconnected
            Yflag=1;
            TutteFlag=1;
            while true
                % get Adj matrix of the cut meshes.
                [A1,A2]=computeAdjacencies;
                N=numel(In1);
                D=inf(N);
                for k=unique(I(:,1))'
                    i=k; j=I(I(:,1)==k,2)';
                    % We need to remove all vertices except i,j
                    J=setdiff(1:N,[i,j]);
                    [A1t,A2t]=RemoveIndFromAdjacencies(A1,A2,J);
                    D(i,j) = graphshortestpath(A1t,In1(i),In1(j),'Directed', false);
                    D(i,j) = D(i,j)+graphshortestpath(A2t,In2(i),In2(j),'Directed', false);
                end

                if isinf(min(D(:)))
                    break;
                end
                [i,j]=find(D==min(D(:)));
                if (obj.Yparam1~=1) && Yflag
                    Dsort=sort(D(:));
                    val=Dsort(obj.Yparam1);
                    [i,j]=find(D==val);
                    Yflag=0;
                end

                J=setdiff(1:N,[i,j]);
                if obj.UseOriginalMetric
                    tempV1=obj.Cutter1.CutM.V;
                    tempV2=obj.Cutter2.CutM.V;
                    obj.Cutter1.CutM.V=obj.M1.V(:,obj.Cutter1.New2Old);
                    obj.Cutter2.CutM.V=obj.M2.V(:,obj.Cutter2.New2Old);
                    [A1,A2]=computeAdjacencies;
                    obj.Cutter1.CutM.V=tempV1;
                    obj.Cutter2.CutM.V=tempV2;
                end
                [A1t,A2t]=RemoveIndFromAdjacencies(A1,A2,J);
                [~,path1,~]=graphshortestpath(A1t,In1(i),In1(j),'Directed', false);
                [~,path2,~]=graphshortestpath(A2t,In2(i),In2(j),'Directed', false);
                i=N2O(i); j=N2O(j);
                obj.UpdatePathsLists(i,j,path1,path2) % New paths between vertex i and j
                
                % update CandidatePairs, depening on how which kind of cut we want
                % to do (patches,tree)
                if isempty(FixedPairs)
                    switch obj.Method
                        case 'Patches'
                            CandidatePairs=setdiff(CandidatePairs,[i,j],'rows'); % remove i,j from the list of landmark pairs
                            ConnectedPairs=[ConnectedPairs;[i,j]];
                        case 'Tree'
                            % Candidate pairs are all which have one connected
                            % vertex and one unconnected
                            ConnectedPairs=[ConnectedPairs;[i,j]];
                            Con=unique(ConnectedPairs);
                            Uncon=setdiff(1:numel(obj.I1),Con);
                            [p,q] = meshgrid(Con, Uncon);
                            CandidatePairs=[p(:),q(:)]
                    end
                else
                    ConnectedPairs=[ConnectedPairs;[i,j]];
                    CandidatePairs=FixedPairs(1,:); %list of yet to be connected pairs in the original mesh
                    FixedPairs(1,:)=[];
                end
                [In1,In2,I,N2O]=obj.AddCut(path1,path2,CandidatePairs);
                if obj.UseTutte && TutteFlag
                    TutteIt(obj.Cutter1.CutM);
                    TutteIt(obj.Cutter2.CutM);
                    TutteFlag=0;
                end
            end
            if obj.UseTutte
                obj.Cutter1.CutM.V=obj.M1.V(:,obj.Cutter1.New2Old);
                obj.Cutter2.CutM.V=obj.M2.V(:,obj.Cutter2.New2Old);
            end
            function TutteIt(M)
                T=triangulation(M.F',M.V');
                Bind=T.freeBoundary; Bind=Bind(:,1);
                NB=numel(Bind);
%                 A=M.computeAdjacencyMatrix; A=diag(sum(A,2))-A;
                A=M.computeLaplacian;
%                 C=sparse(1:NB,Bind,ones(NB,1),NB,M.Nv);
%                 B=zeros(M.Nv,2);
                DD=[cos(2*pi*(0:NB-1)/NB)',sin(2*pi*(0:NB-1)/NB)'];
%                 AA=[A'*A,C';C,zeros(NB)];
%                 xy=AA\[B;DD];
                Lind=zeros(1,M.Nv); Lind(Bind)=1;
                xy=A(~Lind,~Lind)\(-A(~Lind,Bind)*DD);
                M.V(3,:)=[];
                M.V(:,Bind)=DD';
                M.V(:,~Lind)=xy';
            end
            function [A1t,A2t]=RemoveIndFromAdjacencies(A1,A2,J)
                A1t=A1;
                A2t=A2;
                A1t(In1(J),:)=0;
                A1t(:,In1(J))=0;
                A2t(In2(J),:)=0;
                A2t(:,In2(J))=0;
            end
            function [A1,A2]=computeAdjacencies
                CM1=obj.Cutter1.CutM;
                CM2=obj.Cutter2.CutM;
                A1=CM1.computeAdjacencyMatrix('distance');
                A2=CM2.computeAdjacencyMatrix('distance');
                B=CM1.findBoundaries;
                B=setdiff(B,In1);
                A1(B,:)=0;
                A1(:,B)=0;
                B=CM2.findBoundaries;
                B=setdiff(B,In2);
                A2(B,:)=0;
                A2(:,B)=0;
            end
        end
        function UpdatePathsLists(obj,i,j,path1,path2)
            % Update the arrays of paths
            obj.Path1{i,j}=obj.Cutter1.CutM.V(:,path1);
            obj.Path2{i,j}=obj.Cutter2.CutM.V(:,path2);
            obj.PathInd1{i,j}=obj.Cutter1.New2Old(path1);
            obj.PathInd2{i,j}=obj.Cutter2.New2Old(path2);
            notify(obj,'PathsChanged');
        end
        function FindSegments(obj)
            % first we give each segment a number. each segment is a
            % connected component
            Cutter1=obj.Cutter1;
            Cutter2=obj.Cutter2;
            O2N1=Cutter1.Old2New;
            O2N2=Cutter2.Old2New;
            N2O1=Cutter1.New2Old;
            N2O2=Cutter2.New2Old;
            A1=Cutter1.CutM.computeAdjacencyMatrix;
            A2=Cutter2.CutM.computeAdjacencyMatrix;
            [S1,C1]=graphconncomp(A1,'Directed',false);
            [S2,C2]=graphconncomp(A2,'Directed',false);
            assert(S1==S2);
            % iterate over all split landmarks and change the connected
            % component number in the second mesh to match the first mesh.
            for i=1:numel(obj.I1)
                NI1=O2N1{obj.I1(i)};
                NI2=O2N2{obj.I2(i)};
                for j=1:numel(NI1)
                    CI1=C1(NI1(j));
                    CI2=C2(NI2(j));
                    if CI1==CI2
                        continue;
                    end
                    temp=(C2==CI1);
                    C2(C2==CI2)=CI1;
                    C2(temp)=CI2;
                end
            end
            obj.C1=C1;
            obj.C2=C2;
        end
        
        function [I1,I2,I,N2O]=AddCut(obj,P1,P2,I)
            obj.Cutter1.Cut(P1);
            obj.Cutter2.Cut(P2);
            O2N1=obj.Cutter1.Old2New(obj.I1);
            O2N2=obj.Cutter2.Old2New(obj.I2);
            dups=cellfun(@numel,O2N1);
            N2O=[];
            for i=1:length(dups)
                N2O=[N2O,zeros(1,dups(i))+i];
            end
            assert(all(dups==cellfun(@numel,O2N2)));
            I1=cell2mat(O2N1);
            I2=cell2mat(O2N2);
            IO2N=mat2cell(1:sum(dups),1,dups);
            Inew=zeros(0,2);
            for i=1:size(I,1)
                l=IO2N{I(i,1)}';
                r=IO2N{I(i,2)}';
                Inew=[Inew;[kron(l,ones(numel(r),1)),kron(ones(numel(l),1),r)]];
            end
            I=Inew;
        end
        function computePathsVoronoi(obj)
            error('Do not use')
            if Pr.Size~=Pl.Size
                error('Selected a different number of points');
            end
            Dl=zeros(Pl.Size,Ml.Nv);
            Dr=zeros(Pr.Size,Mr.Nv);
            Il=knnsearch(Ml.V',Pl.X');
            Ir=knnsearch(Mr.V',Pr.X');
            
            for i=1:Pr.Size
                Dl(i,:)=Ml.computeShortestPath(Il(i));
                Dr(i,:)=Mr.computeShortestPath(Ir(i));
            end
            [~,Vorl]=min(Dl);
            [~,Vorr]=min(Dr);
            FVl=unique(sort(Vorl(Ml.F)',2),'rows');
            FVr=unique(sort(Vorr(Mr.F)',2),'rows');
            Dell=zeros(Pl.Size);
            Delr=zeros(Pr.Size);
            for i=1:size(FVl,1)
                T=unique(FVl(i,:));
                if length(T)==1
                    continue;
                end
                Dell(T(1),T(2))=1;
                if length(T)==3
                    Dell(T(2),T(3))=1;
                    Dell(T(1),T(3))=1;
                end
            end
            for i=1:size(FVr,1)
                T=unique(FVr(i,:));
                if length(T)==1
                    continue;
                end
                Delr(T(1),T(2))=1;
                if length(T)==3
                    Delr(T(2),T(3))=1;
                    Delr(T(1),T(3))=1;
                end
            end
        end
    end
    methods(Static)
        
    end
    events
        PathsChanged
    end
end

