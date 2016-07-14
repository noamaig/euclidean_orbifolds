classdef SingleMeshCutter < handle
    %COMPATIBLEMESHSEGMENTER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        M1
        I1
        C1
        
        Method='Tree';
        Cutter1
    end
    properties(SetObservable)
        Path1
        Path2
        PathInd1
        PathInd2
        FixedPairs
        Yparam1
        UseTutte=0;
        UseOriginalMetric=0;
    end
    
    methods
        function obj=SingleMeshCutter(M1,I1)
            obj.Yparam1=1;
            obj.M1=M1;
            obj.Cutter1=MeshCutter(M1);
            obj.I1=I1;
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
            M1=obj.GatherResults;
            Segs={M1};
        end
        function [M1]=GatherResults(obj)
            Cutter1=obj.Cutter1;
           
            M1=struct('V',Cutter1.CutM.V,'F',Cutter1.CutM.F,...
                'Old2New',{Cutter1.Old2New},'New2Old',{Cutter1.New2Old},...
                'I',obj.I1,'Paths',{obj.PathInd1},'Segments',obj.C1);
           
        end
    end
    methods(Access=private)
        function CutMeshes(obj)
            if ~obj.UseTutte && obj.UseOriginalMetric
                error('UseOriginalMetric can only be used with UseTutte');
            end
            
            
            FixedPairs=obj.FixedPairs;
            
            obj.Path1=cell(numel(obj.I1));
           
            obj.PathInd1=cell(numel(obj.I1));
           
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
            
            N2O=1:numel(obj.I1);
            % TODO: To save time, don't look at pairs that we already know
            % are disconnected
            Yflag=1;
            TutteFlag=1;
            while true
                % get Adj matrix of the cut meshes.
                [A1]=computeAdjacencies;
                N=numel(In1);
                D=inf(N);
                for k=unique(I(:,1))'
                    i=k; j=I(I(:,1)==k,2)';
                    % We need to remove all vertices except i,j
                    J=setdiff(1:N,[i,j]);
                    [A1t]=RemoveIndFromAdjacencies(A1,J);
                    D(i,j) = graphshortestpath(A1t,In1(i),In1(j),'Directed', false);
                
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
                    
                    obj.Cutter1.CutM.V=obj.M1.V(:,obj.Cutter1.New2Old);
                    
                    [A1]=computeAdjacencies;
                    obj.Cutter1.CutM.V=tempV1;
                    
                end
                [A1t]=RemoveIndFromAdjacencies(A1,J);
                [~,path1,~]=graphshortestpath(A1t,In1(i),In1(j),'Directed', false);
                
                i=N2O(i); j=N2O(j);
                obj.UpdatePathsLists(i,j,path1) % New paths between vertex i and j
                
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
                [In1,I,N2O]=obj.AddCut(path1,CandidatePairs);
                if obj.UseTutte && TutteFlag
                    TutteIt(obj.Cutter1.CutM);
                    
                    TutteFlag=0;
                end
            end
            if obj.UseTutte
                obj.Cutter1.CutM.V=obj.M1.V(:,obj.Cutter1.New2Old);
                
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
            function [A1t]=RemoveIndFromAdjacencies(A1,J)
                A1t=A1;
                
                A1t(In1(J),:)=0;
                A1t(:,In1(J))=0;
                
            end
            function [A1]=computeAdjacencies
                CM1=obj.Cutter1.CutM;
               
                A1=CM1.computeAdjacencyMatrix('distance');
            
                B=CM1.findBoundaries;
                B=setdiff(B,In1);
                A1(B,:)=0;
               
               
            end
        end
        function UpdatePathsLists(obj,i,j,path1,path2)
            % Update the arrays of paths
            obj.Path1{i,j}=obj.Cutter1.CutM.V(:,path1);
            
            obj.PathInd1{i,j}=obj.Cutter1.New2Old(path1);
            
            notify(obj,'PathsChanged');
        end
        function FindSegments(obj)
            % first we give each segment a number. each segment is a
            % connected component
            Cutter1=obj.Cutter1;
          
            O2N1=Cutter1.Old2New;
           
            N2O1=Cutter1.New2Old;
           
            A1=Cutter1.CutM.computeAdjacencyMatrix;
            
            [S1,C1]=graphconncomp(A1,'Directed',false);
           
           
            obj.C1=C1;
            
        end
        
        function [I1,I,N2O]=AddCut(obj,P1,I)
            obj.Cutter1.Cut(P1);
            O2N1=obj.Cutter1.Old2New(obj.I1);
            dups=cellfun(@numel,O2N1);
            N2O=[];
            for i=1:length(dups)
                N2O=[N2O,zeros(1,dups(i))+i];
            end
            I1=cell2mat(O2N1);
            IO2N=mat2cell(1:sum(dups),1,dups);
            Inew=zeros(0,2);
            for i=1:size(I,1)
                l=IO2N{I(i,1)}';
                r=IO2N{I(i,2)}';
                Inew=[Inew;[kron(l,ones(numel(r),1)),kron(ones(numel(l),1),r)]];
            end
            I=Inew;
        end
       
    end
    methods(Static)
        
    end
    events
        PathsChanged
    end
end

