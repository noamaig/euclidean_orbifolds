classdef MeshCutter < handle
    %MESHCUTTER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        CutM
        Old2New
        New2Old
    end
    
    methods
        function obj=MeshCutter(M)
            obj.CutM=TriangleMesh(M);
            obj.Old2New=num2cell(1:M.Nv);
            obj.New2Old=(1:M.Nv)';
        end
        function Cut(obj,P)
            M=obj.CutM;
            V=M.V;
            F=M.F;
            NP=numel(P);
            % first check if one of the vertices of the path is on the
            % boundary
            B=M.findBoundaries;
            PB=ismember(P,B);
            % if any mid point of the path is on the boundary, throw an
            % exception. We can handle these cases in the future
            if any(PB(2:end-1))
                error('Path touches the boundary, and not in a good way');
            end
            
            % if the first or last vertices are not on the boundary, then
            % we don't have to do anything to them. We use that to set the
            % range of the loop
            flag=0; % first vertex is boundary flag
            range=1:NP;
            if ~PB(1)
                range(1)=[];
                flag=1;
            end
            if ~PB(end)
                range(end)=[];
            else
                lastRing=obj.CutM.findOrientedOneRing(P(end));
            end
            for i=range
                V(:,end+1)=V(:,P(i));
                obj.New2Old(end+1)=obj.New2Old(P(i));
                obj.Old2New{obj.New2Old(P(i))}(end+1)=size(V,2);
                ring=obj.CutM.findOrientedOneRing(P(i));
                % the current vertex of the ring is always on the boundary. The last
                % vertex we need is the one on the path.
                if flag && (i==2)
                    ring=circshift(ring,[1-find(ring==P(1)),0]);
                end
                if i<numel(P)
                    ring=ring(1:find(ring==P(i+1)));
                else
                    if find(ring==lastRing(1))==1
                        ring=ring(1:find(ring==P(i)));
                    else
                        ring=ring((find(ring==P(i))+1):end);
                    end
                end
                % Find face ring of vertex
                Fring=find(any(F==P(i)));
                % Relevant faces are those of Fring that have two vertices
                % from ring in them
                Fring=Fring(sum(ismember(F(:,Fring),ring),1)>=2);
                % now we need to take the faces in Fring, and change the
                % old vertex to the new one.
                FFring=F(:,Fring);
                FFring(FFring==P(i))=size(V,2);
                F(:,Fring)=FFring;
                obj.CutM=TriangleMesh('VF',V,F);
            end
        end
    end
end

