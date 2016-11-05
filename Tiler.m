classdef Tiler<handle
    %obnject which generates tiling of a specific orbifold
    
    
    properties
        V_flat;%the image of the vertices
        T;%triangulation
        pathPairs; %cell array of corresponding boundary vertices (see CutMesh.m)
        trans={[1 0;0 1;0 0]};%cell array holding all the transformations generating the tiling,
        %as a cell array where A{i} is a 6x2 matrix of the form [A;v]
        %which represents the affine trans x -> Ax+v
        stack={};%transformations to consider
    end
    
    methods
        function obj=Tiler(V,T,pathPairs)
            %input:
%             V_flat;%the image of the vertices
%             T;%triangulation
%             pathPairs; %cell array of corresponding boundary vertices (see CutMesh.m)
            obj.V_flat=V;
            obj.T=T;
            obj.pathPairs=pathPairs;
            
        end
        function push(obj,A,depth)
            %inserting a transformation A to the stack
            %depth is the recursion depth
            obj.stack{end+1}=[];
            obj.stack{end}.A=A;
            obj.stack{end}.depth=depth;
        end
        function tile(obj,depth)
            %generate the tiling, with a maximum given recursion depth
            obj.push(obj.trans{1},depth);
            while(true)
                item=obj.pop();
                if isempty(item)
                    break;
                end
            obj.addTrans(item.A,item.depth);
            end
        end
        function item=pop(obj)
            %get next transforatmion to check
            if isempty(obj.stack)
                item=[];
                return;
            end
            item=obj.stack{end};
            obj.stack(end)=[];
        end
        function addTrans(obj,A,depth)
            assert(depth>0);
            curV=obj.V_flat*A([1 2],:)'+repmat(A(3,:),length(obj.V_flat),1);
            for i=1:length(obj.pathPairs)
                for j=1:2
                    
                     p1=obj.V_flat(obj.pathPairs{i}(1,j),:);
                     p2=obj.V_flat(obj.pathPairs{i}(end,j),:);
                     q1=curV(obj.pathPairs{i}(1,3-j),:);
                     q2=curV(obj.pathPairs{i}(end,3-j),:);
                     [ R,t ] = computeSimilarity( p1,p2,q1,q2 );
                     t=t';
                     A=[R;t];
                     diffs=cellfun(@(x)norm(A-x,'fro'),obj.trans);
                     if ~isempty(diffs) && min(diffs)<1e-4
                         continue;
                     end
                     obj.trans{end+1}=A;
%                     
                     if depth>1
                        obj.push(A,depth-1);
                     end
                end
            end
        end
    end
    
end

