classdef PosConstraints<handle
    %object to aggregate positional constraints and generate the
    %corresponding matrix and RHS
    properties
        A; %the constraint matrix
        b=[]; %the RHS
        
    end
    
    methods
        function obj=PosConstraints(nvars)
            %nvars - how many vertices are there in the mesh (to set width
            %of the resulting matrix)
            obj.A=sparse(0,nvars*2);
        end
        function addConstraint(obj,ind,w,rhs)
            % adds a positional constraint on a vertex x_ind, 
            % so that x_ind*w=rhs 
            assert(length(rhs)==2);
            obj.A(end+1,ind*2-1)=w;
            obj.A(end+1,ind*2)=w;
            
            obj.b=[obj.b;rhs];
        end
        function newConstraint(obj)
            %add a new line to the constraints matrix
            obj.A(end+1,1)=0;
            obj.b(end+1)=0;
        end
        function addLineConstraint(obj,ind,n,offset)
            %add constraint for x_ind to lie on an infinite line, according
            %to a normal and offset, so that <x_ind,n>=offset
            obj.A(end+1,ind*2-1:ind*2)=n;
            
            obj.b=[obj.b;offset];
        end
        function addTransConstraints(obj,sinds,tinds,T)
            %add constraints so that T*x_sinds=x_tinds where T is a 2X2
            %matrix, and the Transformation T is modified to affine from
            %linear by requiring that T*x_si-x_ti=T*x_s1-x_t1
            
            %iterate over all vertices except first vertex
            for ind=2:length(sinds)
                %iterate on both rows ot the 2X2 matrix T
                for vert_ind=1:2
                    %building up the equations by summing up the terms
                    %
                    % <T(vert_ind,:),x_si>
                    obj.A(end+1,sinds(ind)*2+[-1,0])=T(vert_ind,:);
                    % -<T(vert_ind,:),x_s1> 
                    obj.A(end,sinds(1)*2+[-1,0])=obj.A(end,sinds(1)*2+[-1,0])-T(vert_ind,:);
                    %  -x_ti
                    obj.A(end,tinds(ind)*2+vert_ind-2)=obj.A(end,tinds(ind)*2+vert_ind-2)-1;
                    % +x_t1
                    obj.A(end,tinds(1)*2+vert_ind-2)=obj.A(end,tinds(1)*2+vert_ind-2)+1;
                    %left hand side is zero
                    obj.b=[obj.b;0];
                end
            end
        end
    end
    
end

