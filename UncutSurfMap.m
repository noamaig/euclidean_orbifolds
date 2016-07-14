classdef UncutSurfMap <handle
    % a mapping of the vertices of the uncut meshes
    
    properties
        M1;
        M2;
        barCoords1to2;
        barCoords2to1;
    end
    
    methods
        function obj=UncutSurfMap(M1_uncut,M2_uncut,M1_cut,M2_cut,lifter1to2,lifter2to1)
         
            % BC is for vertices on cut mesh. We care only about original mesh. So we can
            % take the BC of only one of the duplucated vertices for each of the old
            % vertices.
            obj.M1=M1_uncut;
            obj.M2=M2_uncut;
            O2N1=M1_cut.Old2New; 
            O2N2=M2_cut.Old2New;
            N2O1=M1_cut.New2Old; 
            N2O2=M2_cut.New2Old;
            O2N1=cellfun(@(x)x(end),O2N1); % get the first entry of each cell
            O2N2=cellfun(@(x)x(end),O2N2);
            BC1to2=lifter1to2.barCoordsMat;%for each point on 1, barcoords on 2
            BC2to1=lifter2to1.barCoordsMat;
            temp1=BC1to2(O2N1,:);
            temp2=BC2to1(O2N2,:);
            
            
            obj.barCoords1to2=sparse(size(M1_uncut.V,2),size(M2_uncut.V,2));
            obj.barCoords2to1=sparse(size(M2_uncut.V,2),size(M1_uncut.V,2));
            for i=1:size(temp1,1)
                ti=temp1(i,:);
                ind=find(ti);
                obj.barCoords1to2(i,N2O2(ind))=ti(ind);
            end
            for i=1:size(temp2,1)
                ti=temp2(i,:);
                ind=find(ti);
                obj.barCoords2to1(i,N2O1(ind))=ti(ind);
            end
        end
    end
end

