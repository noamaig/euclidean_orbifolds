classdef CutMesh < handle
    %represents a sphere mesh which has been cut open
    
    properties
        V;%vertices of the cut mesh
        T;%triangles of the cut mesh
        pathPairs; %cellarray, where A{i} corresponds to the i'th cut, and 
        % is a m-on-2 matrix where A{i}(m,:) are the indices of the two 
        % corresponding vertices which were duplicated from the m'th vertex
        % on the UNCUT mesh.
        
        cutIndsToUncutInds;% array where A(i) is the index of the original 
        %of vertex i on the uncut mesh
        uncutIndsToCutInds;% cellarry where A{i} contains an array of all 
        %duplicates on the cut mesh, of the vertex i on the uncut mesh
    end
    
    methods
        function obj=CutMesh(V,T,pathPairs,cutIndsToUncutInds,uncutIndsToCutInds)
%             V;%vertices of the cut mesh
%         T;%triangles of the cut mesh
%         pathPairs; %cellarray, where A{i} corresponds to the i'th cut, and 
%         % is a m-on-2 matrix where A{i}(m,:) are the indices of the two 
%         % corresponding vertices which were duplicated from the m'th vertex
%         % on the UNCUT mesh.
%         
%         cutIndsToUncutInds;% array where A(i) is the index of the original 
%         %of vertex i on the uncut mesh
%         uncutIndsToCutInds;% cellarry where A{i} contains an array of all 
%         %duplicates on the cut mesh, of the vertex i on the uncut mesh
            obj.V=V;
            obj.T=T;
            obj.pathPairs=pathPairs;
            obj.cutIndsToUncutInds=cutIndsToUncutInds;
            obj.uncutIndsToCutInds=uncutIndsToCutInds;
        end
    
    end
    
end

