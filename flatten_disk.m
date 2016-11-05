function [V_flat,boundary_segments,flattener]=flatten_disk(V,T,inds,varargin)
%flatten a disk-topology mesh to a plane - to a triangle domain or square,
%depending on the number of prescribed cones
%input 
% V,T - mesh. 
% inds - the indices of the cones
%output 
% V_flat - the embedded vertices positions
% boundary_segments - the indiceds of boundary vertices, grouped into a
% cellarry where A{i} is the array of boundary vertices of i'th edge of the
% domain (tri\square)
% flattener - the Flattener object
p = inputParser;
p.addParameter('verbose',true,@islogical)




if length(inds)==3
    tri_or_square=true;
    
elseif length(inds)==4
    tri_or_square=false;
else 
    error('disk orbifolds need exactly 3 or 4 cones');
end
    

M_orig=[];
M_orig.V=V';
M_orig.F=T';
p.parse(varargin{:});

flattener=Flattener(M_orig,inds,[]);
flattener.verbose=p.Results.verbose;
if tri_or_square
    orbifold_type='freetri';
else
    orbifold_type='freesquare';
end
flattener.flatten(orbifold_type);

V_flat=flattener.flat_V;
boundary_segments=cellfun(@(X)X(:,1),flattener.M_cut.pathPairs,'UniformOutput',false);
end