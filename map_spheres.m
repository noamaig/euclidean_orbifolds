function [ BC1to2,BC2to1,flattener1,flattener2 ] = map_spheres( V1,T1,V2,T2,cones1,cones2,orbifold_type,varargin )
% input 
% V_i,T_i,cones_i - the vertices triangles and cone indices of mesh i
% orbifold_type - the type of orbifold to use
% output  - BCitoj - the barycentric coordinates matrix taking the vertices
% of mesh i to the interior of faces of mesh j
% flattner_i - the flattener of mesh i
p = inputParser;
p.addParameter('verbose',true,@islogical)
p.parse(varargin{:});
[V_flat1,cutMesh1,flattener1]=flatten_sphere(V1,T1,cones1,orbifold_type,varargin{:});
[V_flat2,cutMesh2,flattener2]=flatten_sphere(V2,T2,cones2,orbifold_type,varargin{:});
if p.Results.verbose
    fprintf('Extracting surface mapping 1->2 from orbifold embeddings: ');
    tid=tic;
end
BC1to2=compute_map_from_sphere_embeddings( V_flat1,V_flat2,cutMesh1,cutMesh2 );
if p.Results.verbose
    toc(tid);
    
end
if p.Results.verbose
    fprintf('Extracting surface mapping 2->1 from orbifold embeddings: ');
    tid=tic;
end
BC2to1=compute_map_from_sphere_embeddings( V_flat2,V_flat1,cutMesh2,cutMesh1 );
if p.Results.verbose
    toc(tid);

    
end
end

