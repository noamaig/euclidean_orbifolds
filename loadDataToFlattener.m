function flattener = loadDataToFlattener(fname)
a=load(fname);
[V,T]=read_off(a.mesh_name);
V=V/sqrt(mean(var(V')));
M_orig=[];
M_orig.V=V;
M_orig.F=T;
inds=[length(a.inds)-1:-1:1 length(a.inds)];
flattener=Flattener(M_orig,a.inds(inds),a.singularities);
end

