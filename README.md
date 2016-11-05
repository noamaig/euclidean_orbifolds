# Orbifold Tutte Embeddings, by Noam Aigerman & Yaron Lipman, presented at Siggraph Asia 2015.

Matlab code implementing the [Siggraph Asia 2015 paper, "Orbifold Tutte Embeddings"](http://www.wisdom.weizmann.ac.il/~noamaig/html/projects/orbifold/orbifold_highres.pdf).

An extension of Tutte's embedding to domains with sphere topology, and also to disks, where the boundary is free to move. 

###Main high-level functions:
- `flatten_sphere` - embed a sphere mesh to a sphere orbifold. 
- `flatten_disk` - embed a disk mesh to a disk orbifold. 
- `map_spheres` - compute a map between two sphere meshes using the orbifold embeddings.
- `map_disks` - compute a map between two disk meshes using the orbifold embeddings.

There are 17 Euclidean orbifolds, of which only the 4 sphere orbifolds, and 2 disk orbifolds are implemented here. The 11 other cases can be easily implemented, if one wishes, by extending the cutting mechanism (needs to produce a cut-graph that opens the mesh into a disk) and the appropriate boundary conditions.

###Example scripts:
- `script_embed_sphere` - map a spherical mesh to the Euclidean orbifold of type I.
- `script_embed_sphere_4_points` - map a spherical mesh to the Euclidean orbifold of type IV, first with the initial embedding, and then after modifying the configuration so as to minimize conformal distortion.
- `script_embed_disk` - map a mesh with disk topology to a disk orbifold.
- `script_surface_map_spheres` - use the sphere orbifold embedding to compute a seamless, conformal map between two human heads which interpolates 3 given correspondences, and then visualize the map with a cute animation of a morph between the two meshes.
- `script_surface_map_disk` - use the disk orbifold embedding to compute a conformal map between two disk-meshes which interpolates 3 designated points on the boundary.


The code is provided as-is for academic use only and without any guarantees. Please contact the author to report any bugs.
Written by [Noam Aigerman](http://www.wisdom.weizmann.ac.il/~noamaig/).

 