# Orbifold Tutte Embeddings - Siggraph Asia 2015 paper by Noam Aigerman & Yaron Lipman. 

An extension of Tutte's embedding to domains with sphere topology and also to disks where the boundary is free to move.

Included below is a Matlab implementation of the algorithm

## Example scripts:
`script_embed_sphere` - map a simple sphere mesh to the Euclidean orbifold of type I.

`script_embed_disk` - map a mesh with disk topology to a disk orbifold.

`script_surface_map_spheres` - use the sphere orbifold embedding to compute a seamless, quasi-conformal map between two sphere surfaces which interpolates 4 given correspondences, and then generate a nice animation of a morph between the two:

`script_surface_map_disk` - use the disk orbifold embedding to compute conformal map between two disk meshes which interpolates 3 designated points on the boundary.

The two main high-level functions one needs to use to compute an embedding are `flatten_sphere` and `flatten_disk`.
