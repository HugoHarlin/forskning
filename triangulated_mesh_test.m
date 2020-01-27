%% Developer comments

% Experimental script testing a method of creating a cut mesh.
% the mesh is created using a square mesh that is intersected by some
% function, and the points in the square mesh that lie outside the
% intersection of the square mesh and the intersecting function are
% removed. The intersecting function is also "meshed" meaning that 
% it is represented as a series of points. Finally, the resulting
% collection of points are triangulated using the built ing delaunay triangulation
% function, in preparation for computation. 

% The reason ths is done is to be able to vary the bottom topography
% arbitrarily, with the aim to use hipsographic curves of lakes as the
% bottom topography. Because a finite volume method is used to solve the
% system of PDE:s, the mesh need not be regular and is therefore
% triangulated as stated previously. 

% There are some challenges involved: 
% 1 - remove meshgrid points outside the desired domain
% 2 - once the triangulation is complete, the border elements
%     needs to be located from the resulting mesh somehow.
% 3 - once this is done, the side of the border triangles that face
%     outwards need to be identified systematically (there are bc:s folks)
% 4 - Once the mesh is completed, the finite volume method needs to be
%     formulated for an arbirary triangulated mesh (normal vectors calculated etc.)
%     so that it can be applied to the created mesh with any bottom
%     topography.

% if this can be done, we can use any bottom topography desired with
% arbirary resolution of the mesh.  

%% Code
N = 10; % resolution of the square mesh grid

[X_sq,Y_sq] = meshgrid((0:1/(N-1):1),(0:1/(N-1):1));
x = 0:1/(N-1):1;
bottom_fn = 1-x.^2;
bottom_coord = [x;bottom_fn];





