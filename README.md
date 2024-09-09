# coho

GitLab repository: https://gitlab.tue.nl/20204732/ijnme-computational-homogenization-of-poroelastic-media

dataset DOI: https://data.4tu.nl/datasets/f422cac5-b959-445d-9e2b-ad4067ee9580

### Renan Liupekevicius Carnielli [r.liupekevicius.carnielli@tue.nl](mailto::r.liupekevicius.carnielli@tue.nl)

Data Availability of the IJNME article "Transient computational homogenization of heterogeneous poroelastic media with local resonances"

Author: Renan Liupekevicius Carnielli, TU Eindhoven, April 2024

ASSOCIATED PEER-REVIEWED PUBLICATION
https://onlinelibrary.wiley.com/doi/10.1002/nme.7505

## comsol.zip folder description (DOWNLOAD AVAILABLE https://data.4tu.nl/datasets/f422cac5-b959-445d-9e2b-ad4067ee9580)

This zip folder contains the comsol models for the direct numerical simulations (DNS) and the homogenised macro-scale model (HOM) shown in the article. 

-*dispersion DNS* is the comsol model that computes the Bloch Analysis (figure 4) and the local resonance modes (figure 3);

-*transmission DNS* is the comsol model that computes the DNS of the transmission problem, figure 5(c) dashed curve.

-*transmission homogenized* is the comsol model that computes the homogenised transmission problem, figure 5(c) solid curve.


## Note

1- run `start.m` to include the the path to the folders *fun*, *tensorlab* and *meshes*.

2- export comsol 5.4 mesh file `.mphtxt` of your favorite geometry to the folder *meshes*, or use the mesh file provided.

3- FATAL: remember to rename the exported comsol mesh `.mphtxt` to `.txt`.

4- the files `computedworkspace_hexagon_loc.mat` is the output workspaces of `main.m` for the unit cell, figure 2 of the article.
## `main.m` file description

1- `main.m` computes the homogenized material coefficients of equations (47) and equation (48).

2-  The input for `main.m` are: the mesh file in the folder *meshes* with the list of tags of the corresponding comsol model, and the material properties of the saturating fluid and the solid phases which can be defined within the script `main.m`.

3-  The appropriate mesh is created by exporting a text file `.mphtxt` from comsol 5.4, see 'Note' 2nd and 3rd item.


## *fun* folder
Developed functions to assist the computations in `main.m`.

## *tensorlab* folder
Collection of functions and classes that define a tensor object.

## *meshes* folder
It contains the mesh files of the unit cell shown in the article or any other unit cell you may want to homogenise. 
