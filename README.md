# voronoi_uniform
Uniformly places a random point in all *bounded* cells of a two-dimensional Voronoi/Dirichlet tesselation of some given point pattern, such as a realization of a (random) point process. 

If you use this code in a research publication, please cite this repository:

@misc{keeler2019voronoi,
author = {Keeler, H.P.},
title = {Uniform placement of random points in bounded Voronoi cells},
year = {2019},
publisher = {GitHub},
journal = {GitHub repository},
howpublished = {\url{https://github.com/hpkeeler/voronoi_uniform}},
commit = {bdb3a6239476c6cebc8ef63953a3215027e65d66}
}


The MATLAB and Python code do essentially the same thing, with the key functions respectively being located in funVoronoiUniform.m and funVoronoiUniform.py. More details are found in the comments in  the files. To create a Voronoi/Dirichlet tesselation, either the MATLAB[1] or SciPy[2] function needs to be run, depending on which language is being used.

The random (uniform) placement step is done by first dividing each (bounded) Voronoi cell (ie an irregular polygon) with, say, m sides into m scalene triangles. The i th triangle is then randomly chosen based on the ratio of areas. A point is then uniformly placed on the i-th triangle (via eq. 1 in [3]). The random (uniform) placement step is repeated for all *bounded* Voronoi cells.

To test the functions, run VoronoiUniformTest.m or VoronoiUniformTest.py (with the aforementioned files in the same directory/folder). For a given point pattern (for example, a realization of a Poisson point process), these files run the function funVoronoiUniform repeatedly over number of simulations. For each simulation, a single random  point is uniformly placed in each *bounded* cell. Then the averages  of the Cartesian components of the uniformly placed points are then taken.  As the number of simulations (of placing single points) increases, these averages should converge to the centroids (or geometric centres) of all the *bounded* Voronoi cells.

Author: H.Paul Keeler, 2019 
[1] http://www.mathworks.com.au/help/matlab/ref/voronoin.html
[2]  http://scipy.github.io/devdocs/generated/scipy.spatial.Voronoi.html
[3] Osada, R., Funkhouser, T., Chazelle, B. and Dobkin. D., "Shape distributions", ACM Transactions on Graphics, vol 21, issue 4,
 2002
