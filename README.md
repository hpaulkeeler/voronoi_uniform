# voronoi_uniform
Uniformly places a random point on all *bounded* cells of a two-dimensional Voronoi/Dirichlet tesselation. 

If you use this code in a research publication, please cite this repository.

The MATLAB and Python code do essentially the same thing, with the key functions respectively being located in funVoronoiUniform.m and funVoronoiUniform.py. More details are found in the comments in  the files.

The uniform placement step is done by first dividing each (bounded) Voronoi cell (ie an irregular polygon) with, say, m sides into m scalene triangles. The i th triangle is then randomly chosen based on the ratio of areas. A point is then uniformly placed on the i th triangle (via eq. 1 in [1]). The placement step is repeated for all bounded Voronoi cells.

To test the functions, run VoronoiUniformTest.m or VoronoiUniformTest.py (with the previous files in the same directory/folder). For a given point pattern (for example, a realization of a Poisson point process), these files run the function funVoronoiUniform repeatedly, placing a single uniform points on each cell over many simulations. Then the averages are taken and compared to the centroids of each (bounded) cell. For large number of simulations, exact and empirical centroids should agree.

Author: H.Paul Keeler, 2019 
[1] Osada, R., Funkhouser, T., Chazelle, B. and Dobkin. D., "Shape distributions", ACM Transactions on Graphics, vol 21, issue 4,
 2002
