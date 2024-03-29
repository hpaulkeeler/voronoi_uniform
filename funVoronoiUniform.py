# uu,vv, indexBound=funVoronoiUniform(voronoiData)
# This code places uniformly points on *bounded* cells of a Voronoi
# tesselation (also called a Voronoi diagram or Dirichlet tesselation).
# This code uses a Voronoi tessellation based on an artibrary point
# two-dimensional pattern. The Voronoi tesselation is first found using the
# SciPy function[1], which is based on the Qhull project[2].
#
# INPUTS:
#
# voronoiData is data structure, created by the SciPy function Voronoi, that
# describes the Voronoi tesselation; see SciPy function Voronoi[1].
#
# OUTPUTS:
# uu and vv are vectors corresponding to the Cartesian coordinates of the
# uniformly placed points.
#
# indexBound is an index array for the bounded cells
#
# EXAMPLE: consider a point pattern described by 1-D arrays xx and yy,
# correspondong to the Cartesian coordinates. Then run the SciPy function
#
# voronoiData=Voronoi(np.stack((xx,yy), axis=1));
#
# Then the uniform points on bounded cells are obtained with:
#
# uu,vv, indexBound=funVoronoiUniform(voronoiData);
#
# All points (or cells) of the point process are numbered arbitrarily.
# In each *bounded* Voronoi cell a new point is uniformly placed.
#
# The placement step is done by first dividing each (bounded) Voronoi cell
# (ie an irregular polygon) with, say, m sides into m scalene triangles.
# The i th triangle is then randomly chosen based on the ratio of areas.
# A point is then uniformly placed on the i th triangle (via eq. 1 in [3]).
# The random placement step is repeated for all bounded Voronoi cells.
#
# Author: H.Paul Keeler, 2019.
# hpaulkeeler.com
# github.com/hpaulkeeler/voronoi_uniform
#
# References:
# [1] http://scipy.github.io/devdocs/generated/scipy.spatial.Voronoi.html
# [2] http://www.qhull.org/
# [3] Osada, R., Funkhouser, T., Chazelle, B. and Dobkin. D.,
# "Shape distributions", ACM Transactions on Graphics, vol 21, issue 4,
# 2002

import numpy as np  # NumPy package for arrays, random number generation, etc
from scipy.spatial import Voronoi  # for voronoi tessellation


def funVoronoiUniform(voronoiData):
    #retrieve points of underlying point pattern
    xx = voronoiData.points[:, 0]  # component
    yy = voronoiData.points[:, 1]  # component

    numbCells = len(xx)  # number of Voronoi cells (including unbounded ones)

    vertexAll = voronoiData.vertices  # retrieve x/y coordiantes of all vertices
    cellAll = voronoiData.regions  # may contain empty array/set

    indexP2C = voronoiData.point_region  # index mapping between cells and points

    ##initialize  arrays
    booleBound = np.zeros(numbCells, dtype=bool)
    uu = np.zeros(numbCells)
    vv = np.zeros(numbCells)

    ##Loop through for all Voronoi cells
    for ii in range(numbCells):
        booleBound[ii] = not(any(np.array(cellAll[indexP2C[ii]]) == -1))
        #checks if the Voronoi cell is bounded. if bounded, calculates its area
        #and assigns a single point uniformly in the Voronoi cell.

        ### START -- Randomly place a point in a Voronoi cell -- START###
        if booleBound[ii]:
            xx0 = xx[ii]
            yy0 = yy[ii]  # the (Poisson) point of the Voronoi cell

            #print(jj,xx0,yy0)
            #x/y values for current cell
            xxCell = vertexAll[cellAll[indexP2C[ii]], 0]
            yyCell = vertexAll[cellAll[indexP2C[ii]], 1]

            ### START -- Caclulate areas of triangles -- START###
            #calculate area of triangles of bounded cell (ie irregular polygon)
            numbTri = len(xxCell)  # number of triangles

            #create some indices for calculating triangle areas
            indexVertex = np.arange(numbTri+1)  # number vertices
            indexVertex[-1] = 0  # repeat first index (ie returns to the start)
            indexVertex1 = indexVertex[np.arange(
                numbTri)]  # first vertex index
            indexVertex2 = indexVertex[np.arange(
                numbTri)+1]  # second vertex index
            #calculate areas of triangles using shoelace formula
            areaTri = abs((xxCell[indexVertex1]-xx0)*(yyCell[indexVertex2]-yy0)
                          - (xxCell[indexVertex2]-xx0)*(yyCell[indexVertex1]-yy0))/2
            areaPoly = sum(areaTri)  # total area of cell/polygon
            ###END-- Caclulate areas of triangles -- END###

            ###START -- Randomly placing point -- START###
            ### place a point uniformaly in the (bounded) polygon
            #randomly choose the triangle (from the set that forms the polygon)
            cdfTri = np.cumsum(areaTri)/areaPoly  # create triangle CDF
            # use CDF to choose #
            indexTri = (np.random.rand() <= cdfTri).argmax()

            indexVertex1 = indexVertex[indexTri]  # first vertex index
            indexVertex2 = indexVertex[indexTri+1]  # second vertex index
            #for all triangles except the last one
            xxTri = [xx0, xxCell[indexVertex1], xxCell[indexVertex2]]
            yyTri = [yy0, yyCell[indexVertex1], yyCell[indexVertex2]]

            #create two uniform random variables on unit interval
            uniRand1 = np.random.rand()
            uniRand2 = np.random.rand()

            #x coordinate (via eq. 1 in [3])
            uu[ii] = (1-np.sqrt(uniRand1))*xxTri[0]\
                + np.sqrt(uniRand1)*(1-uniRand2)*xxTri[1]\
                + np.sqrt(uniRand1)*uniRand2*xxTri[2]
            #y coordinate (via eq. 1 in [3])
            vv[ii] = (1-np.sqrt(uniRand1))*yyTri[0]\
                + np.sqrt(uniRand1)*(1-uniRand2)*yyTri[1]\
                + np.sqrt(uniRand1)*uniRand2*yyTri[2]
            ###END -- Randomly placing point -- END###

    ### END -- Randomly place a point in a Voronoi cell -- END###

    indexBound = np.arange(numbCells)[booleBound]  # find bounded cells
    #remove unbounded cells
    uu = uu[indexBound]
    vv = vv[indexBound]
    
    #return results
    return(uu, vv, indexBound)
