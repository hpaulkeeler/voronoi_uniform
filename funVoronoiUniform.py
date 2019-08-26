# This code places uniformly points on cells of a Voronoi tesselation. 
# This code generates a Voronoi tessellation based on an artibrary point
# pattern. The Voronoi tesselation is found using the SciPy function[1], 
# which is based on the Qhull project[2]. 
#
# INPUTS:
# xx and yy are vectors correspondong to the Cartesian 
# coordinates of the points in a point pattern. xx and yy must have the
# same length. 
# 
# voronoiData is data structure that describes the Voronoi
# tesselation; see SciPy function voronoin[1].
# 
# OUTPUTS: 
# uu and vv are vectors corresponding to the Cartesian coordinates of the
# uniformly placed points. 
#
# indexBounded is an index array for the bounded cells
#
# EXAMPLE: consider a point pattern described by 1-D arrays xx and yy,
# correspondong to the Cartesian coordinates. Then run the SciPy function 
# 
# voronoiData=Voronoi(np.stack((xx,yy), axis=1));
#
# Then the uniform points on bounded cells are obtained with:
# uu,vv, indexBounded=funVoronoiUniform(xx,yy,voronoiData);
# 
# uu and vv are vectors corresponding to the Cartesian coordinates of the
# uniformly placed points. 
#
# indexBounded is an index array for the bounded cells
#
# All points (or cells) of the point process are numbered arbitrarily.
# In each *bounded* Voronoi cell a new point is uniformly placed.
#
# The placement step is done by first dividing each (bounded) Voronoi cell
# (ie an irregular polygon) with, say, m sides into m scalene triangles.
# The i th triangle is then randomly chosen based on the ratio of areas.
# A point is then uniformly placed on the i th triangle (via eq. 1 in [3]).
# The placement step is repeated for all bounded Voronoi cells.
#
# Author: H.Paul Keeler, 2019
#
# [1] http://scipy.github.io/devdocs/generated/scipy.spatial.Voronoi.html
# [2] http://www.qhull.org/
# [2] Osada, R., Funkhouser, T., Chazelle, B. and Dobkin. D.,
# "Shape distributions", ACM Transactions on Graphics, vol 21, issue 4,
# 2002

import numpy as np;  # NumPy package for arrays, random number generation, etc
from scipy.spatial import Voronoi, voronoi_plot_2d #for voronoi tessellation

def funVoronoiUniform(voronoiData,xx,yy):
    #reshape arrays (possibly not necessary if they are already row arrays)
    xx=xx.flatten(); yy=yy.flatten(); 
    
    numbCells=len(xx); #number of Voronoi cells (including unbounded)
    
    vertexAll=voronoiData.vertices; #retrieve x/y coordiantes of all vertices
    cellAll=voronoiData.regions; #may contain empty array/set
    
    indexP2C=voronoiData.point_region; #index mapping between cells and points
    
    ##initiate arrays
    booleBounded=np.zeros(numbCells,dtype=bool);
    uu=np.zeros(numbCells);
    vv=np.zeros(numbCells);
    
    ##Loop through for all Voronoi cells
    for ii in range(numbCells):        
        booleBounded[ii]=not(any(np.array(cellAll[indexP2C[ii]])==-1));
        #checks if the Voronoi cell is bounded. if bounded, calculates its area 
        #and assigns a single point uniformally in the Voronoi cell.
        
        ### START -- Randomly place a point in a Voronoi cell -- START###
        if booleBounded[ii]:               
            xx0=xx[ii];yy0=yy[ii]; #the (Poisson) point of the Voronoi cell
            
            #print(jj,xx0,yy0)        
            #x/y values for current cell
            xxCell=vertexAll[cellAll[indexP2C[ii]],0];        
            yyCell=vertexAll[cellAll[indexP2C[ii]],1];
                    
            ### START -- Caclulate areas of triangles -- START###
            #calculate area of triangles of bounded cell (ie irregular polygon)
            numbTri=len(xxCell); #number of triangles        
            
            #create some indices for calculating triangle areas
            indexVertex= np.arange(numbTri+1); #number vertices
            indexVertex[-1]=0; #repeat first index (ie returns to the start)
            indexVertex1=indexVertex[np.arange(numbTri)]; #first vertex index
            indexVertex2=indexVertex[np.arange(numbTri)+1];  #second vertex index
            #using area equation for a triangle
            areaTri=abs((xxCell[indexVertex1]-xx0)*(yyCell[indexVertex2]-yy0)\
                        -(xxCell[indexVertex2]-xx0)*(yyCell[indexVertex1]-yy0))/2;            
            areaPoly=sum(areaTri);        
            ###END-- Caclulate areas of triangles -- END###
            
            ###START -- Randomly placing point -- START###
            ### place a point uniformaly in the (bounded) polygon
            #randomly choose the triangle (from the set that forms the polygon)
            cdfArea=np.cumsum(areaTri)/areaPoly; #create triangle CDF
            indexTri=(np.random.rand() <= cdfArea).argmax(); #use CDF to choose #
                        
            indexVertex1=indexVertex[indexTri]; #first vertex index
            indexVertex2=indexVertex[indexTri+1]; #second vertex index
            #for all triangles except the last one
            xxTri=[xx0, xxCell[indexVertex1],xxCell[indexVertex2]];
            yyTri=[yy0, yyCell[indexVertex1],yyCell[indexVertex2]];
            
            #create two uniform random variables on unit interval
            uniRand1=np.random.rand(); uniRand2=np.random.rand();
            
            #x coordinate
            uu[ii]=(1-np.sqrt(uniRand1))*xxTri[0]\
            +np.sqrt(uniRand1)*(1-uniRand2)*xxTri[1]\
            +np.sqrt(uniRand1)*uniRand2*xxTri[2]
            #y coordinate
            vv[ii]=(1-np.sqrt(uniRand1))*yyTri[0]\
            +np.sqrt(uniRand1)*(1-uniRand2)*yyTri[1]\
            +np.sqrt(uniRand1)*uniRand2*yyTri[2];
            ###END -- Randomly placing point -- END###        

    ### END -- Randomly place a point in a Voronoi cell -- END###
    
    indexBounded=np.arange(numbCells)[booleBounded]; #find bounded cells
    uu=uu[indexBounded]; vv=vv[indexBounded]; #remove unbounded cells
    return(uu,vv,indexBounded)