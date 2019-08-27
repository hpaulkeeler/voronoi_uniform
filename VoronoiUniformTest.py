# This code tests the function funVoronoiUniform, which places a single
# random point uniformly on *boundeded* Voronoi cells.
#
# For the test, a single realization of a (homogeneous) Poisson point
# process (PPP) is first generated on a rectangle (but any (finite) point
# pattern can be used). Then the Voronoi tesselation is created by using
# the SciPy function[1], which is based on the Qhull project[2]. Then the
# function funVoronoiUniform is (repeatedly) used to uniformly place a
# single random point in each *bounded* Voronoi cell. The averages of these
# random single points are then taken. As the number of simulations (of
# placing single points) increases, these averages should converge to the
# centroids (or geometric centres) of all the *bounded* Voronoi cells.
#
# The placement step is done by first dividing each *bounded* Voronoi cell
# (ie an irregular polygon) with, say, m sides into m scalene triangles.
# The i-th triangle is then randomly chosen based on the ratio of areas.
# A point is then uniformly placed on the i th triangle (via eq. 1 in [3]).
# The placement step is repeated for all bounded Voronoi cells.
#
# All the cells and points of the underlying point pattern are numbered or 
# indexed arbitrarily. The numbering/indexing of cells and points do NOT agree. 
#
# A Voronoi diagram is displayed over the underlying point pattern.
# Points of the underlying point pattern are marked green and blue if they
# are located respectively in bounded and unbounded Voronoi cells.
# The uniformly placed points in the bounded cells are marked red.
#
# If there are no bounded Voronoi cells, no diagram is created.
#
# Author: H.Paul Keeler, 2019.
# hpaulkeeler.com
# github.com/hpaulkeeler/voronoi_uniform
#
# References:
# [1] http://scipy.github.io/devdocs/generated/scipy.spatial.Voronoi.html
# [2] http://www.qhull.org/
# [2] Osada, R., Funkhouser, T., Chazelle, B. and Dobkin. D.,
# "Shape distributions", ACM Transactions on Graphics, vol 21, issue 4,
# 2002

import numpy as np;  # NumPy package for arrays, random number generation, etc
import matplotlib.pyplot as plt  # for plotting
from scipy.spatial import Voronoi, voronoi_plot_2d #for voronoi tessellation
from funVoronoiUniform import funVoronoiUniform #for placing uniform points


plt.close('all');  # close all figures

boolePlot=True; # set to True for plot, False for no plot
numbSim=10**3; #number of random point placements (given one point pattern)

# Simulation window parameters
xMin = 0;
xMax = 1;
yMin = 0;
yMax = 1;
 # rectangle dimensions
xDelta = xMax - xMin; #width
yDelta = yMax - yMin; #height
areaTotal = xDelta * yDelta;

# Point process parameters
lambda0 = 12;  # intensity (ie mean density) of the Poisson point process

# Simulate a Poisson point process
numbPoints = np.random.poisson(lambda0 * areaTotal);  # Poisson number of points
xx = xDelta * np.random.uniform(0, 1, numbPoints) + xMin;  # x coordinates of Poisson points
yy = yDelta * np.random.uniform(0, 1, numbPoints) + yMin;  # y coordinates of Poisson points

##TEMP: A non-random example of a point pattern can be used.
#xx=np.array([1,3,5,6,3,2,12,1,3,8])/15;
#yy=np.array([2,3,5,1,7,4,3,14,14,13])/15;
#numbPoints=len(xx);

numbCells=numbPoints; #number of Voronoi cells (including unbounded)

#reshape arrays (possibly not necessary if they are already row vectors)
xx=xx.flatten(); yy=yy.flatten(); 
xxyy=np.stack((xx,yy), axis=1); #combine x and y coordinates
##Perform Voroin tesseslation using built-in function
voronoiData=Voronoi(xxyy);
vertexAll=voronoiData.vertices; #retrieve x/y coordiantes of all vertices
cellAll=voronoiData.regions; #may contain empty array/set
numbCells=numbPoints; #number of Voronoi cells (including unbounded)
indexP2C=voronoiData.point_region; #index mapping between cells and points

#initiate arrays for empirical estimates of centroids
xCentEmp=np.zeros(numbCells); #x component 
yCentEmp=np.zeros(numbCells); #y component 
#Loop through for multiple simulations 
for ss in range(numbSim):
    #randomly place a point on all the bounded Voronoi cells
    [uu, vv,indexBounded]=funVoronoiUniform(voronoiData);
    xCentEmp[indexBounded]=xCentEmp[indexBounded]+uu;
    yCentEmp[indexBounded]=yCentEmp[indexBounded]+vv;

### END -- Randomly place a point in a Voronoi cell -- END###
numbBounded=len(indexBounded); #number of bounded cells

#estimate empirical centroids
xCentEmp=xCentEmp[indexBounded]/numbSim; 
yCentEmp=yCentEmp[indexBounded]/numbSim;

#create a function for calculating the centroid of a polgyon. 
#For the formula, see, for example:
#https://en.wikipedia.org/wiki/Centroid#Of_a_polygon
#http://demonstrations.wolfram.com/CenterOfMassOfAPolygon/
def funCentroid(x,y):          
    indexCent=np.arange(len(x)); #create index 
    #loop back arrays tos start
    x=np.append(x,x[0]); 
    y=np.append(y,y[0]);
    xyStep=x[indexCent]*y[indexCent+1]-x[indexCent+1]*y[indexCent];                
    areaPoly=sum(xyStep)/2
    xCentroid = sum((x[indexCent]+x[indexCent+1])*xyStep)/areaPoly/6;
    yCentroid = sum((y[indexCent]+y[indexCent+1])*xyStep)/areaPoly/6; 
    return(xCentroid, yCentroid)
    

#initiate arrays for (analtic) calculations of centroids
xCentExact=np.zeros(numbBounded); #x component 
yCentExact=np.zeros(numbBounded); #x component 
#loop through for all bounded cells and calculate centroids
for ii in range(numbBounded):
    xxCell=vertexAll[cellAll[indexP2C[indexBounded[ii]]],0];
    yyCell=vertexAll[cellAll[indexP2C[indexBounded[ii]]],1];
    xCent,yCent = funCentroid(xxCell,yyCell);
    xCentExact[ii]=xCent; yCentExact[ii]=yCent;

####START -- Plotting section -- START###
if (numbBounded>0) and (boolePlot):
    #create voronoi diagram on the point pattern
    voronoi_plot_2d(voronoiData, show_points=False,show_vertices=False); 
    #plot the underlying point pattern
    plt.scatter(xx, yy, edgecolor='b', facecolor='none');
    #put a red o uniformly in each bounded Voronoi cell
    plt.scatter(uu, vv, edgecolor='r', facecolor='none');
    #put a green star on the base station of each Voronoi bounded cell
    plt.scatter(xx[indexBounded], yy[indexBounded], color='g', marker='*');
    #number the points
    for ii in range(numbPoints):
        plt.text(xx[ii]+xDelta/50, yy[ii]+yDelta/50, ii);   
   
####END -- Plotting section -- END###