# This code generates a Voronoi-Poisson tessellation and places a point
# uniformly on each cell. 
# A (homogeneous) Poisson point process (PPP) is first created on a 
# rectangle.Then the Voronoi tesselation is found using the MATLAB 
# function[1], which is based on the Qhull project[2] .
# All points (or cells) of the PPP are numbered arbitrarily.
# In each *bounded* Voronoi cell a new point is uniformly placed.
#
# The placement step is done by first dividing each (bounded) Voronoi cell
# (ie an irregular polygon) with, say, m sides into m scalene triangles.
# The i th triangle is then randomly chosen based on the ratio of areas.
# A point is then uniformly placed on the i th triangle (via eq. 1 in [3]).
# The placement step is repeated for all bounded Voronoi cells.
# A Voronoi diagram is displayed over the PPP. Points of the underlying PPP
# are marked green and blue if they are located respectively in bounded and
# unbounded Voronoi cells. The uniformally placed points in the bounded
# cells are marked red.
#
# If there are no bounded Voronoi cells, no diagram is created.
# Author: H.Paul Keeler, 2019
#
# [1] http://scipy.github.io/devdocs/generated/scipy.spatial.Voronoi.html
# [2] http://www.qhull.org/
# [2] Osada, R., Funkhouser, T., Chazelle, B. and Dobkin. D.,
# "Shape distributions", ACM Transactions on Graphics, vol 21, issue 4,
# 2002

import numpy as np;  # NumPy package for arrays, random number generation, etc
import matplotlib.pyplot as plt  # for plotting
from scipy.spatial import Voronoi, voronoi_plot_2d #for voronoi tessellation
from funVoronoiUniform import funVoronoiUniform


plt.close('all');  # close all figures

boolePlot=True; # set to True for plot, False for no plot
numbSim=10**4; #number of random point placements (given one point pattern)

# Simulation window parameters
xMin = 0;
xMax = 1;
yMin = 0;
yMax = 1;
xDelta = xMax - xMin;
yDelta = yMax - yMin;  # rectangle dimensions
areaTotal = xDelta * yDelta;

# Point process parameters
lambda0 = 12;  # intensity (ie mean density) of the Poisson process

# Simulate a Poisson point process
numbPoints = np.random.poisson(lambda0 * areaTotal);  # Poisson number of points
xx = xDelta * np.random.uniform(0, 1, numbPoints) + xMin;  # x coordinates of Poisson points
yy = yDelta * np.random.uniform(0, 1, numbPoints) + yMin;  # y coordinates of Poisson points

#TEMP: a non-random example
#xx=np.array([1,3,5,6,3,2,12,1,3,8]);
#yy=np.array([2,3,5,1,7,4,3,14,14,13]);
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

#initiate arrays
xCentEmp=np.zeros(numbCells);
yCentEmp=np.zeros(numbCells);
#Loop through for multiple simulations 
for ss in range(numbSim):
    #randomly place a point on all the bounded Voronoi cells
    [uu, vv,indexBounded]=funVoronoiUniform(voronoiData,xx,yy);
    xCentEmp[indexBounded]=xCentEmp[indexBounded]+uu;
    yCentEmp[indexBounded]=yCentEmp[indexBounded]+vv;

### END -- Randomly place a point in a Voronoi cell -- END###
numbBound=len(indexBounded); #number of bounded cells

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
    
#calculate centroids
#initiate arrays    
xCentExact=np.zeros(numbBound); yCentExact=np.zeros(numbBound);
for ii in range(numbBound):
    xxCell=vertexAll[cellAll[indexP2C[indexBounded[ii]]],0];
    yyCell=vertexAll[cellAll[indexP2C[indexBounded[ii]]],1];
    xCent,yCent = funCentroid(xxCell,yyCell);
    xCentExact[ii]=xCent; yCentExact[ii]=yCent;

####START -- Plotting section -- START###
if (numbBound>0) and (boolePlot):
#    figure; grid; hold on;
    #plotting points in PPP
    voronoi_plot_2d(voronoiData); #create voronoi diagram on the PPP
    #put a red o uniformally in each bounded Voronoi cell
    plt.scatter(uu, vv, edgecolor='r', facecolor='none');
    #put a green star on the base station of each Voronoi bounded cell
    plt.scatter(xx[indexBounded], yy[indexBounded], edgecolor='g', marker='*');
    #number the points
    for ii in range(numbPoints):
        plt.text(xx[ii], yy[ii], ii);   
   
####END -- Plotting section -- END###