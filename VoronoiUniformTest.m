% This code generates a Voronoi-Poisson tessellation and places a point
% uniformly on each cell. 
% A (homogeneous) Poisson point process (PPP) is first created on a 
% rectangle.Then the Voronoi tesselation is found using the MATLAB 
% function[1], which is based on the Qhull project[2] .
% All points (or cells) of the PPP are numbered arbitrarily.
% In each *bounded* Voronoi cell a new point is uniformly placed.
%
% The placement step is done by first dividing each (bounded) Voronoi cell
% (ie an irregular polygon) with, say, m sides into m scalene triangles.
% The i th triangle is then randomly chosen based on the ratio of areas.
% A point is then uniformly placed on the i th triangle (via eq. 1 in [3]).
% The placement step is repeated for all bounded Voronoi cells.
% A Voronoi diagram is displayed over the PPP. Points of the underlying PPP
% are marked green and blue if they are located respectively in bounded and
% unbounded Voronoi cells. The uniformally placed points in the bounded
% cells are marked red.
%
% If there are no bounded Voronoi cells, no diagram is created.
%
% Author: H.Paul Keeler, 2019
%
% [1] http://www.mathworks.com.au/help/matlab/ref/voronoin.html
% [2] http://www.qhull.org/
% [2] Osada, R., Funkhouser, T., Chazelle, B. and Dobkin. D.,
% "Shape distributions", ACM Transactions on Graphics, vol 21, issue 4,
% 2002

clearvars; close all; clc;

boolePlot=1; % set to 1 for plot, 0 for no plot
numbSim=10^3; %number of random point placements (given one point pattern)

%Simulation window parameters
xMin=0;xMax=1;
yMin=0;yMax=1;
xDelta=xMax-xMin;yDelta=yMax-yMin; %rectangle dimensions
areaTotal=xDelta*yDelta; %area of rectangle

%Point process parameters
lambda=12; %intensity (ie mean density) of the Poisson process

%Simulate Poisson point process
numbPoints=poissrnd(areaTotal*lambda);%Poisson number of points
xx=xDelta*(rand(numbPoints,1))+xMin;%x coordinates of Poisson points
yy=xDelta*(rand(numbPoints,1))+yMin;%y coordinates of Poisson points

%TEMP: a non-random example
%xx=[1,3,5,6,3,2,12,1,3,8];
%yy=[2,3,5,1,7,4,3,14,14,13];
%numbPoints=length(xx);

numbCells=numbPoints; %number of Voronoi cells (including unbounded)

%reshape vectors (possibly not necessary if they are already column vectors)
xx=xx(:); yy=yy(:);
xxyy=[xx yy]; %combine x and y coordinates
%Perform Voronoi tesseslation using built-in function
[vertexAll,cellAll]=voronoin(xxyy);

%initiate arrays
xCentEmp=zeros(numbCells,1);
yCentEmp=zeros(numbCells,1);
%Loop through for multiple simulations
for ss=1:numbSim
    %randomly place a point on all the bounded Voronoi cells
    [uu, vv,indexBound]=funVoronoiUniform(vertexAll,cellAll,xx,yy);
    xCentEmp(indexBound)=xCentEmp(indexBound)+uu;
    yCentEmp(indexBound)=yCentEmp(indexBound)+vv;
end
numbBound=length(indexBound); %number of bounded cells

%estimate empirical centroids
xCentEmp=xCentEmp(indexBound)/numbSim;
yCentEmp=yCentEmp(indexBound)/numbSim;

%calculate centroids
xCentExact=zeros(numbBound,1); yCentExact=zeros(numbBound,1);
for ii=1:numbBound
    xxCell=vertexAll(cellAll{indexBound(ii)},1);
    yyCell=vertexAll(cellAll{indexBound(ii)},2);
    polyin = polyshape(xxCell,yyCell);
    [xCent,yCent] = centroid(polyin);
    xCentExact(ii)=xCent; yCentExact(ii)=yCent;
end

%%%START -- Plotting section -- START%%%
if (numbBound>0)&&boolePlot
    figure; grid; hold on;
    %plotting points in PPP
    scatter(xx,yy);
    voronoi(xx,yy); %create voronoi diagram on the PPP
    %put a red o uniformally in each bounded Voronoi cell
    scatter(uu,vv,'ro');
    % put a green * on the base station of each Voronoi bounded cell
    scatter(xx(indexBound),yy(indexBound),'g*');
    
    %number the points/cells
    labels=cellstr(num2str((1:numbPoints)'));%labels correspond to their order
    text(xx, yy, labels, 'VerticalAlignment','bottom', ...
        'HorizontalAlignment','right');
    
end
%%%END -- Plotting section -- END%%%