% This code tests the function funVoronoiUniform, which places a single
% random point uniformly on *boundeded* Voronoi cells in two dimensions.
%
% For the test, a single realization of a (homogeneous) Poisson point
% process (PPP) is first generated on a rectangle (but any (finite)
% two-dimensional point pattern can be used). Then the Voronoi tesselation
% is created by using the MATLAB function[1], which is based on the Qhull
% project[2]. Then the function funVoronoiUniform is (repeatedly) used to
% uniformly place a single random point in each *bounded* Voronoi cell.
% The averages of these random single points are then taken. As the number
% of simulations (of placing single points) increases, these averages
% should converge to the centroids (or geometric centres) of all the
% *bounded* Voronoi cells.
%
% The placement step is done by first dividing each *bounded* Voronoi cell
% (ie an irregular polygon) with, say, m sides into m scalene triangles.
% The i-th triangle is then randomly chosen based on the ratio of areas.
% A point is then uniformly placed on the i th triangle (via eq. 1 in [3]).
% The placement step is repeated for all bounded Voronoi cells.
%
% All the cells and points of the underlying point pattern are numbered or
% indexed arbitrarily. The numbering/indexing of cells and points agree.
%
% A Voronoi diagram is displayed over the underlying point pattern.
% Points of the underlying point pattern are marked green and blue if they
% are located respectively in bounded and unbounded Voronoi cells.
% The uniformly placed points in the bounded cells are marked red.
%
% If there are no bounded Voronoi cells, no diagram is created.
%
% Author: H.Paul Keeler, 2019.
% hpaulkeeler.com
% github.com/hpaulkeeler/voronoi_uniform
%
% References:
% [1] http://www.mathworks.com.au/help/matlab/ref/voronoin.html
% [2] http://www.qhull.org/
% [2] Osada, R., Funkhouser, T., Chazelle, B. and Dobkin. D.,
% "Shape distributions", ACM Transactions on Graphics, vol 21, issue 4,
% 2002

clearvars; close all; clc;

boolePlot=1; % set to 1 for plot, 0 for no plot
numbSim=10^3; %number of random point placements (given one point pattern)

%Simulation window parameters
xMin=0;
xMax=1;
yMin=0;
yMax=1;
%rectangle dimensions
xDelta=xMax-xMin; %width
yDelta=yMax-yMin; %height
areaTotal=xDelta*yDelta; %area of rectangle

%Point process parameters
lambda=12; %intensity (ie mean density) of the Poisson point process

%Simulate Poisson point process
numbPoints=poissrnd(areaTotal*lambda);%Poisson number of points
xx=xDelta*(rand(numbPoints,1))+xMin;%x coordinates of Poisson points
yy=xDelta*(rand(numbPoints,1))+yMin;%y coordinates of Poisson points

%%TEMP: A non-random example of a point pattern can be used.
%xx=[1,3,5,6,3,2,12,1,3,8]/15;
%yy=[2,3,5,1,7,4,3,14,14,13]/15;
%numbPoints=length(xx);

numbCells=numbPoints; %number of Voronoi cells (including unbounded)

%reshape vectors (possibly not necessary if they are already column vectors)
xx=xx(:); yy=yy(:);
xxyy=[xx yy]; %combine x and y coordinates
%Perform Voronoi tesseslation using built-in function voronoin
[vertexAll,cellAll]=voronoin(xxyy);

%initialize  arrays for empirical estimates of centroids
xCentEmp=zeros(numbCells,1); %x component
yCentEmp=zeros(numbCells,1); %y component
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

%initialize  arrays for (analtic) calculations of centroids
xCentExact=zeros(numbBound,1); %x component
yCentExact=zeros(numbBound,1); %y component
%loop through for all bounded cells and calculate centroids
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
    %create voronoi diagram on the point pattern
    voronoi(xx,yy);
    %plot the underlying point pattern
    scatter(xx,yy,'bo');
    %put a red o uniformly in each bounded Voronoi cell
    scatter(uu,vv,'ro');
    % put a green * on the base station of each Voronoi bounded cell
    scatter(xx(indexBound),yy(indexBound),'g*');
    
    %number the points/cells
    labels=cellstr(num2str((1:numbPoints)'));%labels correspond to their order
    text(xx, yy, labels, 'VerticalAlignment','bottom', ...
        'HorizontalAlignment','right');
end
%%%END -- Plotting section -- END%%%