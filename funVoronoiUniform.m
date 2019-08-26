% This code places uniformly points on cells of a Voronoi tesselation. 
% This code generates a Voronoi tessellation based on an artibrary point
% pattern. The Voronoi tesselation is found using the MATLAB function[1], 
% which is based on the Qhull project[2]. 
%
% INPUTS:
% xx and yy are vectors correspondong to the Cartesian 
% coordinates of the points in a point pattern. xx and yy must have the
% same length. 
% 
% vertexAll and cellAll are data structures that describe the Voronoi
% tesselation; see MATLAB function voronoin[1].
% 
% OUTPUTS: 
% uu and vv are vectors corresponding to the Cartesian coordinates of the
% uniformly placed points. 
%
% indexBounded is an index array for the bounded cells
%
% EXAMPLE: consider a point pattern described by 1-D arrays xx and yy,
% correspondong to the Cartesian coordinates. Then run the MATLAB function 
% 
% [vertexAll,cellAll]=voronoin([xx(:) yy(:)]);
%
% Then the uniform points on bounded cells are obtained with:
% [uu,vv]=funVoronoiUniform(xx,yy,vertexAll,cellAll);
% 
%
% All points (or cells) of the point process are numbered arbitrarily.
% In each *bounded* Voronoi cell a new point is uniformly placed.
%
% The placement step is done by first dividing each (bounded) Voronoi cell
% (ie an irregular polygon) with, say, m sides into m scalene triangles.
% The i th triangle is then randomly chosen based on the ratio of areas.
% A point is then uniformly placed on the i th triangle (via eq. 1 in [3]).
% The placement step is repeated for all bounded Voronoi cells.
%
% Author: H.Paul Keeler, 2019
%
% [1] http://www.mathworks.com.au/help/matlab/ref/voronoin.html
% [2] http://www.qhull.org/
% [2] Osada, R., Funkhouser, T., Chazelle, B. and Dobkin. D.,
% "Shape distributions", ACM Transactions on Graphics, vol 21, issue 4,
% 2002

function [uu,vv, indexBounded]=funVoronoiUniform(vertexAll,cellAll,xx,yy)
%reshape arrays (possibly not necessary if they are already column arrays)
xx=xx(:); yy=yy(:);

numbCells=length(xx); %number of Voronoi cells (including unbounded)

%initiate arrays
booleBounded=zeros(numbCells,1);
uu=zeros(numbCells,1);
vv=zeros(numbCells,1);

%Loop through for all Voronoi cells
for ii=1:numbCells
    %check for unbounded cells (index 1 corresponds to a point at infinity)
    booleBounded(ii)=~(any((cellAll{ii})==1));
    
    %checks if the Voronoi cell is bounded. if bounded, calculates its area
    %and assigns a single point uniformally in the Voronoi cell.
    
    %%% START -- Randomly place a point in a Voronoi cell -- START%%%
    if booleBounded(ii)
        xx0=xx(ii);yy0=yy(ii); %the (Poisson) point of the Voronoi cell
        %x/y values for current cell
        xxCell=vertexAll(cellAll{ii},1);
        yyCell=vertexAll(cellAll{ii},2);
        
        %%% START -- Caclulate areas of triangles -- START%%%
        %calculate area of triangles of bounded cell (ie irregular polygon)
        numbTri=length(xxCell); %number of triangles
        %create some indices for calculating triangle areas
        indexVertex=1:(numbTri+1); %number vertices
        indexVertex(end)=1; %repeat first index (ie returns to the start)
        indexVertex1=indexVertex(1:numbTri); %first vertex index
        indexVertex2=indexVertex(2:numbTri+1);  %second vertex index
        %using area equation for a triangle
        areaTri=abs((xxCell(indexVertex1)-xx0).*(yyCell(indexVertex2)-yy0)...
            -(xxCell(indexVertex2)-xx0).*(yyCell(indexVertex1)-yy0))/2;
        areaPoly=sum(areaTri);
        %%%END-- Caclulate areas of triangles -- END%%%
        
        %%%START -- Randomly placing point -- START%%%
        %%% place a point uniformaly in the (bounded) polygon
        %randomly choose the triangle (from the set that forms the polygon)
        cdfArea=cumsum(areaTri)/areaPoly; %create triangle CDF
        indexTri=find(rand(1,1)<=cdfArea,1); %use CDF to choose #
        
        indexVertex1=indexVertex(indexTri); %first vertex index
        indexVertex2=indexVertex(indexTri+1); %second vertex index
        %for all triangles
        xxTri=[xx0, xxCell(indexVertex1),xxCell(indexVertex2)];
        yyTri=[yy0, yyCell(indexVertex1),yyCell(indexVertex2)];
        
        %create two uniform random variables on unit interval
        uniRand1=rand(1,1); uniRand2=rand(1,1);
        
        %point is uniformally placed in the triangle via equation (1)in [2]
        %x coordinate
        uu(ii)=(1-sqrt(uniRand1))*xxTri(1)...
            +sqrt(uniRand1)*(1-uniRand2)*xxTri(2)...
            +sqrt(uniRand1)*uniRand2*xxTri(3);
        %y coordinate
        vv(ii)=(1-sqrt(uniRand1))*yyTri(1)...
            +sqrt(uniRand1)*(1-uniRand2)*yyTri(2)...
            +sqrt(uniRand1)*uniRand2*yyTri(3);
        %%%END -- Randomly placing point -- END%%%
        
    end
    %%% END -- Randomly place a point in a Voronoi cell -- END%%%
end

indexBounded=find(booleBounded==1); %find bounded cells
uu=uu(indexBounded); vv=vv(indexBounded); %remove unbounded cells
end