function nurbs = BspGcMak(coefs,knots) 
% 
% Function Name: 
%  
%   BspGcMak - Construct the B-Spline for generalized coordinate space
%  
% Calling Sequence: 
%  
%   nurbs   = nrbmak(cntrl,knots); 
%  
% Parameters: 
%  
%   cntrl       : Control points, these can be either Cartesian or 
% 		homogeneous coordinates. 
%  
% 		For a curve the control points are represented by a 
% 		matrix of size (dim,nu) and for a surface a multidimensional 
% 		array of size (dim,nu,nv). Where nu is number of points along 
% 		the parametric U direction, and nv the number of points 
%               along the V direction. Dim is the dimension valid options 
% 		are 
% 		2 .... (x,y)        2D Cartesian coordinates 
% 		3 .... (x,y,z)      3D Cartesian coordinates 
% 		4 .... (wx,wy,wz,w) 4D homogeneous coordinates 
%  
%   knots	: Non-decreasing knot sequence spanning the interval 
%               [0.0,1.0]. It's assumed that the curves and surfaces 
%               are clamped to the start and end control points by knot 
%               multiplicities equal to the spline order. 
%               For curve knots form a vector and for a surface the knot 
%               are stored by two vectors for U and V in a cell structure. 
%               {uknots vknots} 
%                
%   nurbs 	: Data structure for representing a NURBS curve. 
%  
% NURBS Structure: 
%  
%   Both curves and surfaces are represented by a structure that is 
%   compatible with the Spline Toolbox from Mathworks 
%  
% 	nurbs.form   .... Type name 'B-NURBS' 
% 	nurbs.dim    .... Dimension of the control points 
% 	nurbs.number .... Number of Control points 
%       nurbs.coefs  .... Control Points 
%       nurbs.order  .... Order of the spline 
%       nurbs.knots  .... Knot sequence 
%  
%   Note: the control points are always converted and stored within the 
 
%  D.M. Spink 
%  Copyright (c) 2000. 
 
nurbs.form   = 'BspGc'; 
np = size(coefs); 
dim = np(1);
% nurbs.dim = int16(dim);  
nurbs.dim = dim; 
 
% constructing a curve 
nurbs.coefs = coefs; 
% nurbs.order = int16(size(knots,2)-np(2)); 
nurbs.order = size(knots,2)-np(2); 
nurbs.knots = knots;
% knots = sort(knots); 
% nurbs.knots = (knots-knots(1))/(knots(end)-knots(1)); 
end

 
 
 
