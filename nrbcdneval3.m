function vals = nrbcdneval3(nurbs, tt, selorder) 
% Evaluation of the derivative NURBS curve. 
% 
%     [pnt, jac] = nrbdeval(crv, dcrv, tt) 
% 
% INPUTS: 
% 
%   crv    - original NURBS curve. 
% 
%   srf    - original NUBRS surface 
% 
%   dcrv   - NURBS derivative represention of crv 
% 
%   dsrf   - NURBS derivative represention of surface 
% 
%   tt     - parametric evaluation points 
%            If the nurbs is a surface then tt is a cell 
%            {tu, tv} are the parametric coordinates 
%   3      - Derivative order
% 
% 
% Examples: 
%  
%   // Determine the first derivatives a NURBS curve at 9 points for 0.0 to 
%   // 1.0 
%   tt = linspace(0.0, 1.0, 9); 
%   dcrv = nrbderiv(crv); 
%   [pnts,jac] = nrbdeval(crv, dcrv, tt); 
 
%  D.M. Spink 
%  Copyright (c) 2000. 

lent = length(tt);
vals = zeros(3,lent,4);
ws = zeros(4,lent);
[vals(:,:,1),ws(1,:)] = nrbeval(nurbs{1}, tt);  
% NURBS is a curve 

%   temp = cw(ones(3,1),:); 
%   pnt = cp./temp; 
vals(:,:,1) = vals(:,:,1)./(ones(3,1)*ws(1,:)); 

pasc = [ones(2,1); zeros(3,1)];
for i = 2:4
  [a,ws(i,:)] = nrbeval(nurbs{i},tt); 
  w0 = bsxfun(@times,pasc(2:i),ws(2:i,:))';
  wtmp = reshape(kron(w0,ones(3,1)),3,[],i-1);
  v0 = bsxfun(@times,wtmp,vals(:,:,i-1:-1:1));
  vals(:,:,i) = bsxfun(@rdivide,a-sum(v0,3),ws(1,:)); 
  pasc1 = [pasc(1:i);0]; pasc2 = [0;pasc(1:i)];
  pasc(1:i+1) = pasc1+pasc2;
end
if ~isempty(selorder)
  vals = vals(:,:,selorder);
end
 
