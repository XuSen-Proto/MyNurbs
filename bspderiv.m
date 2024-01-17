function   [dc,dk] = bspderiv(d,c,k) 
%  
% Function Name: 
%  
%   bspdeval - Evaluate the control points and knot sequence of the derivative 
%              of a univariate B-Spline. 
%  
% Calling Sequence: 
%  
%   [dc,dk] = bspderiv(d,c,k) 
%  
% Parameters: 
%  
%   d	: Degree of the B-Spline. 
%  
%   c	: Control Points, matrix of size (dim,nc). 
%  
%   k	: Knot sequence, row vector of size nk. 
%  
%   dc	: Control points of the derivative 
%  
%   dk	: Knot sequence of the derivative 
%  
% Description: 
%  
%   Evaluate the derivative of univariate B-Spline, which is itself a B-Spline. 
%   This function provides an interface to a toolbox 'C' routine. 
[mc,nc] = size(c); 
% nc = size(c,2);
nk = numel(k); 
                                                     % 
                                                     % int bspderiv(int d, double *c, int mc, int nc, double *k, int nk, double *dc, 
                                                     %              double *dk) 
                                                     % { 
                                                     %   int ierr = 0; 
                                                     %   int i, j, tmp; 
                                                     % 
                                                     %   // control points 
                                                     %   double **ctrl = vec2mat(c,mc,nc); 
                                                     % 
                                                     %   // control points of the derivative 
% dc = zeros(mc,nc-1);                                 %   double **dctrl = vec2mat(dc,mc,nc-1); 
                                                     % 
% for i=0:nc-2                                         %   for (i = 0; i < nc-1; i++) { 
%    tmp = d / (k(i+d+2) - k(i+2));                    %     tmp = d / (k[i+d+1] - k[i+1]); 
   % for j=0:mc-1                                      %     for (j = 0; j < mc; j++) { 
   %     dc(j+1,i+1) = tmp*(c(j+1,i+2) - c(j+1,i+1));  %       dctrl[i][j] = tmp * (ctrl[i+1][j] - ctrl[i][j]); 
   % end                                               %     } 
%    dc(：,i+1) = tmp*(c(：,i+2) - c(：,i+1));
% end                                                  %   } 
%                                                      % 
% dk = zeros(1,nk-2);                                  %   j = 0; 
% for i=1:nk-2                                         %   for (i = 1; i < nk-1; i++) 
%    dk(i) = k(i+1);                                   %     dk[j++] = k[i]; 
% end                                                  % 
%                                                      %   freevec2mat(dctrl); 
%                                                      %   freevec2mat(ctrl); 
%                                                      % 
%                                                      %   return ierr; 
%                                                      % }</pre>
% for i=2:nc
%   tmp = d / (k(i+d) - k(i));
%   dc(:,i-1) = tmp*(c(:,i) - c(:,i-1));
% end
dc = bsxfun(@times,(bsxfun(@rdivide,d,(k((2:nc)+d) - k(2:nc)))),...
    (c(:,2:nc) - c(:,1:nc-1)));
dk = k(2:nk-1);