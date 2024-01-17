function p = bspevals(d,c,k,us) 
%  
% Function Name: 
%  
%   bspeval - Evaluate a univariate B-Spline. 
%  
% Calling Sequence: 
%  
%   p = bspeval(d,c,k,u) 
%  
% Parameters: 
%  
%   d	: Degree of the B-Spline. 
%  
%   c	: Control Points, matrix of size (dim,nc). 
%  
%   k	: Knot sequence, row vector of size nk. 
%  
%   us	: Parametric evaluation points, u1<u2<...un, row vector of size nu. 
%  
%   p	: Evaluated points, matrix of size (dim,nu) 
%  
% Description: 
%  
%   Evaluate a univariate B-Spline. This function provides an interface to 
%   a toolbox 'C' routine. 
nu = numel(us); 
% [mc,nc] = size(c); 
mc = size(c,1); nc = size(c,2);
                                                %   int bspeval(int d, double *c, int mc, int nc, double *k, int nk, double *u,int nu, double *p){ 
                                                %   int ierr = 0; 
                                                %   int i, s, tmp1, row, col; 
                                                %   double tmp2; 
                                                % 
                                                %   // Construct the control points 
                                                %   double **ctrl = vec2mat(c,mc,nc); 
                                                % 
                                                %   // Contruct the evaluated points 
p = zeros(mc,nu);                               %   double **pnt = vec2mat(p,mc,nu); 
                                                % 
                                                %   // space for the basis functions 
N = zeros(1,d+1);                               %   double *N = (double*) mxMalloc((d+1)*sizeof(double)); 
s = findspan(nc-1, d, us(1), k);                % 
int = s(1);                                     %   // for each parametric point i 
for col=1:nu                                    %   for (col = 0; col < nu; col++) { 
                                                %     // find the span of u[col] 
%     s = findspan(nc-1, d, us(col), k);           %     s = findspan(nc-1, d, u[col], k); 
    if us(col)>=k(int+2)
        int = int+1;
    end
    N(1:d+1) = basisfun(s(1),us(col),d,k);                 %     basisfun(s, u[col], d, k, N); 
                                                % 
    tmp1 = s(1) - d + 1;                           %     tmp1 = s - d; 
%     c
%     for row=1:mc                                %     for (row = 0; row < mc; row++)  { 
%         tmp2 = 0;                               %       tmp2 = 0.0; 
%         for i=0:d                               %       for (i = 0; i <= d; i++) 
%            tmp2 = tmp2 + N(i+1)*c(row,tmp1+i);  % 	tmp2 += N[i] * ctrl[tmp1+i][row]; 
%         end                                     % 
%         p(row,col) = tmp2;                      %       pnt[col][row] = tmp2; 
%     end                                         %     } 
    p(:,col) = c(:,tmp1:tmp1+d)*(N(1:d+1)');
%     p(:,col) = sum(bsxfun(@times,c(1:4,tmp1:tmp1+d),N(1:d+1)),2);
end                                             %   } 
                                                % 
                                                %   mxFree(N); 
                                                %   freevec2mat(pnt); 
                                                %   freevec2mat(ctrl); 
                                                % 
                                                %   return ierr; 
                                                %   } 
