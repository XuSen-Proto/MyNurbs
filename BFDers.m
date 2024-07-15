function Ders = BFDers(i,u,p,U,ord)                 
% BASISFUN  Basis function for B-Spline 
% ------------------------------------------------------------------------- 
% ADAPTATION of BASISFUN from C Routine 
% ------------------------------------------------------------------------- 
% 
% Calling Sequence: 
%  
%   N = basisfun(i,u,p,U,ord) 
%    
%    INPUT: 
%    
%      i - knot span  ( from FindSpan() ) 
%      u - parametric point 
%      p - spline degree 
%      U - knot sequence 
%      ord - derivative order
%    
%    OUTPUT: 
%    
%      N - Basis functions and its 0th to ordth derivative (p+1)*ord matrix
%    
%    Algorithm A2.2 from 'The NURBS BOOK' pg70. 
                                                                                                                                                   
                                                   
i = i + 1; 
                                                  
left = zeros(p+1,1,'like',u);                              
right = zeros(p+1,1,'like',u);                              

N = zeros(p+1,'like',u);
N(1,1) = 1;                                        
for j=1:p                                          
    jc = j+1;
    left(jc) = u - U(i+1-j);                     
    right(jc) = U(i+j) - u;                      
    saved = 0;                                    
 
    for r=0:j-1                                   
        rc = r+1;
        temp = right(rc+1) + left(jc-r);
        N(jc,rc) = temp;
        temp = N(rc,jc-1)/temp;
        N(rc,jc) = saved + right(rc+1)*temp;        
        saved = left(jc-r)*temp;                  
    end                                           
 
    N(jc,jc) = saved;                            
end                                             
   
Ders = zeros(ord+1,p+1,'like',u);
Ders(1,:) = N(:,end)';
if ord>0
    as = zeros(2,p+1,'like',u);
    for r = 0:p
        s1 = int16(1);
        s2 = int16(2);
        rc = r+1;
        as(1,1) = 1;
        for k = 1:ord
            d = 0;
            kc = k+1;
            rk = rc-k;
            pk = p-k+1;
            if r >= k
                temp = as(s1,1)/N(pk+1,rk);
                as(s2,1) = temp;
                d = temp*N(rk,pk);
            end
            if (rk >= 0)
                j1 = int16(2);
            else
                j1 = 2-rk;
            end
            if (r < pk)
                j2 = k;
            else
                j2 = p-r+1;
            end
            for j = j1:j2
                jl = j-1;
                as(s2,j) = (as(s1,j)-as(s1,jl))/N(pk+1,rk+jl);
                d = d + as(s2,j)*N(rk+jl,pk);
            end
            if (r <= pk-1)
                as(s2,kc) = -as(s1,kc-1)/N(pk+1,rc);
                d = d + as(s2,kc)*N(rc,pk);
            end
            Ders(kc,rc) = d;
            j = s1;
            s1 = s2;
            s2 = j;
        end
    end
    r = p;
    for k = 1:ord
        Ders(k+1,:) = double(r)*Ders(k+1,:);
        r = r*(p-k);
    end
end
end


