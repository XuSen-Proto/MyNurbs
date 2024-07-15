function [Ders,dDers] = BFDersDu(i,u,p,U,ord,du,dU)                 
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
nx = length(du);

left = zeros(p+1,1,'like',u);     
right = zeros(p+1,1,'like',u);
dleft = zeros(p+1,nx,'like',u);
dright = zeros(p+1,nx,'like',u);
dN = zeros(p+1,p+1,nx,'like',u);

N = zeros(p+1,'like',u);
N(1,1) = 1;                                        
for j=1:p                                          
    jc = j+1;
    left(jc) = u - U(i+1-j);                     
    right(jc) = U(i+j) - u;
    dleft(jc,:) = du - dU(i+1-j,:);                     
    dright(jc,:) = dU(i+j,:) - du;
    saved = 0;
    dsaved = zeros(1,1,nx);
    for r = 0:j-1                                   
        rc = r+1;
        temp = right(rc+1) + left(jc-r);
        dtemp = reshape(dright(rc+1,:) + dleft(jc-r,:), 1,1,[]);

        N(jc,rc) = temp;
        dN(jc,rc,:) = dtemp;
        
        temp2 = N(rc,jc-1)/temp;
        dtemp = (dN(rc,jc-1,:) - temp2.*dtemp)./temp;
        temp = temp2;

        N(rc,jc) = saved + right(rc+1)*temp;        
        saved = left(jc-r)*temp; 
        dN(rc,jc,:) = dsaved + reshape(dright(rc+1,:),1,1,[]).*temp + right(rc+1)*dtemp;
        dsaved = reshape(dleft(jc-r,:),1,1,[]).*temp + left(jc-r)*dtemp;
    end
    N(jc,jc) = saved;
    dN(jc,jc,:) = dsaved; 
end                                             
   
Ders = zeros(ord+1,p+1,1,'like',u);
Ders(1,:,:) = N(:,end)';
dDers = zeros(ord+1,p+1,nx,'like',u);
dDers(1,:,:) = permute(dN(:,end,:), [2 1 3]);
if ord>0
    as = zeros(2,p+1,'like',u);
    das = zeros(2,p+1,nx,'like',u);
    for r = 0:p
        s1 = int16(1);
        s2 = int16(2);
        rc = r+1;
        as(1,1) = 1;
        das(1,1,:) = 0;
        for k = 1:ord
            d = 0;
            dd = zeros(1,1,nx,'like',u);
            kc = k+1;
            rk = rc-k;
            pk = p-k+1;
            if r >= k
                Ntmp = 1/N(pk+1,rk);
                temp = as(s1,1)*Ntmp;
                dtemp = (das(s1,1,:) - temp.*dN(pk+1,rk,:)).*Ntmp;
                as(s2,1) = temp;
                das(s2,1,:) = dtemp;
                Ntmp = N(rk,pk);
                d = temp*Ntmp;
                dd = temp*dN(rk,pk,:) + dtemp*Ntmp;
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
                Ntmp = 1/N(pk+1,rk+jl);
                a_N = (as(s1,j)-as(s1,jl))*Ntmp;
                as(s2,j) = a_N;
                das(s2,j,:) = (das(s1,j,:)-das(s1,jl,:) - a_N.*dN(pk+1,rk+jl,:)).*Ntmp;
                Ntmp = N(rk+jl,pk);
                d = d + as(s2,j)*Ntmp;
                dd = dd + das(s2,j,:)*Ntmp+ a_N.*dN(rk+jl,pk,:);
            end
            if (r <= pk-1)
                Ntmp = 1/N(pk+1,rc);
                a_N = as(s1,kc-1)*Ntmp;
                as(s2,kc) = -a_N;
                das(s2,kc,:) = -(das(s1,kc-1,:) - a_N.*dN(pk+1,rc,:)).*Ntmp;
                Ntmp = N(rc,pk);
                d = d - a_N*Ntmp;
                dd = dd + Ntmp.*das(s2,kc,:) - a_N.*dN(rc,pk,:);
            end
            Ders(kc,rc) = d;
            dDers(kc,rc,:) = dd;
            j = s1;
            s1 = s2;
            s2 = j;
        end
    end
    r = p;
    for k = 1:ord
        Ders(k+1,:,:) = double(r).*Ders(k+1,:,:);
        dDers(k+1,:,:) = double(r).*dDers(k+1,:,:);
        r = r*(p-k);
    end
end
end


