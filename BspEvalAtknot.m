function pts = BspEvalAtknot(Gc_Crv,ord)
% 左闭右开区间
p = Gc_Crv.order-1;
U = Gc_Crv.knots;
ii = p-1; % basis function interval index
ich = 0; % control points index
us = Gc_Crv.knots(Gc_Crv.order:end-p);
lu = int16(length(us));
pts = zeros(Gc_Crv.dim,lu,ord+1,'like',Gc_Crv.knots(1));
for i = 1:lu-1
    ii = ii+1;
    ich = ich+1;
    N = BFDers(ii,us(i),p,U,ord);
    pts(:,i,:) = reshape(Gc_Crv.coefs(:,ich:ich+p)*N',Gc_Crv.dim,1,[]);
end
N = BFDers(ii,us(lu),p,U,ord);
pts(:,lu,:) = reshape(Gc_Crv.coefs(:,ich:ich+p)*N',Gc_Crv.dim,1,[]);
end