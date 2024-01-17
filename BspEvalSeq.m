function pts = BspEvalSeq(us,Gc_Crv,ord)
lu = int16(length(us));
p = Gc_Crv.order-1;
U = Gc_Crv.knots;
ii = p; % basis function interval index
iu = p+2; % interval right value index
ich = 1; % control points index
cpts = Gc_Crv.coefs(:,ich:ich+p);
pts = zeros(Gc_Crv.dim,lu,ord+1);
for i = 1:lu
    utmp = us(i);
    if utmp>U(iu)
        while utmp>U(iu)
            ii = ii+1;
            iu = iu+1;
            ich = ich+1;
        end
        cpts = Gc_Crv.coefs(:,ich:ich+p);
    end
    N = BFDers(ii,utmp,p,U,ord);
    pts(:,i,:) = reshape(cpts*N',Gc_Crv.dim,1,[]);
end
end