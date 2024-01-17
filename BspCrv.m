classdef BspCrv < handle
    properties
        knots; coefs; order; dim;
        kid; coefsc; 
        Basis; GetSpan;
    end
    methods
        function obj = BspCrv(knots,coefs)
            obj.knots = knots;
            obj.coefs = coefs;
            obj.dim = size(coefs,1); nc = size(coefs,2);
            order = size(knots,2)-nc; obj.order = order;
            obj.Basis = @(i,u)basisfun(i,u,order-1,knots);
            obj.GetSpan = @(u)findspan(nc-1,order-1,u,knots);
            kid = obj.GetSpan(0); obj.kid = kid(1); obj.Select(kid(1));
        end
        
        function Select(obj,k)
            tmp = k-obj.order+2;
            obj.coefsc = obj.coefs(:,tmp:tmp+obj.order-1);
        end
        
        function [val,k] = Eval(obj,u)
            ki = obj.GetSpan(u);
            k = ki(1); obj.Select(k);
            N = obj.Basis(k,u);
            val = obj.coefsc*N';
        end

        function val = Eval_km(obj,u)
            if u>obj.knots(obj.kid+2)
                iadd = 1; itmp = obj.kid+2;
                while u>obj.knots(itmp+iadd)
                    iadd = iadd+1;
                end
                obj.kid = obj.kid+iadd;
                obj.Select(obj.kid);
            end
            N = obj.Basis(obj.kid,u);
            val = obj.coefsc*N';
        end
        
        function vals = Eval_us(obj,us)
            len = length(us);
            vals = zeros(obj.dim,len);
%             [vals(:,1),obj.kid] = obj.Eval(us(1));
            for i = 1:len
                vals(:,i) = obj.Eval_km(us(i));
            end
        end
        
        function crv = Val(obj)
            crv = nrbmak(obj.coefs,obj.knots);
        end
    end
end