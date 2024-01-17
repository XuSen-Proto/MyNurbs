function pts = SolBdBspCptsS(pvs,ukrefs,scales)
p = size(ukrefs,1); dim = size(pvs,1); num = size(pvs,2);
% pts = [pvs(:,:,1) zeros(length(pvs),order-1)];
pts = zeros(dim,p,num);
pts(:,1,:) = pvs(:,:,1);
ur = reshape(ukrefs(2:end,:)-ones(p-1,1)*ukrefs(1,:),[],1,num);
Q = zeros(dim*(p-2),p-2,num);
sclflag = ~isempty(scales);
for k = 1:num
    qtmp = zeros(dim,p);
    for i=2:p
        coeftmp0 = bsxfun(@rdivide,ur(1:i-1,:,k),(p-i+2:p)'); coeftmp = coeftmp0;
        for j =i-2:-1:1
            coeftmp(j,:) = bsxfun(@times,coeftmp(j,:),coeftmp(j+1,:));
        end
        if sclflag
            pvs(:,k,i:end) = bsxfun(@rdivide,pvs(:,k,i:end),scales);
        end
        intmp = reshape(pvs(:,k,i:-1:max([i-1 2])),3,[],1);
        if size(intmp,2)==size(coeftmp,1)
            addtmp = intmp*coeftmp;
        else
            addtmp = [intmp qtmp(:,1:size(coeftmp,1)-size(intmp,2))]*coeftmp;
        end
        pts(:,i,k) = pts(:,i-1,k) + addtmp;
        if i>2
            Q(dim*(i-3)+(1:dim),end,k) = bsxfun(@rdivide,addtmp,coeftmp(end,:));
            for j = 1:i-3
                Q(dim*(i-3)+(1:dim),end-j,k) = bsxfun(@rdivide,(Q(dim*(i-3)+(1:dim),end-j+1,k)-...
                    Q(dim*(i-4)+(1:dim),end-j+1,k)),coeftmp0(end-j+1,:));
            end
            % for codegen
            qtmp(:,1:i-2) = Q(dim*(i-3)+(1:dim),end-i+3:end,k);
        end
    %     scales = scales.*scales;
    end
end
pts = pts(:,2:end,:);
end