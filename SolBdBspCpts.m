% function pts = SolBdBspCpts(pvs,ukrefs,scales)
function pts = SolBdBspCpts(pvs,ukrefs)
order = int8(size(ukrefs,1)); 
dim = int8(size(pvs,1)); 
num = int8(size(pvs,2));
% pts = [pvs(:,:,1) zeros(length(pvs),order-1)];
pts = zeros(dim,order,num);
pts(:,1,:) = pvs(:,:,1);
for i=2:order
    idxtmp = order-i+2;
%     pvs(:,:,idxtmp) = pvs(:,:,idxtmp)./(scales.^(idxtmp-1));
    pts(:,idxtmp:end,:) = reshape(ukrefs(2:i,:)-ukrefs(1,:),1,[],num)...
        .*(pts(:,idxtmp:end,:) + reshape(pvs(:,:,idxtmp),dim,1,[]))/double(i);
end
for i = 2:order
    pts(:,i,:) = pts(:,i-1,:) + pts(:,i,:);
end
pts = pts(:,2:end,:);
end