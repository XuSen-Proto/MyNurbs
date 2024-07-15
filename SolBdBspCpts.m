% function pts = SolBdBspCpts(pvs,ukrefs,scales)
function pts = SolBdBspCpts(pvs,ukrefs)
% pvs: dim *num*order
% pts: dim*order*num
order = int8(size(ukrefs,1)); 
dim = int8(size(pvs,1)); 
num = int8(size(pvs,2));
% pts = [pvs(:,:,1) zeros(length(pvs),order-1)];
pts = zeros(dim,order,num,'like',pvs(1));
pts(:,1,:) = pvs(:,:,1);
for i=2:order
    idxtmp = order-i+2;

    % 用于测试
    % disp('a')
    % size(reshape(ukrefs(2:i,:)-ukrefs(1,:),1,[],num))
    % disp('b')
    % size(pts(:,idxtmp:order,:))
    % disp('c')
    % size(reshape(pvs(:,:,idxtmp),dim,1,[]))
%     pvs(:,:,idxtmp) = pvs(:,:,idxtmp)./(scales.^(idxtmp-1));

    %用于数值运算
    % pts(:,idxtmp:order,:) = reshape(ukrefs(2:i,:)-ukrefs(1,:),1,[],num)...
    %     .*(pts(:,idxtmp:order,:) + reshape(pvs(:,:,idxtmp),dim,1,[]))/double(i);

    %用于符号运算
    utmp = reshape(ukrefs(2:i,:)-ukrefs(1,:),1,[],num);
    atmp = reshape(pvs(:,:,idxtmp),dim,1,[]);
    for k = idxtmp:order
        pts(:,k,:) = utmp(:,k-idxtmp+1,:).*(pts(:,k,:)+atmp)/double(i);
    end
end

for i = 2:order
    pts(:,i,:) = pts(:,i-1,:) + pts(:,i,:);
end

pts = pts(:,2:end,:);
end