function [des] = LGPC(im, kps,eo, patch_size, s,o,M)

KPS=kps';
[yim,xim,~] = size(im);

CS = zeros(yim, xim, o);
for j=1:o
    for i=1:s
        CS(:,:,j)=CS(:,:,j)+abs(eo{i,j});
    end
end
[~, lgpc] = max(CS,[],3);

des = zeros(M*M*o, size(KPS,2));
kps_to_ignore = zeros(1,size(KPS,2));

for k = 1: size(KPS,2)
    x = round(KPS(1, k));
    y = round(KPS(2, k));
    
    x1 = max(1,x-floor(patch_size/2));
    y1 = max(1,y-floor(patch_size/2));
    x2 = min(x+floor(patch_size/2),size(im,2));
    y2 = min(y+floor(patch_size/2),size(im,1));
    
    if y2-y1 ~= patch_size || x2-x1 ~= patch_size
        kps_to_ignore(k)=1;
        continue;
    end
    
    patch = lgpc(y1:y2,x1:x2);
    [ys,xs]=size(patch);
    
    ns=M;
    LGPC_des = zeros(ns,ns,o);
    
    for j = 1:ns
        for i = 1:ns
            clip = patch(round((j-1)*ys/ns+1):round(j*ys/ns),round((i-1)*xs/ns+1):round(i*xs/ns));
            LGPC_des(j,i,:) = permute(hist(clip(:), 1:o), [1 3 2]);
        end
    end
    
    LGPC_des=LGPC_des(:);
    
    if norm(LGPC_des) ~= 0
        LGPC_des = LGPC_des /norm(LGPC_des);
    end

    
    des(1:M*M*o,k)=LGPC_des;
end
des = struct('kps', KPS(:,kps_to_ignore ==0)', 'des', des(:,kps_to_ignore==0)');


