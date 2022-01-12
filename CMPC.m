function [cmpc] = CMPC(eo,o,s,d)

[r,c,~] = size(eo{1,1});
temp = zeros(r, c, o);
for j=1:o
    for i=1:s
        temp(:,:,j)=temp(:,:,j)+abs(eo{i,j});
    end
end

g = fspecial('gaussian',d,d/6);

for i = 1:o
    temp(:,:,i) = conv2(temp(:,:,i), g, 'same'); 
end

temp_array = zeros(r,c,o+2);
for i = 2:o+1
    temp_array(:, :, i) = temp(:, :, i - 1);
end

delta = 0.0000000001;
denominator_normalized = zeros(r,c,1) + delta;

cmpc = zeros(r,c,o);
for i = 1:o
    cmpc(:, :, i) = temp_array(:, :, i) + 3 * temp_array(:, :, i + 1) + temp_array(:, :, i + 2);
    denominator_normalized = denominator_normalized + cmpc(:, :, i).^(2);
end
denominator_normalized = denominator_normalized.^(0.5);

for i = 1:o
    cmpc(:, : , i) = cmpc(:, : , i)./denominator_normalized;
end