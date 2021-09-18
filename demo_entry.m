clc;
clear;
close all;
warning('off');

% read two images 
file_image='.\CoFSM\';
[filename,pathname]=uigetfile({'*.*','All Files(*.*)'},'Select reference image',...
                          file_image);
im1=imread(strcat(pathname,filename));
[filename,pathname]=uigetfile({'*.*','All Files(*.*)'},'Select the image to be registered',...
                          file_image);
im2=imread(strcat(pathname,filename));

if size(im1,3)==1
    temp=im1;
    im1(:,:,1)=temp;
    im1(:,:,2)=temp;
    im1(:,:,3)=temp;
end

if size(im2,3)==1
    temp=im2;
    im2(:,:,1)=temp;
    im2(:,:,2)=temp;
    im2(:,:,3)=temp;
end

tic;
disp('RIFT coarse matching')
% RIFT feature detection and description
[des_m1,des_m2,eo1,eo2,m1] = RIFT_no_rotation_invariance(im1,im2,4,6,72);

% nearest matching
[indexPairs,matchmetric] = matchFeatures(des_m1.des,des_m2.des,'MaxRatio',1,'MatchThreshold', 100);
matchedPoints1 = des_m1.kps(indexPairs(:, 1), :);
matchedPoints2 = des_m2.kps(indexPairs(:, 2), :);
[matchedPoints2,IA]=unique(matchedPoints2,'rows');
matchedPoints1=matchedPoints1(IA,:);

% FSC outlier removal
H=FSC(matchedPoints1,matchedPoints2,'affine',2);
Y_=H*[matchedPoints1';ones(1,size(matchedPoints1,1))];
Y_(1,:)=Y_(1,:)./Y_(3,:);
Y_(2,:)=Y_(2,:)./Y_(3,:);
E=sqrt(sum((Y_(1:2,:)-matchedPoints2').^2));
inliersIndex=E<3;
cleanedPoints1 = matchedPoints1(inliersIndex, :);
cleanedPoints2 = matchedPoints2(inliersIndex, :);

disp('CMPC fine matching')
% reuse log-Gabor convolution sequences to construct
% Channel Maps of Phase Congruency (CMPC)
CMPC1 = CMPC(eo1,4,6);
CMPC2 = CMPC(eo2,4,6);

% template matching with CMPCs as inputs and Three-dimensional phase
% correlation as similarity measure
[matchedPoints1,matchedPoints2] = CMPC_match(m1,CMPC1,CMPC2,H,100);

% FSC outlier removal
H=FSC(matchedPoints1,matchedPoints2,'affine',2);
Y_=H*[matchedPoints1';ones(1,size(matchedPoints1,1))];
Y_(1,:)=Y_(1,:)./Y_(3,:);
Y_(2,:)=Y_(2,:)./Y_(3,:);
E=sqrt(sum((Y_(1:2,:)-matchedPoints2').^2));
inliersIndex=E<1;
cleanedPoints3 = matchedPoints1(inliersIndex, :);
cleanedPoints4 = matchedPoints2(inliersIndex, :);

% print and show the results
[NCM,~] = size(cleanedPoints3);
fprintf('number of correct matchs:%d\n',NCM);
fprintf('running time:%.3fs\n',toc);

disp('Show matches')
figure;
imshow(im1),hold on;
plot(cleanedPoints3(:,1),cleanedPoints3(:,2),'+','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10);hold on;
title('reference image');
figure;
imshow(im2),hold on;
plot(cleanedPoints4(:,1),cleanedPoints4(:,2),'+','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',10);hold on;
title('sensed image');

disp('Show registration and fusion result')
% registration
image_fusion(im2,im1,double(H));

