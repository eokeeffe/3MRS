% This is a fully implementation of the proposed 3MRS algorithm. 

clc;clear;close all;
warning('off')
% read two images 
file_image='.\CoFSM\';
[filename,pathname]=uigetfile({'*.*','All Files(*.*)'},'Select reference image',...
                          file_image);
im1=imread(strcat(pathname,filename));

[filename,pathname]=uigetfile({'*.*','All Files(*.*)'},'Select the sensed image',...
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

%parameters set
n=4;
k=5;
t=0.5;
N=96;
M=8;
d=3;
L=96;

tic;
%%%%    Coarse matching    %%%%
disp('Coarse matching')
% feature detection and description
[des_m1,des_m2,cmpc1,cmpc2,M1] = FaetureCreation(im1,im2,n,k,t,N,M,d);

% nearest matching
[indexPairs,matchmetric] = matchFeatures(des_m1.des,des_m2.des,'MaxRatio',1,'MatchThreshold', 100);
matchedPoints1 = des_m1.kps(indexPairs(:, 1), :);
matchedPoints2 = des_m2.kps(indexPairs(:, 2), :);
[matchedPoints2,IA]=unique(matchedPoints2,'rows');
matchedPoints1=matchedPoints1(IA,:);

%outlier removal
H=FSC(matchedPoints1,matchedPoints2,'affine',2);
Y_=H*[matchedPoints1';ones(1,size(matchedPoints1,1))];
Y_(1,:)=Y_(1,:)./Y_(3,:);
Y_(2,:)=Y_(2,:)./Y_(3,:);
E=sqrt(sum((Y_(1:2,:)-matchedPoints2').^2));
inliersIndex=E<3;
cleanedPoints1 = matchedPoints1(inliersIndex, :);
cleanedPoints2 = matchedPoints2(inliersIndex, :);

%%%%    Fine matching    %%%%
disp('Fine matching')
%template matching
[matchedPoints1,matchedPoints2] = CMPC_match(M1,cmpc1,cmpc2,H,L);

%outlier removal
H=FSC(matchedPoints1,matchedPoints2,'affine',2);
Y_=H*[matchedPoints1';ones(1,size(matchedPoints1,1))];
Y_(1,:)=Y_(1,:)./Y_(3,:);
Y_(2,:)=Y_(2,:)./Y_(3,:);
E=sqrt(sum((Y_(1:2,:)-matchedPoints2').^2));
inliersIndex=E<1;
cleanedPoints3 = matchedPoints1(inliersIndex, :);
cleanedPoints4 = matchedPoints2(inliersIndex, :);

fprintf('t:%.3fs\n',toc);
[NCM,~] = size(cleanedPoints3);
fprintf('NCM:%d\n',NCM);

disp('Show matches')
% Show results
figure; showMatchedFeatures(im1, im2, cleanedPoints3, cleanedPoints4, 'montage');

disp('registration result')
% registration
image_fusion(im2,im1,double(H));


