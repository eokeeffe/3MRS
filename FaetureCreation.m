function [des_m1,des_m2,cmpc1,cmpc2,M1] = FaetureCreation(im1,im2,s,o,t,patch_size,M,d)

[M1,m1,~,~,~,eo1,~,~] = phasecong3(im1,s,o,3,'mult',1.6,'sigmaOnf',0.75,'g', 3, 'k',1);
[M2,m2,~,~,~,eo2,~,~] = phasecong3(im2,s,o,3,'mult',1.6,'sigmaOnf',0.75,'g', 3, 'k',1);

%PC-FAST detector 
M1=t*M1+(t-1)*m1;M2=t*M2+(t-1)*m2;
a=max(M1(:)); b=min(M1(:)); M1=(M1-b)/(a-b);
a=max(M2(:)); b=min(M2(:)); M2=(M2-b)/(a-b);

[M1H,M1W] = size(M1);
marg= round(patch_size/2)+2;          
M11 = M1(marg:M1H-marg,marg:M1W-marg);
m1_points = detectFASTFeatures(M11,'MinContrast',0.05);
m1_points = m1_points.selectStrongest(5000); 
points1 =m1_points.Location + marg - 1;

[M2H,M2W] = size(M2);
marg= round(patch_size/2)+2;          
M22 = M2(marg:M2H-marg,marg:M2W-marg);
m2_points = detectFASTFeatures(M22,'MinContrast',0.05);
m2_points = m2_points.selectStrongest(5000); 
points2 =m2_points.Location + marg - 1;

%LGPC feature descriptor
[des_m1] = LGPC(im1, points1,eo1, patch_size, s,o,M);
[des_m2] = LGPC(im2, points2,eo2, patch_size, s,o,M);

%CMPC template feature
cmpc1 = CMPC(eo1,o,s,d);
cmpc2 = CMPC(eo2,o,s,d);






