function [CP_Ref,CP_Sen] = CMPC_match(m1,cmpc1,cmpc2,H,templateSize)

[im_RefH,im_RefW] = size(m1);

templateRad = round(templateSize/2);    
marg=templateRad+2;

matchRad = templateRad;

im1 = m1(marg:im_RefH-marg,marg:im_RefW-marg);
m1_points = detectFASTFeatures(im1,'MinContrast',0.05);
m1_points = m1_points.selectStrongest(5000); 
points1 =[m1_points.Location(:,2),m1_points.Location(:,1)] + marg - 1;
pNum = size(points1,1); 

C = 1;

for n = 1: pNum
    X_Ref=points1(n,2);
    Y_Ref=points1(n,1);
    
    tempCo = [X_Ref;Y_Ref;1];
    tempCo1 = H*tempCo;
    
    X_Sen_c1 = round(tempCo1(1));
    Y_Sen_c1 = round(tempCo1(2)); 

    if (X_Sen_c1 < marg+1 || X_Sen_c1 > size(cmpc2(:,:,1),2)-marg || Y_Sen_c1<marg+1 || Y_Sen_c1 > size(cmpc2(:,:,1),1)-marg)
        continue;
    end            
    
    featureSub_Ref = single(cmpc1(Y_Ref-matchRad:Y_Ref+matchRad,X_Ref-matchRad:X_Ref+matchRad,:));
    featureSub_Sen = single(cmpc2(Y_Sen_c1-matchRad:Y_Sen_c1+matchRad,X_Sen_c1-matchRad:X_Sen_c1+matchRad,:));

    [max_i, max_j] = PhaseCorrelation3D(featureSub_Ref, featureSub_Sen, matchRad);
    [a,~] = size(max_i);
    if a>1
        continue;
    end
     
    Y_match = Y_Sen_c1 + max_i;
    X_match = X_Sen_c1 + max_j;
    
    C = C+1;
    CP_Ref(C,:) = [X_Ref,Y_Ref];
    CP_Sen(C,:) = [X_match,Y_match];
end




