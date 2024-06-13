function [P,D,Ps,p,d,ps] = CP(MM,mask)

a= size(MM);
row = a(1);
colum = a(2);

M00 = squeeze(MM(:,:,1,1));
M01 = squeeze(MM(:,:,1,2));
M02 = squeeze(MM(:,:,1,3));
M03 = squeeze(MM(:,:,1,4));

M10 = squeeze(MM(:,:,2,1));
M11 = squeeze(MM(:,:,2,2));
M12 = squeeze(MM(:,:,2,3));
M13 = squeeze(MM(:,:,2,4));

M20 = squeeze(MM(:,:,3,1));
M21 = squeeze(MM(:,:,3,2));
M22 = squeeze(MM(:,:,3,3));
M23 = squeeze(MM(:,:,3,4));

M30 = squeeze(MM(:,:,4,1));
M31 = squeeze(MM(:,:,4,2));
M32 = squeeze(MM(:,:,4,3));
M33 = squeeze(MM(:,:,4,4));

M00 = squeeze(MM(:,:,1,1));
M01 = squeeze(MM(:,:,1,2))./M00;
M02 = squeeze(MM(:,:,1,3))./M00;
M03 = squeeze(MM(:,:,1,4))./M00;

M10 = squeeze(MM(:,:,2,1))./M00;
M11 = squeeze(MM(:,:,2,2))./M00;
M12 = squeeze(MM(:,:,2,3))./M00;
M13 = squeeze(MM(:,:,2,4))./M00;

M20 = squeeze(MM(:,:,3,1))./M00;
M21 = squeeze(MM(:,:,3,2))./M00;
M22 = squeeze(MM(:,:,3,3))./M00;
M23 = squeeze(MM(:,:,3,4))./M00;

M30 = squeeze(MM(:,:,4,1))./M00;
M31 = squeeze(MM(:,:,4,2))./M00;
M32 = squeeze(MM(:,:,4,3))./M00;
M33 = squeeze(MM(:,:,4,4))./M00;

M00 = M00./M00;


for ii=1:row
    for jj=1:colum

        P(ii,jj) = sqrt((M10(ii,jj)^2+M20(ii,jj)^2+M30(ii,jj)^2)/M00(ii,jj)^2);
        D(ii,jj) = sqrt((M01(ii,jj)^2+M02(ii,jj)^2+M03(ii,jj)^2)/M00(ii,jj)^2);
        Ps(ii,jj) = sqrt((M11(ii,jj)^2+M12(ii,jj)^2+M13(ii,jj)^2+M21(ii,jj)^2+M22(ii,jj)^2+M23(ii,jj)^2+M31(ii,jj)^2+M32(ii,jj)^2+M33(ii,jj)^2)/(3*M00(ii,jj)^2));

    end
end

P(isnan(P))=0;
p=sum(sum(P.*mask))/sum(sum(mask));

D(isnan(D))=0;
d=sum(sum(D.*mask))/sum(sum(mask));

Ps(isnan(Ps))=0;
ps=sum(sum(Ps.*mask))/sum(sum(mask));



