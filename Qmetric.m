function [Q,q] = Qmetric(MM,PD,D,mask)
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

Q=(3*PD.^2-D.^2)./(1+D.^2);

Q(isnan(Q))=0;
q=sum(sum(Q))/sum(sum(mask));

