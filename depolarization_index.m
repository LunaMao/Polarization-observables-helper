function [FD,PD,fd,pd] = depolarization_index(MM,mask)

a = size(MM);
row = a(1);
colum = a(2);


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



T_mb = (M00.^2+M01.^2+M02.^2+M03.^2) + (M10.^2+M11.^2+M12.^2+M13.^2) + (M20.^2+M21.^2+M22.^2+M23.^2) + (M30.^2+M31.^2+M32.^2+M33.^2);
FD = (4*M00.^2 - T_mb)/3;
PD = sqrt((1-FD./(M00.^2)));

fd=sum(sum(FD.*mask))/sum(sum(mask));

PD(isnan(PD))=0;
pd=sum(sum(PD.*mask))/sum(sum(mask));

% PD = PD*mask;
% FD =FD.*mask;
% PD=PD.*mask;
