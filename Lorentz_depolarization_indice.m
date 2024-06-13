function [L1,L2,l1,l2] = Lorentz_depolarization_indice(MM,mask)

a = size(MM);
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



G = diag([1 -1 -1 -1]);
for i = 1 : row
    for j = 1 : colum
        M = squeeze(MM(i,j,:,:));
        MT = M';
        N=G*MT*M*G;
        eigenvalues = eig(N);
        [eigenvalues_sorted, Ori_indices] = sort(eigenvalues, 'descend');
        rho_1 = eigenvalues_sorted(1);
        rho_2 = eigenvalues_sorted(2);
        rho_3 = eigenvalues_sorted(3);
        rho_4 = eigenvalues_sorted(4);
        
        L1(i,j)=sqrt((rho_2+rho_3+rho_4)/(3*rho_1));
        tup = 4*(rho_1^2+rho_2^2+rho_3^2+rho_4^2)-(rho_1+rho_2+rho_3+rho_4)^2;
        tdown = 3*(rho_1+rho_2+rho_3+rho_4)^2;
        L2(i,j) = sqrt(tup/tdown);

        % L2(i,j)=sqrt((4*(rho_1^2+rho_2^2+rho_3^2+rho_4^2)-(rho_1+rho_2+rho_3+rho_4)^2))/(3*(rho_1+rho_2+rho_3+rho_4)^2);        
    end
end

L1(isnan(L1))=0;
l1=sum(sum(L1))/sum(sum(mask));

L2(isnan(L2))=0;
l2=sum(sum(L2))/sum(sum(mask));
