function [P1,P2,P3,p1,p2,p3] = IPP(MM,mask)

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


pauli0 = ([1 0; 0 1]);
pauli1 = ([1 0; 0 -1]);
pauli2 = ([0 1; 1 0]);
pauli3 = ([0 -1i;1i 0]);

for ii=1:row
    for jj=1:colum

        h11 = M00(ii,jj)*kron(pauli0,conj(pauli0));
        h12 = M01(ii,jj)*kron(pauli0,conj(pauli1));
        h13 = M02(ii,jj)*kron(pauli0,conj(pauli2));
        h14 = M03(ii,jj)*kron(pauli0,conj(pauli3));

        h21 = M10(ii,jj)*kron(pauli1,conj(pauli0));
        h22 = M11(ii,jj)*kron(pauli1,conj(pauli1));
        h23 = M12(ii,jj)*kron(pauli1,conj(pauli2));
        h24 = M13(ii,jj)*kron(pauli1,conj(pauli3));

        h31 = M20(ii,jj)*kron(pauli2,conj(pauli0));
        h32 = M21(ii,jj)*kron(pauli2,conj(pauli1));
        h33 = M22(ii,jj)*kron(pauli2,conj(pauli2));
        h34 = M23(ii,jj)*kron(pauli2,conj(pauli3));

        h41 = M30(ii,jj)*kron(pauli3,conj(pauli0));
        h42 = M31(ii,jj)*kron(pauli3,conj(pauli1));
        h43 = M32(ii,jj)*kron(pauli3,conj(pauli2));
        h44 = M33(ii,jj)*kron(pauli3,conj(pauli3));


        % 
        % H(1,1) = M00(ii,jj)+ M01(ii,jj)+ M10(ii,jj)+ M11(ii,jj);
        % H(1,2) = M02(ii,jj)+ M12(ii,jj)+ 1i*(M03(ii,jj)+ M13(ii,jj));
        % H(1,3) = M20(ii,jj)+ M21(ii,jj)- 1i*(M30(ii,jj)+ M31(ii,jj));
        % H(1,4) = M22(ii,jj)+ M33(ii,jj)+ 1i*(M23(ii,jj)- M32(ii,jj));
        % 
        % H(2,1) = M02(ii,jj)+ M12(ii,jj)- 1i*(M03(ii,jj)- M13(ii,jj));
        % H(2,2) = M00(ii,jj)- M01(ii,jj)+ M10(ii,jj)- M11(ii,jj);
        % H(2,3) = M22(ii,jj)- M33(ii,jj)- 1i*(M23(ii,jj)+ M32(ii,jj));
        % H(2,4) = M20(ii,jj)- M21(ii,jj)- 1i*(M30(ii,jj)- M31(ii,jj));
        % 
        % H(3,1) = M20(ii,jj)+ M21(ii,jj)+ 1i*(M30(ii,jj)+ M31(ii,jj));
        % H(3,2) = M22(ii,jj)- M33(ii,jj)+ 1i*(M23(ii,jj)+ M32(ii,jj));
        % H(3,3) = M00(ii,jj)+ M01(ii,jj)- M10(ii,jj)+ M11(ii,jj);
        % H(3,4) = M02(ii,jj)- M12(ii,jj)+ 1i*(M03(ii,jj)- M13(ii,jj));
        % 
        % H(4,1) = M22(ii,jj)+ M33(ii,jj)- 1i*(M23(ii,jj)- M32(ii,jj));
        % H(4,2) = M20(ii,jj)- M21(ii,jj)+ 1i*(M30(ii,jj)- M31(ii,jj));
        % H(4,3) = M02(ii,jj)- M12(ii,jj)-1i*(M03(ii,jj)- M13(ii,jj));
        % H(4,4) = M00(ii,jj)- M01(ii,jj)+ M10(ii,jj) + M11(ii,jj);
        % 
        H = (h11+h12+h13+h14) + (h21+h22+h23+h24) + (h31+h32+h33+h34) + (h41+h42+h43+h44);
        H =H/4;
        
        eigenvalues = eig(H);
        [eigenvalues_sorted, Ori_indices] = sort(eigenvalues, 'descend');
        lambda0 = eigenvalues_sorted(1);
        lambda1 = eigenvalues_sorted(2);
        lambda2 = eigenvalues_sorted(3);
        lambda3 = eigenvalues_sorted(4);
        
        P1(ii,jj) = (1*(lambda0-lambda1))/M00(ii,jj);
        P2(ii,jj) = (1*(lambda0-lambda1)+2*(lambda1-lambda2))/M00(ii,jj);
        P3(ii,jj) = (1*(lambda0-lambda1)+2*(lambda1-lambda2)+3*(lambda2-lambda3))/M00(ii,jj);

    end
end

P1(isnan(P1))=0;
p1=sum(sum(P1))/(row*colum);

P2(isnan(P2))=0;
p2=sum(sum(P2))/(row*colum);

P3(isnan(P3))=0;
p3=sum(sum(P3))/(row*colum);

