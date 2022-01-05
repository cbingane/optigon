function [B,C] = difference_of_psd(A)
if ishermitian(A)
    [V,D] = eig(A);
    D = sparse(D);
    P = (D+abs(D))/2;
    B = V*P*V'; B = (B+B')/2;
    C = B-A;
end
end