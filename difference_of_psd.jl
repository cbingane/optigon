function difference_of_psd(A)
if ishermitian(A)
    S = eigen(A);
    D = S.values;
    V = S.vectors;
    P = (D+abs.(D))/2;
    B = V*Diagonal(P)*V'; B = (B+B')/2;
    C = B-A;
    return(B,C)
end
end