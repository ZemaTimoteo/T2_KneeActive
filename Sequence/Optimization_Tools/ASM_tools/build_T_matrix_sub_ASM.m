% reduced version
function T =  build_T_matrix_sub_ASM(AA,nn)
T=zeros([nn nn]);

ind = 1:(nn/3);
T(3*nn*(ind-1)+3*(ind-1)+1)  = AA(1);
T(3*nn*(ind-1)+3*(ind-1)+2)  = AA(2);
T(3*nn*(ind-1)+3*(ind-1)+3)  = AA(3);
T(3*nn*ind-2*nn+3*(ind-1)+1) = AA(4);
T(3*nn*ind-2*nn+3*(ind-1)+2) = AA(5);
T(3*nn*ind-2*nn+3*(ind-1)+3) = AA(6);
T(3*nn*ind-nn+3*(ind-1)+1)   = AA(7);
T(3*nn*ind-nn+3*(ind-1)+2)   = AA(8);
T(3*nn*ind-nn+3*(ind-1)+3)   = AA(9);

T = sparse(T);
%        T = sparse(kron(eye(nn/3),AA));% slow
end