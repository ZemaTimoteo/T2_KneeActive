
function w = blockdiag_mult_ASM(A,v)
    %given A square matrix, it returns   ' w =  kron(A,eye(nA))*v  '   efficiently
    nA = size(A,1); %square
    V = reshape(v,nA,length(v)/nA);
    W = A*V;
    w = W(:);
end