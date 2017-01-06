function X = kronProd(A,B)
%KRON Kronecker product.
%   kron(A,B) returns the Kronecker product of two matrices A and B, of 
%   dimensions I-by-J and K-by-L respectively. The result is an I*K-by-J*L
%   block matrix in which the (i,j)-th block is defined as A(i,j)*B.

[I, J] = size(A);
[K, L] = size(B);

% Both matrices are dense.
A = reshape(A,[1 I 1 J]);
B = reshape(B,[K 1 L 1]);
X = reshape(bsxfun(@times,A,B),[I*K J*L]);

end