function [A] = coef_matrix(n,X)
A=zeros(size(X,2),n+1);
for i=1:size(X,1)
 for j=1:n+1
 A(i,j)= X(i)^(n+1-j);
 end
end
end