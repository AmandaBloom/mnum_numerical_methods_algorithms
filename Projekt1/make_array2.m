function [A,B] = make_array2(n)
A=zeros(n,n); % tworzenie nowej pustej macierzy
B=zeros(n,1);
for (i=1:n)
    for(j=1:n)
        A(i,j)=4*(i-j)+2; % wypełnianie macierzy danymi z treści zadania
    end
    A(i,i)=1/3;
    B(i,1)=3.5-0.4*i;
end
end