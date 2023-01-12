function [A,B] = make_array1(n)
A=zeros(n,n); % tworzenie nowej pustej macierzy
B=zeros(n,1);
for (i=1:n)
    for(j=1:n)
        if (j-i==-2 | j-i==2)
            A(i,j)=3-j/n; % wypełnianie macierzy danymi z treści zadania
        end
    end
    A(i,i)=-17;
    B(i,1)=2.5+0.5*i;
end
end