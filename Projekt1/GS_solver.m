% Solver funkcji Gaussa-Seidel'a 
% x0    - aproksymacja wektora x - może być wektor l-b losowych
% tol   - oznaczone jako delta - tolerancja, kryterium stopu
% X     - najlepsze rozwiązanie
% eps_  - bład rozwiązania
% k     - liczba operacji poprawiania wyniku
function [x, eps_, k] = GS_solver(A,b,x0,tol) 
n=size(A,2);    % rozmiar macierzy
D=zeros(n,n);   % inicjacja pustych macierzy 
U=zeros(n,n);
L=zeros(n,n);
for(i=1:n) % generacja podziału A = L+D+U
    for(j=i+1:n)
        U(i,j)=A(i,j); 
    end
    D(i,i)=A(i,i);
end
L=A-U-D; % generacja macierzy L zgodnie ze wzorem 

Xj=x0; % alokacja pamięci na wektory X wejściowych iteracji (i)
Xi=zeros(n,1); % oraz następnych iteracji (i+1)
euk=norm(Xj,2); % norma euklidesowa dla Xj-Xi (poprawianie rozwiązania)
k=1;
while(euk > tol) % iter. poprawianie rozwiązań aż do osiągnięcia tolerancji
    Xi=Xj; % (i+1) rozwiązanie staje się (i) rozwiązaniem
    Xj=zeros(n,1);
    w=U*Xi-b;
    for(i=1:n) % wylicznie i-tej wart. nowego wektora rozwiązań
        for(j=1:i)
            Xj(i,1)=Xj(i,1)-(L(i,j)*Xj(j,1));
        end
        Xj(i,1)=(Xj(i,1)-w(i,1))/D(i,i);
    end
    euk=norm(Xj-Xi,2); % norma euklidesowa dla wektora xj-xi
    k=k+1;
end
eps_=norm(A*Xj-b,2);
x=Xj;
end
