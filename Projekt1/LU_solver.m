function [X, L, U, b, Prow, Pcol] = LU_solver(A, b)
% L,U - macierze trójkątne rozkładu LU 
% X - wektor rozwiązań
n = size(A, 2);  % pobieranie rozmiaru macierzy
L = zeros(n, n); % inicjowanie wartości macierzy trójkątnych
Prow = eye(n);   % macierz przeksztłceń wierszy
Pcol = eye(n);   % macierz przekształceń kolumn

for(i=1:n)
    maxA=abs(A(i,i)); % pobranie modułu wartości 1-ego elementu podmacierzy
    maxrow=i; % pobranie numeru 1-ej kolumny i wiersza
    maxcol=i; 
    for(i_n=i:n) % wyznaczanie elemantu o największym module
        for(j_n=i:n)
            if(abs(A(i_n,j_n))>maxA)
                maxA=abs(A(i_n,j_n));
                maxrow=i_n;
                maxcol=j_n;
            end
        end
    end

    A(:,[i maxcol]) = A(:,[maxcol i]); % podmiana wierszy i kolumn macierzy
    A([i maxrow],:) = A([maxrow i],:);

    b([i maxrow],:) = b([maxrow i],:);
  
    L(:,[i maxcol]) = L(:,[maxcol i]);
    L([i maxrow],:) = L([maxrow i],:);

    Pcol(:,[i maxcol]) = Pcol(:,[maxcol i]); 
    Prow([i maxrow],:) = Prow([maxrow i],:);

    for(i_n=i+1:n)
        if(A(i,i)~=0)
            L(i_n,i) = A(i_n,i)/A(i,i); % parametry macierzy trójkątnej L
            A(i_n,:) = A(i_n,:)-(L(i_n,i)*A(i,:)); % w_j = w_j - l_ji*w_i
            b(i_n,:) = b(i_n,:)-(L(i_n,i)*b(i,:)); 
        end
    end
end
U = A;  % przekształcenia A są równe U

X = zeros(n,1); % inicjacja wektora rozwiązań
for (i=n:-1:1) % wyznaczanie wektora X
    for (j=n:-1:i+1)
        X(i,1)=X(i,1)-(U(i,j)*X(j,1));
    end
    X(i,1)=(X(i,1)+b(i,1))/U(i,i);
end
X = Pcol*X; % przekształcenie wektora rozwiązań (poprawna kolejność)

L = L + eye(n); % dodanie do L macierzy jednostkowej (zgodność L*U = P*A)
end