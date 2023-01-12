function [tt, X1, X2, T, err, hh] = rk4z(dx1,dx2,xs1,xs2,h0,b)
%
%   AUTOR
%       Tomasz Pawlak, 304104
%
%   CEL
%       Obliczanie trajektorii ruchu punktu na przedziale dla początkowych
%       współrzędnych wykorzystując 
%
%   PARAMETRY WEJŚCIOWE
%       dx1,dx2 -   równania różniczkowe opisujące ruch
%       xs1,xs2 -   współrzędne startowe układu (x1,x2)
%       h0      -   krok początkowy
%       b       -   koniec obserwowanego przedziału <0, b>
%   PARAMETRY WYJŚCIA
%       tt      -   czas operacji
%       X1, X2  -   wektory wartości funkcji w kolejnych krokach
%       T       -   wektor wartości zmiennej niezależnej między iteracjami
%       err     -   wektor najbardziej znaczącego błędu w każdej iteracji
%       hh      -   wektor długości kolejnych kroków
%   PRZYKŁADOWE WYWOŁANIE
%       >> [ttr,X1,X2,T,err,hh]=rk4z(@(x1,x2)x2+x1*(0.3-(x1)^2-(x2)^2),@(x1,x2)-x1+x2*(0.3-(x1)^2-(x2)^2),0.001,-0.02,0.01,20)
%
PRINT=true;
tic;

% Inicjacja kroku, stałych korekty długości kroku oraz dokładności i
% minimlanego kroku
h=h0;
hh=h;
s = 0.9;
beta_ = 5;
E_n = 1e-4;
E_b = 1e-8;
hmin = 1e-12;
%RelTol',E_n,'AbsTol',E_b;

% Inicjacja wartości zmiennej niezależnej
t=h;

% Inicjacja wektorów wartości zmiennych zależnych oraz wyniku
X1=xs1;
X2=xs2;
T=0;

% Początkowe wartości kroku oraz półkroku
x1=xs1;
x2=xs2;
x1h=xs1;
x2h=xs2;
x1d=xs1;
x2d=xs2;

% Inicjacja wektora błędów
err=0;

% Wskaźnik iteracji
n=0;
while true
    n=n+1;
    if(mod(n,2) == 1)
        % Równania półkroku 1
        k_1_1d=dx1(x1,x2);
        k_1_2d=dx2(x1,x2);
        k_2_1d=dx1((x1+0.25*h*k_1_1d),(x2+0.25*h*k_1_2d));
        k_2_2d=dx2((x1+0.25*h*k_1_1d),(x2+0.25*h*k_1_2d));
        k_3_1d=dx1((x1+0.25*h*k_2_1d),(x2+0.25*h*k_2_2d));
        k_3_2d=dx2((x1+0.25*h*k_2_1d),(x2+0.25*h*k_2_2d));
        k_4_1d=dx1((x1+0.5*h*k_3_1d), (x2+0.5*h*k_3_2d));
        k_4_2d=dx2((x1+0.5*h*k_3_1d), (x2+0.5*h*k_3_2d));
    
        x1d=x1+(1/12)*h*(k_1_1d+2*k_2_1d+2*k_3_1d+k_4_1d);
        x2d=x2+(1/12)*h*(k_1_2d+2*k_2_2d+2*k_3_2d+k_4_2d);
        
        % Równania półkroku 2
        k_1_1h=dx1(x1d,x2h);
        k_1_2h=dx2(x1d,x2h);
        k_2_1h=dx1((x1d+0.25*h*k_1_1h),(x2h+0.25*h*k_1_2h));
        k_2_2h=dx2((x1d+0.25*h*k_1_1h),(x2h+0.25*h*k_1_2h));
        k_3_1h=dx1((x1d+0.25*h*k_2_1h),(x2h+0.25*h*k_2_2h));
        k_3_2h=dx2((x1d+0.25*h*k_2_1h),(x2h+0.25*h*k_2_2h));
        k_4_1h=dx1((x1d+0.5*h*k_3_1h), (x2h+0.5*h*k_3_2h));
        k_4_2h=dx2((x1d+0.5*h*k_3_1h), (x2h+0.5*h*k_3_2h));
    
        x1h=x1d+(1/12)*h*(k_1_1h+2*k_2_1h+2*k_3_1h+k_4_1h);
        x2h=x2d+(1/12)*h*(k_1_2h+2*k_2_2h+2*k_3_2h+k_4_2h);
    end
    if(mod(n,2) == 0)
        % Równania kroku całkowitego
        k_1_1=dx1(x1,x2);
        k_1_2=dx2(x1,x2);
        k_2_1=dx1((x1+0.5*h*k_1_1),(x2+0.5*h*k_1_2));
        k_2_2=dx2((x1+0.5*h*k_1_1),(x2+0.5*h*k_1_2));
        k_3_1=dx1((x1+0.5*h*k_2_1),(x2+0.5*h*k_2_2));
        k_3_2=dx2((x1+0.5*h*k_2_1),(x2+0.5*h*k_2_2));
        k_4_1=dx1((x1+h*k_3_1),(x2+h*k_3_2));
        k_4_2=dx2((x1+h*k_3_1),(x2+h*k_3_2));

        x1=x1+(1/6)*h*(k_1_1+2*k_2_1+2*k_3_1+k_4_1);
        x2=x2+(1/6)*h*(k_1_2+2*k_2_2+2*k_3_2+k_4_2);
        X1=[X1; x1h];
        X2=[X2; x2h];
        
        T=[T; t];
        % Obliczanie błędów na podstawie kroku
        
        delta1=(1/15)*abs(x1h-x1);
        delta2=(1/15)*abs(x2h-x2);
        epse1=abs(x1h)*E_n + E_b;
        epse2=abs(x2h)*E_n + E_b;

        err=[err; min(epse1, epse2)];
        
        alfa=min(epse1/delta1, epse2/delta2);
        alfa=alfa^(1/5);

        % Ustalenie nowej wartości kroku dla kolejnej pętli
        h_new = s*alfa*h;

        if s*alfa >= 1
            % Następny Punkt
            h=min([h_new, beta_*h, b-t]);
            t=t+h;
                
        elseif h_new >= hmin
            % Nowa wartość kroku
            h=h_new;
            n=n-2;
            % Zakładamy powtórzenie kroku
            err = err(1:length(err)-1);
            X1 = X1(1:length(X1)-1);
            X2 = X2(1:length(X2)-1);
            x1 = X1(length(X1));
            x2 = X2(length(X2));
            T = T(1:length(T)-1);
            hh = hh(1:length(hh)-1);
        else
            % Przypadek niemożliwy do rozwiązania - za duże E
            error('Niemożliwe rozwiązanie z zadaną dokładnością');
        end
        hh=[hh; h];
    end
    if t+h >= b
        break
    end
end
tt = toc;

if PRINT
    plot(T,X1);
    title('RK4 Trajektoria x1(t)');
    xlabel('Czas t');
    ylabel('Rozwiązanie x1');
    w=waitforbuttonpress;

    plot(T,X2);
    title('RK4 Trajektoria x2(t)');
    xlabel('Czas t');
    ylabel('Rozwiązanie x2');
    w=waitforbuttonpress;

    plot(X1, X2);
    title('RK4 Trajektoria Przestrzeń Fazowa');
    xlabel('Rozwiązanie x1');
    ylabel('Rozwiązanie x2');
    w=waitforbuttonpress;
    clf;
    close all;
end

