function [tt, X1, X2, T] = ode(xs1,xs2,b)
%
%   AUTOR
%       Tomasz Pawlak, 304104
%
%   CEL
%       Obliczanie trajektorii ruchu punktu na przedziale dla początkowych
%       współrzędnych na podstawie funkcji z funkcji zad1_func oraz metody
%       ode45 matlaba.
%
%   PARAMETRY WEJŚCIOWE
%       xs1,xs2 - współrzędne startowe układu (x1,x2)
%       b - koniec obserwowanego przedziału <0, b>  
%
%   PARAMETRY WYJŚCIA
%       tt  -   czas operacji
%       X1, X2 - wektory wartości funkcji w kolejnych krokach
%       T   -   wektor kolejnych kroków (jednostek czasu)
%
%   PRZYKŁADOWE WYWOŁANIE
%       >> [tto, X1, X2, T] = ode(0.001,-0.02,20)
%

% Czy chcemy wygenerować wykresy z tego pliku
PRINT=true;

tic;
% Przedział czasu
tspan=[0 b];
% Punkty startowe
x0=[xs1;xs2];
% Zadana tolerancja
opts = odeset('RelTol',1e-4,'AbsTol',1e-8);
% Generacja celu zadania
[t, x]=ode45(@zad1_func, tspan, x0, opts);

tt=toc;

% Generacja wykresów funkcji
if PRINT
    plot(t,x(:,1));
    title('ODE45 Trajektoria x1(t)');
    xlabel('Czas t');
    ylabel('Rozwiązanie x1');
    w=waitforbuttonpress;

    plot(t,x(:,2));
    title('ODE45 Trajektoria x2(t)');
    xlabel('Czas t');
    ylabel('Rozwiązanie x2');
    w=waitforbuttonpress;

    plot(x(:,1),x(:,2));
    title('ODE45 Trajektoria Przestrzeń Fazowa');
    xlabel('Rozwiązanie x1');
    ylabel('Rozwiązanie x2');
    w=waitforbuttonpress;
end
X1=x(:,1);
X2=x(:,2);
T=t;

end

function dxdt = zad1_func(t,x)
%
%   AUTOR
%       Tomasz Pawlak, 304104
%
%   CEL
%       Generacja funkcji dx1/dxt oraz dx2/dt dla metody ode45 zawartej w
%       pliku ode - liczącej trajektorie na podstawie wzorów.
%
%   PARAMETRY WEJŚCIOWE
%       t - czas, parametr niezależny
%       x(1), x(2) - położenie początkowe, parametr zależny
%
%   PARAMETRY WYJŚCIA
%       dxdt - macierz dwóch funkcji pochodnych określających położenie
%       punktów na wykresie
%
%   PRZYKŁADOWE WYWOŁANIE
%       zawarte w pliku ode.m, standardowa funkcja:
%       dxdt = [x(2)+x(1)*(0.3-(x(1))^2-(x(2))^2);...
%       -x(1)+x(2)*(0.3-(x(1))^2-(x(2))^2)];
%
dxdt = [x(2)+x(1)*(0.3-(x(1))^2-(x(2))^2);...
    -x(1)+x(2)*(0.3-(x(1))^2-(x(2))^2)];
end