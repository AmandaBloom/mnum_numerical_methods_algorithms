function [] = comparison(xs1,xs2,h0,b)
%
%   AUTOR
%       Tomasz Pawlak, 304104
%
%   CEL
%       Porównanie solverów matlaba ode45 z pliku ode oraz autorskiego
%       solwera z pliku rk4.m w obliczaniu trajektorii ruchu punktu na 
%       przedziale dla początkowych współrzędnych
%
%   PARAMETRY WEJŚCIOWE
%       dx1,dx2 -   równania różniczkowe opisujące ruch
%       xs1,xs2 -   współrzędne startowe układu (x1,x2)
%       h0      -   krok początkowy dla solwera rk4
%       b       -   koniec obserwowanego przedziału <0, b>  
%
%   PRZYKŁADOWE WYWOŁANIE
%       >> comparison(0.001,-0.02,0.01,20)
%
[ttr, X1r, X2r, Tr, err, hh] = rk4z(xs1,xs2,h0,b);
[tto, X1o, X2o, To] = ode(xs1,xs2,b);
    plot(Tr,X1r,To,X1o);
    title('Porównanie Trajektori x1(t)');
    xlabel('Czas t');
    ylabel('Rozwiązanie x1');
    legend('RK4', 'ODE45');
    w=waitforbuttonpress;
    clf;
    plot(Tr,X2r,To,X2o);
    title('Porównanie Trajektorii x2(t)');
    xlabel('Czas t');
    ylabel('Rozwiązanie x2');
    legend('RK4', 'ODE45');
    w=waitforbuttonpress;
    clf;
    plot(X1r,X2r,X1o,X2o);
    title('Porównanie Trajektorii Przestrzeń Fazowa (x1,x2)');
    xlabel('Rozwiązanie x1');
    ylabel('Rozwiązanie x2');
    legend('RK4', 'ODE45');
    w=waitforbuttonpress;
    clf;
    barX = categorical({'RK4','ODE45'});
    barX = reordercats(barX,{'RK4','ODE45'});
    bar(barX, [ttr, tto]);
    title('Porównanie czasów wykonania programów [s]');
    w=waitforbuttonpress;
    clf;
    plot(Tr, hh);
    title('Zależność długości kroku od czasu dla metody RK4');
    xlabel('Czas t [s]');
    ylabel('h [s]');
    w=waitforbuttonpress;
    clf;
    plot(Tr, err);
    title('Zależność estymaty błędu od czasu dla metody RK4');
    xlabel('Czas t [s]');
    ylabel('Epsilon');
    w=waitforbuttonpress;
    clearAllMemoizedCaches;
    clf;
    clear all;
    close all;
end