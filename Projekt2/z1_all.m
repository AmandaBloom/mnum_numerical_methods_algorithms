function [tabS, tabN, t] = z1_all(f,s,ab_S,x_N,delta,imax)
%
%   AUTOR
%       Tomasz Pawlak, 304104
%
%   CEL
%       Generacja wszystkich miejsc zerowych funkcji za pomocą algorytmów
%       wykorzystujących metodę siecznych oraz stycznych (Newtona)
%
%   PARAMETRY WEJŚCIOWE
%       f       - funkcja dana jako wyrażenie
%       s       - przedział obserwacji dany jako wektor 2-arg
%       ab_S    - przedziały początkowe dane jako tabela
%       x_N     - punkty początkowe dla metody Newtona dane jako tabela
%       delta   - dokładność rozwiązania
%       imax    - maksymalna l-a iteracji dla obu metod
%
%   PARAMETRY WYJŚCIA
%       tabS, tabN      - macierze rozwiązań punktów zerowych obu metod
%
%   PRZYKŁADOWE WYWOŁANIE
%   >> [tabS, tabN, t] = z1_all(@(x)1.2*sin(x)+2*log(x+2)-5,[2 12],[6 8 8 10],[7 9],1e-8,15)
%   >> [tabS, tabN, t] = z1_all(@(x)2*x.^(4)+3*x.^(3)-6*x.^(2)+4*x+7,[-5 5],[-3 -2 -2 0],[-2.5 -0.5],1e-8,15)
%
samples=(s(1):(s(2)-s(1))/100:s(2))';
Y = f(samples);
xf_s = [];
ff_s = xf_s;
xf_n = [];
ff_n = xf_n;
t = zeros(2,1);

%%% Metoda Siecznych
tic;
for i=1:length(ab_S)/2
    [xf, ff] = secant(f, ab_S(2*i-1), ab_S(2*i), delta, imax);
    if xf>s(1) & xf<s(2)
        xf_s = [xf_s; xf];
        ff_s = [ff_s; ff];
    else
        disp("Nie znalezienio rozwiązania metodą siecznych w przedziale ",...
            ab_S(2*i-1), " ", ab_S(2*i));
    end
end
t(1)=toc;
tabS = [xf_s ff_s];

plot(samples, Y);
ylabel("Y");
xlabel("X");
title("Metoda siecznych");
hold on;
plot(xf_s, ff_s,'.',...
    "MarkerSize",12,"Color","black");
legend("y=f(x)", "Rozwiązania metodą siecznych",...
    "Location", "northeast");
ax = gca;
ax.XAxisLocation = 'origin';
hold off;
w=waitforbuttonpress;
clf;


%%% Metoda Newtona
tic;
for i=1:length(x_N)
    [xf, ff] = newton(f, x_N(i), delta, imax);
    if xf>s(1) & xf<s(2)
        xf_n = [xf_n; xf];
        ff_n = [ff_n; ff];
    else
        disp("Nie znalezienio rozwiązania metodą Newtona w otoczeniu punktu x0=", ...
            x_N(i), " z poprawnego przedziału wejściowego");
    end
end
t(2)=toc;
tabN = [xf_n ff_n];


plot(samples, Y);
ylabel("Y");
xlabel("X");
title("Metoda newtona");
hold on;
plot(xf_n, ff_n,'.',...
    "MarkerSize",12,"Color","black");
legend("y=f(x)", "Rozwiązania metodą Newtona",...
    "Location", "northeast");
ax = gca;
ax.XAxisLocation = 'origin';
hold off;
w=waitforbuttonpress;
clf;
close all;