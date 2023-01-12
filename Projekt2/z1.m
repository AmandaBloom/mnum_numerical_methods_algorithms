function [tabS, tabN, imaxs, imaxn] = z1(f,s,a,b,x0,delta,imax)
%
%   AUTOR
%       Tomasz Pawlak, 304104
%
%   CEL
%       Porównanie kolejnych iteracji znajdywania wybranego miejsca
%       zerowych metodami Newtona oraz Siecznych oraz gene-
%       racja tablic z tymi wynikami i wykresów
%
%   PARAMETRY WEJŚCIOWE
%       f       - funkcja dana jako wyrażenie
%       s       - przedział obserwacji dany jako wektor 2-arg
%       a,b     - przedział początkowy obserwacji dla m.siecznych
%       x0      - punkt początkowy dla m. Newtona
%       delta   - dokładność rozwiązania
%       imax    - maksymalna l-a iteracji dla obu metod
%
%   PARAMETRY WYJŚCIA
%       tabS, tabN      - tablice kolejnych iteracji przybliżeń
%                       punktów zerowych
%       imaxs, imaxn    - liczba wykonanych iteracji obu metod
%
%   PRZYKŁADOWE WOŁANIE
%       [tabS, tabN, imaxs, imaxn] = z1(@(x)1.2*sin(x)+2*log(x+2)-5,[2 12],3,11,9,1e-8,15)
%
samples=(s(1):(s(2)-s(1))/100:s(2))';
Y = f(samples);

plot(samples, Y);
ylabel("Y");
xlabel("X");
title("Metoda siecznych");
hold on;

Xn=zeros(imax,1);
Xs=Xn;
Yn=zeros(imax,1);
Ys=Yn;

i = 0; xf=a;
while abs(f(xf)) > delta && i < imax
    xf = (a*f(b)-b*f(a))/(f(b)-f(a));
    a = b;
    b = xf;
    i = i + 1;
    Xs(i) = xf;
    Ys(i) = f(xf);
end
imaxs = i;
fprintf("Secant (x,y)=%d,%d\n",Xs(imaxs),Ys(imaxs));

plot(Xs(1:imaxs,1), Ys(1:imaxs,1),'.',...
    "MarkerSize",12,"Color","black");
legend("y=f(x)", "Iteracje metody siecznych",...
    "Location", "southwest");
ax = gca;
ax.XAxisLocation = 'origin';
hold off;
w=waitforbuttonpress;
clf;

plot(samples, Y);
ylabel("Y");
xlabel("X");
title("Metoda newtona");
hold on;

syms X
df = matlabFunction(diff(f(X), X));

i = 0; x=x0; fx=feval(f,x); 
while abs(fx) > delta && i < imax
     i = i + 1;
     x = x - fx/df(x);
     fx=feval(f,x);
     Xn(i) = x;
     Yn(i) = f(x);
end
imaxn = i;
fprintf("Newton (x,y)=%d,%d\n",Xn(imaxn),Yn(imaxn));

plot(Xn(1:imaxn,1), Yn(1:imaxn,1),'.',...
    "MarkerSize",12,"Color","black");
legend("y=f(x)", "Iteracje metody newtona",...
    "Location", "southwest");
ax = gca;
ax.XAxisLocation = 'origin';
hold off;
w=waitforbuttonpress;
clf;
close all;

tabS = [Xs(1:imaxs,1) Ys(1:imaxs,1)];
tabN = [Xn(1:imaxn,1) Yn(1:imaxn,1)];
