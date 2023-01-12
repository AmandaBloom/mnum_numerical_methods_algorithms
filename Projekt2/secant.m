function [xf, ff, iexe, texe] = secant(f, a, b, delta, imax)
%
%   AUTOR
%       Tomasz Pawlak, 304104
%
%   CEL
%       Poszukiwanie pierwiastka funkcji jednej zmiennej
%       metodą siecznych 
%  
%   PARAMETRY WEJSCIOWE
%       f      -  funkcja dana jako wyrazenie  
%       a,b    -  przedział poczatkowy izolacji
%       delta  -  dokladnosc
%       imax   -  maksymalna liczba iteracji
%
%   PARAMETRY WYJSCIOWE
%       xf     -  rozwiazanie 
%       ff     -  wartosc funkcji w xf
%       iexe   -  liczba iteracji wykonanych
%       texe   -  czas obliczen [s]
%
%   PRZYKLADOWE WYWOLANIE
%       >> [xf, ff, iexe, texe] = secant(@ (x) sin(x), 3, 3.5, 1e-8, 100)
%
tic; 
i = 0; xf=a;
while abs(f(xf)) > delta && i < imax
    xf = (a*f(b)-b*f(a))/(f(b)-f(a));
    a = b;
    b = xf;
    i = i + 1;
end
texe=toc; iexe=i;
ff=f(xf);
