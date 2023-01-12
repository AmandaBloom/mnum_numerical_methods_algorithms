function [] = plot_function(f,a,b)
%
%   AUTOR
%       Tomasz Pawlak, 304104
%
%   CEL
%       Narysowanie przebiegu funkcji w danym przedziale
%
%   PARAMETRY WEJŚCIOWE
%       f       - funkcja dana jako wyrażenie
%       a,b     - przedział próbkowania funkcji
%
%   PRZYKŁADOWE WYWOŁANIE
%   >> plot_function(@(x) sin(x), 0, 6.28)
%   >> plot_function(@(x) 2*x.^(4)+3*x.^(3)-6*x.^(2)+4*x+7, -4, 3)

samples=(a:(b-a)/100:b)';
y = f(samples);

plot(samples, y);
ylabel("Y");
xlabel("X");
title("Metoda siecznych");
ax = gca;
ax.XAxisLocation = 'origin';