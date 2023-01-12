function [roots, t] = z2(p, x_k, tol)
%
%   AUTOR
%       Tomasz Pawlak, 304104
%
%   CEL
%       Poszukiwanie wszystkich pierwiastka rzeczywistych oraz
%       zespolonych Drugą Metodą Muellera
%  
%   PARAMETRY WEJSCIOWE
%       p      -  Wielomian dany jako tabela
%       x_k    -  Przybliżenie wybranego miejsca zerowego
%       tol    -  Tolerancja docelowa wyniku
%
%   PARAMETRY WYJSCIOWE
%       roots  -  Tablica pierwiastków wielomianu
%
%   PRZYKLADOWE WYWOLANIE
%       >> roots = z2([2 3 -6 4 7], -2.5, 1e-8)
%

poly_ = p;
tic;
% Stopień wielomianu jednoznacznie definiuje liczbę pierwiastków
deg=length(p)-1;
% Inicjacja tablicy pierwiastków
roots=zeros(deg,1);
for i=1:deg
    %Wykorzystanie metody MM2
    roots(i)=MM2(p, x_k, tol);
    
    n=length(p);
    new_p=zeros(n,1);
    % Deflacja czynnikiem liniowym - schemat Hornera
    for j = 2:n
        new_p(j)=p(j-1)+(new_p(j-1)*roots(i));
    end
    % Wielomian zostaje uproszczony po podzieleniu go przez (x-roots(i))
    p = new_p;
end
t = toc;
% Porządkowanie pierwiastków przed returnem funkcji
roots = sort(roots);

% Rysowanie wykresu
% Inicjacja tabeli pierwiastków rzeczywistych
plot_roots = [];
for i=1:length(roots)
    if abs(roots(i)) == abs(real(roots(i))) | abs(imag(roots(i))) < tol
        %Dopisywanie pierwiastków rzeczywistych do tablicy
        plot_roots = [plot_roots, real(roots(i))];
    end
end
% Najmniejszy i największy pierwiastek rzeczywisty
a = min(plot_roots);
b = max(plot_roots);
% Określenie zakresu próbek do wyświetlenia
if a==b
    samples = a-5:0.2:a+5;
else
    samples = a-0.1*(b-a):0.01*(b-a):b+0.1*(b-a);
end
if isempty(plot_roots)
    samples = -30:1:30;
end

% Rysowanie wykresu
plot(samples, polyval(poly_, samples), plot_roots, polyval(poly_, plot_roots), '.', ...
    "MarkerSize", 12, "Color", "black");
title_1 = num2str(reshape(poly_, 1, []));
title_2 = "Pierwiastki wielomianu [";
title_3 = "]";
title(strcat(title_2, title_1, title_3));

% Dodatkowe rysowanie osi rzędnych i odciętych
ax = gca;
ax.XAxisLocation = 'origin';
if a<0 & b>0
    ax.YAxisLocation = 'origin';
end

end