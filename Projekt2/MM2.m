function [r] = MM2(p, x_k, tol)
%
%   AUTOR
%       Tomasz Pawlak, 304104
%
%   CEL
%       Poszukiwanie jednego pierwiastka wielomianu
%       Drugą Metodą Muellera
%  
%   PARAMETRY WEJSCIOWE
%       p      -  Wielomian dany jako tabela
%       x_k    -  Przybliżenie zera
%       tol    -  Tolerancja docelowa wyniku
%
%   PARAMETRY WYJSCIOWE
%       r      -  Pierwiastek rozwiązania z dokładnością +- tol
%
%   PRZYKLADOWE WYWOLANIE
%       >> r = MM2([2 3 -6 4 7], -2.5, 1e-8)
%

%while (abs(f(x_k)) > 1e-8)
while abs(polyval(p,x_k)) > tol
    %sqrt_ = sqrt(df(x_k)^2-2*f(x_k)*d2f(x_k));
    sqrt_ = sqrt(polyval(polyder(p),x_k)^2-2*polyval(p,x_k)*polyval(polyder(polyder(p)),x_k));
    
    %z1 = -2*f(x_k)/(df(x_k)+sqrt_);
    z1 = -2*polyval(p,x_k)/(polyval(polyder(p),x_k)+sqrt_);
    
    %z2 = -2*f(x_k)/(df(x_k)-sqrt_);
    z2 = -2*polyval(p,x_k)/(polyval(polyder(p),x_k)-sqrt_);
    
    %Wybór elementu paraboli o mniejszym module
    %if abs(df(x_k)+sqrt_) > abs(df(x_k)-sqrt_)
    if abs(polyval(polyder(p),x_k)+sqrt_) > abs(polyval(polyder(p),x_k)-sqrt_)
        z_min = z1;
    else
        z_min = z2;
    end
    r = x_k+z_min;
    x_k = r;
end
end
