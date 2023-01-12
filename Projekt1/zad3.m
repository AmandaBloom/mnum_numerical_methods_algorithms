function [] = zad3()
x_data=(-10:2:10)'; % Współrzędna x macierzy danych
y_data=[-3.578,-5.438,-4.705,-3.908,-2.069,0.942,...
    -0.725,-4.128,-11.160,-23.440,-42.417]'; % Współrzędna y mac. danych
samples=-10:0.2:10; % Wektor 1 wart. próbek danych dla aproks. funkcji

poly_grades=[3,5,7,9,10]'; % Wielomiany stopni 3,5,7,9,10
y_normal=zeros(5,size(samples,2)); % Inicjacja wektorów y, gdzie y=f(x)
y_svd=zeros(5,size(samples,2)); % celem próbkowania
eps_normal = zeros(2,5); % Inicjacja macierzy błędów aproksymacji dla 
eps_svd = zeros(2,5);    % obu metod

% Układ r-ń normalnych
plot(x_data,y_data,'.',"MarkerSize",12,"Color","black");
set(gcf, 'WindowState', 'maximized');
ylabel("Y");
xlabel("X");
title("Układ równań normalnych");
hold on;

for i=1:5 
    % Generacja macierzy współczynników dla danego st. wielomianu
    A=coef_matrix(poly_grades(i,1), x_data); 
    % Obliczanie wekt. współ. wielomianu
    x_normal=linsolve(A' * A,A' * y_data);
    % Obliczanie wart. wielomianu dla próbek
    y_normal(i,:)=polyval(x_normal',samples);
    % Nanoszenie na wykres
    plot(samples,y_normal(i,:));  
    % Obliczanie wartości błędów
    eps_normal(1,i)=norm(y_data - A * x_normal,2);
    eps_normal(2,i)=norm(y_data - A * x_normal,inf);
end
legend("y=f(x)","wielomian 3 stopnia","wielomian 5 stopnia",...
    "wielomian 7 stopnia","wielomian 9 stopnia","wielomian 10 stopnia",...
    "Location","southwest");
hold off;
w=waitforbuttonpress;
clf;

% Rozkład SVD 
plot(x_data,y_data,'.',"MarkerSize",12,"Color","black");
ylabel("Y");
xlabel("X");
title("Rozkład SVD");
hold on;

for i=1:5
    % Generacja macierzy współczynników dla danego st. wielomianu
    A=coef_matrix(poly_grades(i,1), x_data);
    % Wzory 3.39-40 (macierz pseudoodwrotna)
    [U,E,V]=svd(A);
    E=E(1:size(A,2),:); % Skracanie rozmiaru do nxn, aby wykonać inwersję
    Eplus=[inv(E) zeros(size(A,2),size(A,1)-size(A,2))];
    pinvA=V*Eplus*U.';
    %Obliczanie wektora współcz. wielomianu
    %Takie same wyniki dla wbudowanej f-kcji x_svd=pinv(A)*y_data;
    x_svd=pinvA*y_data;
    % Obliczanie wart. wielomianu dla próbek
    y_svd(i,:)=polyval(x_svd',samples);
    % Nanoszenie wartości próbek na wykres
    plot(samples,y_svd(i,:));
    % Obliczanie wartości błędów
    eps_svd(1,i)=norm(y_data - A * x_svd,2);
    eps_svd(2,i)=norm(y_data - A * x_svd,inf);
end
legend("y=f(x)","wielomian 3 stopnia","wielomian 5 stopnia",...
    "wielomian 7 stopnia","wielomian 9 stopnia","wielomian 10 stopnia",...
    "Location","southwest");
hold off;
w=waitforbuttonpress;
clf;


% Bład aproksymacji euklides
plot(poly_grades,eps_normal(1,:),poly_grades,eps_svd(1,:));
ylabel("Błąd");
xlabel("Stopień Wielomianu");
title("Zależność normy euklidesowej aproksymacji od stopnia wielomianu");
legend("Układ r-ń normalnych","Rozkład SVD");
w=waitforbuttonpress;
clf;

% Błąd aproksymacji inf
plot(poly_grades,eps_normal(2,:),poly_grades,eps_svd(2,:));
ylabel("Błąd");
xlabel("Stopień Wielomianu");
title("Zależność normy maksimum aproksymacji od stopnia wielomianu");
legend("Układ r-ń normalnych","Rozkład SVD");
w=waitforbuttonpress;
clf;

% Porównanie błędów
loglog(poly_grades,abs(eps_normal(1,:)-eps_svd(1,:)),poly_grades,...
    abs(eps_normal(2,:)-eps_svd(2,:)),"Marker",".","MarkerSize",20 );
ylabel("Różnica błędów");
xlabel("Stopień Wielomianu");
title("Zależność bezwzg. różnicy błędów od stopnia wielomianu");
legend("Norma 2","Norma Inf");
w=waitforbuttonpress;
clf;

close all;
end
