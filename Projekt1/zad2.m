function [] = zad2()

t = [5,10,25,50,100,200];
eps_ = zeros(3,size(t,2));

for(i=1:size(t,2))
    [A,b] = make_array1(t(1,i));
    eps1(1,i)=t(1,i);
    [~,eps_,k]=GS_solver(A,b,rand(t(1,i),1),1e-8);
    eps1(2,i)=eps_;
    eps1(3,i)=k;
end

time_calc = zeros(3,6);
for(i=1:size(t,2))
    [A,b] = make_array1(t(1,i));
    time_calc(1,i)=t(1,i);
    f = @() GS_solver(A,b,rand(t(1,i),1),1e-8);
    time_calc(2,i)=timeit(f);
    f = @() LU_solver(A,b);
    time_calc(3,i)=timeit(f);
end

plot(eps1(1,:), eps1(2,:));
ylabel("Błąd [eps]");
xlabel("Wymiar Macierzy 1 [n]");
title("Zależność błędu od wymiaru macierzy");
legend("Macierz 1");

w=waitforbuttonpress;

plot(eps1(1,:), eps1(3,:));
ylabel("Liczba iteracji");
xlabel("Wymiar Macierzy 1 [n]");
title("Zależność liczby iteracji od wymiaru macierzy");
legend("Macierz 1");

w=waitforbuttonpress;

plot(time_calc(1,:),time_calc(2,:),time_calc(1,:),time_calc(3,:));
ylabel("Czas obliczenia Ax=b [s]");
xlabel("Wymiar Macierzy 1 [n]");
title("Zależność czasu rozwiązania od wymiaru macierzy");
legend("GS\_solver","LU\_solver");

w=waitforbuttonpress;
close all;
return;

end