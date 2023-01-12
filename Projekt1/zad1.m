function [] = zad1()

t = [5,10,25,50,100,200];
epsA1 = zeros(3,size(t,2));
epsA2 = epsA1;

for(i=1:size(t,2))
    epsA1(1,i)=t(1,i);
    [A,b] = make_array1(t(1,i));
    
    X = LU_solver(A,b);
    eps_i = norm(A*X-b,2);
    epsA1(2,i)=eps_i;
    epsA1(3,i)=eps_i/eps;
end

for(i=1:size(t,2))
    epsA2(1,i)=t(1,i);
    [A, b] = make_array2(t(1,i));

    X = LU_solver(A,b);
    eps_i = norm(A*X-b,2);
    epsA2(2,i)=eps_i;
    epsA2(3,i)=eps_i/eps;
end
plot(epsA1(1,:),epsA1(2,:),epsA2(1,:),epsA2(2,:));
ylabel("Błąd");
xlabel("Wymiar macierzy [n]");
title("Zależność błędu od wymiaru macierzy");
legend("Macierz 1", "Macierz 2");

w=waitforbuttonpress;

plot(epsA1(1,:),epsA1(3,:),epsA2(1,:),epsA2(3,:));
ylabel("Błąd [eps]");
xlabel("Wymiar macierzy [n]");
title("Zależność stosunku błąd/eps od wymiaru macierzy");
legend("Macierz 1", "Macierz 2");

w=waitforbuttonpress;
close all;
return

end