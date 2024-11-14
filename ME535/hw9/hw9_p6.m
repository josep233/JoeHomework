clear
clc
close all

E = 210e9;
L = 1;
A = 0.0004;
rho = 7800;
N = 5;
k1 = 2 * N * E * A / L;
k = N * E * A / L;
m = rho * A * L / N;

syms x

phia = [(x/L)^(1),(x/L)^(2),(x/L)^(3)];
alpha = [((2*1-1)/2)*pi,((2*2-1)/2)*pi,((2*3-1)/2)*pi];
phib = [sin(alpha(1)*x/L),sin(alpha(2)*x/L),sin(alpha(3)*x/L)];

for i = 1:3
    for j = 1:3
        Ma(i,j) = int(phia(i)*phia(j) * rho * A, x, 0, L);
        Ka(i,j) = int(diff(phia(i),x,1)*diff(phia(j),x,1) * E * A, x, 0, L);
        Mb(i,j) = int(phib(i)*phib(j) * rho * A, x, 0, L);
        Kb(i,j) = int(diff(phib(i),x,1)*diff(phib(j),x,1) * E * A, x, 0, L);
    end
end

M = eye(5) * m;
K = [k1 + k, -k, 0, 0, 0; -k, 2*k, -k, 0, 0; 0, -k, 2*k, -k, 0; 0, 0, -k, 2*k, -k; 0, 0, 0, -k, k];
[vectors, values] = eig(K,M);
natural_frequencies = [sqrt(values(1,1)); sqrt(values(2,2)); sqrt(values(3,3)); sqrt(values(4,4)); sqrt(values(5,5))];

[vectorsa, valuesa] = eig(double(Ka),double(Ma));
natural_frequenciesa = [sqrt(valuesa(1,1)); sqrt(valuesa(2,2)); sqrt(valuesa(3,3))];
[vectorsb, valuesb] = eig(double(Kb),double(Mb));
natural_frequenciesb = [sqrt(valuesb(1,1)); sqrt(valuesb(2,2)); sqrt(valuesb(3,3))];