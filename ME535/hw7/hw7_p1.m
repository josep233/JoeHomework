clear
clc
close all

syms lambda

poly = (200 - 4*lambda) * (800 - 2*lambda) - 200^2;
r = root(poly,lambda);

sqrt(double(r(1)))
sqrt(double(r(2)))

M = [4, 0; 0, 2];
K = [200, 200; 200, 800];
[vectors, values] = eig(K,M);

natural_frequencies = [sqrt(values(1,1)), sqrt(values(2,2))];
eigenvectors1 = vectors(:,1) / vectors(1,1);
eigenvectors2 = vectors(:,2) / vectors(1,2);
