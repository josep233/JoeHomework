clear
clc
close all

m1 = 100;
m2 = 200;
m3 = 300;
omega1 = 40;
omega2 = 50;
omega3 = 60;
k1 = omega1^2 * m1;
k2 = omega2^2 * m2;
k3 = omega3^2 * m3;

M = [m1, 0, 0; 0, m2, 0; 0, 0, m3];
K = [k1, -k1, 0; -k1, k2 + k1, -k2; 0, -k2, k3 + k2];

[vectors, values] = eig(K,M);

natural_frequencies = [sqrt(values(1,1)), sqrt(values(2,2)), sqrt(values(3,3))];
eigenvectors1 = vectors(:,1) / vectors(1,1);
eigenvectors2 = vectors(:,2) / vectors(1,2);
eigenvectors3 = vectors(:,3) / vectors(1,3);

