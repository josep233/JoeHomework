clear
clc
close all

m = 1;
g = 9.8;
L = 1;
k = 0.05 * m * g / L;

M = [m, 0, 0; 0, m, 0; 0, 0, m];
K = [k, -k, 0; -k, k, -k; 0, -k, k];

[vectors, values] = eig(K,M);

natural_frequencies = [sqrt(values(1,1)),sqrt(values(2,2)),sqrt(values(3,3))];
vector1 = vectors(:,1) / vectors(1,1);
vector2 = vectors(:,2) / vectors(1,2);
vector3 = vectors(:,3) / vectors(1,3);

Phi1 = vector1 ./ (sqrt(transpose(vector1) * M * vector1));
Phi2 = vector2 ./ (sqrt(transpose(vector2) * M * vector2));
Phi3 = vector3 ./ (sqrt(transpose(vector3) * M * vector3));
