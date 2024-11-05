clear
clc
close all

m = 1;
kx = (200 * pi)^2;
ky = kx;
ktheta = 0.1 * kx;
theta_deg = 30;
theta_rad = deg2rad(theta_deg);
beta = 3e-5;

M = [m, 0; 0, m];
K = [kx + ktheta * cos(theta_rad)^2, ktheta * cos(theta_rad) * sin(theta_rad); ktheta * cos(theta_rad) * sin(theta_rad), ky + ktheta * sin(theta_rad)^2];
C = 3e-5 * K;

[vectors, values] = eig(K,M);

natural_frequencies = [sqrt(values(1,1)),sqrt(values(2,2))];
vector1 = vectors(:,1) / vectors(1,1);
vector2 = vectors(:,2) / vectors(1,2);

Phi1 = vector1 ./ (sqrt(transpose(vector1) * M * vector1));
Phi2 = vector2 ./ (sqrt(transpose(vector2) * M * vector2));

zeta1 = (beta / 2) * natural_frequencies(1);
zeta2 = (beta / 2) * natural_frequencies(2);