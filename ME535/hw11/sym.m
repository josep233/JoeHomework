clear
clc
close all

dCL_dalpha = 2 * pi;
x_alpha = 0.05;
r_alpha = 0.5;

syms m b K_alpha K_h e rho S omega_h omega_alpha v

I_alpha = m * b^2 * r_alpha^2;
S_alpha = m * b * x_alpha;
q = (1/2) * rho * v^2;

eq1 = omega_h / omega_alpha == 0.5;
eq2 = e / b == 0.4;
eq3 = 2 * m / (pi * rho * b * S) == 7.5;
eq4 = omega_h^2 == K_h/m;
eq5 = omega_alpha^2 == K_alpha/I_alpha;

M11 = m * b;
M12 = S_alpha;
M21 = S_alpha;
M22 = I_alpha;

K11 = K_h * b;
K12 = q * S * dCL_dalpha;
K21 = 0;
K22 = K_alpha - q * S * e * dCL_dalpha;

M = [M11,M12;M21,M22];
K = [K11,K12;K21,K22];


