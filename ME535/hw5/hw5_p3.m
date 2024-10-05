clear
clc
close all

d = 1e-2;
A = pi * (d / 2)^2;
L = 1;
rho = 7850;
m = rho * A * L;
g = 9.8;
E = 210e9;
k = 2 * E * A / L;
delta_s = m * g / k;
natural_frequency = sqrt(k/m);


sigma_a = rho * L * g;

amplitude_b = (L * m * g / (2 * E * A)) - (m * natural_frequency^2)^(-1);
max_deflection_b = amplitude_b + (m * natural_frequency^2)^(-1);

amplitude_c = abs((L * m * g / (2 * E * A)) - (sqrt(2 * 9.8 * 0.5)/natural_frequency)*1i - (m * natural_frequency^2)^(-1));
max_deflection_c = amplitude_c + (m * natural_frequency^2)^(-1);
stress_c = k * max_deflection_c / A;