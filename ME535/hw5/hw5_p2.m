clear
clc
close all

d = 2e-2;
L = 0.5;
rho = 1000;
E = 9e9;
failure_strain = 0.01/100;
gamma = 0.04;
m = pi * (d/2)^2 * L * rho;
I = (pi/64)*d^4;

xmax = L^3 / (3 * E * I * gamma);
stressmax = 2 * E * (d/2) / (2 * L^2);

k = 3 * E * I / L^3;
resonant_frequency = sqrt(k/m);
resonant_frequency_hz = (resonant_frequency / (2 * pi));

p_actual = 55;

p0 = ((0.01 / 100) * 2 * L^2 / (3 * d)) * 3 * E * I * gamma / (.001 * L^3);