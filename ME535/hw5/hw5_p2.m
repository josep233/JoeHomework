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

k = 3 * E * I / L^3;
resonant_frequency = sqrt(k/m);
resonant_frequency_hz = (resonant_frequency / (2 * pi));

p0 = 0.03 * E * I / (.001 * L^2);