clear
clc
close all

rho = 1.2;
m = 6.47;
chord = 0.36;
span = 1;
I = (1/12) * m * chord^2;
k = 97.7;
natural_frequency = 5.95 * 2 * pi;
K_alpha = natural_frequency^2 * I;
dCL_da = 0.071 * 180 / pi; %from paper?

S = chord * span;
slope = 2 * pi;

e = chord * (0.374 - 0.227);

qD = K_alpha / (S * e * dCL_da);

speed = sqrt(qD * 2 / rho);

