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

S = chord * span;
slope = 2 * pi;

v = sqrt((k/(S * slope) / (0.5 * rho)));

