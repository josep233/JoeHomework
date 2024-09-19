clear
clc
close all

F0 = 1;
beta = 0.1;
T = 5;
t = linspace(0,10,100);
exponential = exp(-beta * (t - 2*T)) - F0;
combined = t;
combined(t < T) = F0;
combined(t > T) = F0 * exp(-beta * (t(t > T) - T));

hold on
plot(t,combined,'black')
plot(t,exponential,'red')