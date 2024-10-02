clear
clc
close all

m1 = 0.5;
m2 = 1;
k = 3.2e3;
mu = 40;
A = 20e-3;
omega = 75;

a = k - m1 * omega^2;
b = 1i*omega*mu - m2*omega^2;
F = -A*1i/(a/b);

abs(F)