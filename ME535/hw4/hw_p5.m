clear
clc
close all

m1 = 0.5;
m2 = 1;
k1 = 3.2e3;
k2 = 0;
mu1 = 0;
mu2 = 40;
A = 20e-3;
omega = 75;

F = (-A * 1i) * (k2 + 1i * omega * mu2 - (m1+m2) * omega^2);

abs(F)
atan2(imag(F)/real(F),1)