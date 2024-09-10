clear
clc
close all

%time vector
t = linspace(0,0.5,10000);

q_real = 0.01 * sin(50 * t) - 0.02 * cos(50 * t - 0.3 * pi);

%complex amplitude, calculated on paper
A_hat = 0.01 / 1i - 0.02 * exp(-1i * 0.3 * pi);

%function parameters
omega = 50;
tau = 2 * pi / omega;
tau / 2

q = real(A_hat * exp(1i * omega * t));

hold on
plot(t,q)
yline(0)