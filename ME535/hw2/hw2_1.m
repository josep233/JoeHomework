clear
clc
close all

%time vector
t = linspace(0,0.5,10000);

%complex amplitude, calculated on paper
A_hat = 0.01 / 1i - 0.02 * exp(-1i * 0.3 * pi);

%function parameters
omega = 50;
tau = 2 * pi / omega;

%answer to part b
tau / 2

%function plotting
q = real(A_hat * exp(1i * omega * t));

%calculate earliest theta when q = 0
theta = atan(0.0062 / -0.0118);

%calculate the first time when q = 0
first_time = (pi / 2 - theta) / 50;

%calculate amplitude magnitude
A_mag = sqrt((-0.0118)^2 + (0.0062)^2);

hold on
plot(t,q)
yline(0)