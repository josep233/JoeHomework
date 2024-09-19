clear
clc
close all

%known parameters
omega_n = 1000;
zeta = 0.015;
omega_d = omega_n * sqrt(1 - zeta^2);
M = 1;
T = 0.01;

%step and impulse functions from book
x1 = @(t) (1./(M.*omega_n)) .* exp(-zeta .* omega_n .* t) .* sin(omega_d .* t) .* heaviside(t);
x2 = @(t) (1./(M.*omega_n.^2)) .* (1 - exp(-zeta .* omega_n .* (t-T)) .* (cos(omega_d .* (t-T)) + (omega_n .* zeta ./ omega_d) .* sin(omega_d .* (t-T)))) .* heaviside((t-T));
x =  @(t) x1(t) + 1000 .* x2(t);

%define time
t = linspace(0,0.07,100);

%plotting
hold on
plot(t,x(t))
plot(t,x1(t),'g')
plot(t,1000*x2(t),'red')
yline(0,'black')

%part b
x_guess = 0.022;
time_actual = fzero(x,x_guess);