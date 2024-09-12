clear
clc
close all

%known parameters
t = linspace(0,5,100);
time_at_first_zero = 1;
tau = 2 * 3;
amplitude_magnitude = 20;

%calculated with omega = 2 pi / tau
omega = (2 * pi) / tau;

%calculate first angle at t = 0
first_theta = pi / 2 - time_at_first_zero * omega;

%calculate real / imaginary parts of complex amplitude
imaginary_part = sin(first_theta) * amplitude_magnitude;
real_part = cos(first_theta) * amplitude_magnitude;

%calculate complex amplitude
A_hat = real_part + imaginary_part * 1i;

%calculate x. Assume cosine function because function starts high
func = real(A_hat * exp(1i * omega * t));

%calculate x dot with calculus
func_prime = real(A_hat * omega * 1i * exp(1i * omega * t));

%calculate x dot dot with calculus
func_prime_prime = real(A_hat * omega * 1i * omega * 1i * exp(1i * omega * t));

%calculate time at which first minimum value of x occurs
time_first_minimum = time_at_first_zero + tau / 4;

%calculate time at which first maximum value of x_dot occurs
time_first_maximum_xdot = time_at_first_zero + tau / 2;

%calculate time at which first maximum value of x_dot_dot occurs
time_first_maximum_xdot_dot = time_at_first_zero + tau / 4;

hold on
plot(t,func,'blue')
plot(t,func_prime,'red')
yline(0)
xlabel('time (s)')
ylabel('x (mm)')
legend(["${x}$","${\dot{x}}$"],'interpreter',['latex'],'location','northwest','FontSize',12)

%part A
disp("x(0) = "+func(1)+" mm")
disp("x_dot(0) = "+func_prime(1)+" mm/ms")

%part b
disp("time of first minimum x value = "+time_first_minimum+" ms")

%part c
disp("time of first maximum x_dot value = "+time_first_maximum_xdot+" ms")

%part d
disp("time of first maximum x_dot_dot value = "+time_first_maximum_xdot_dot+" ms")