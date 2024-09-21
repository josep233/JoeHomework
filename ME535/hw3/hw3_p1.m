clear
clc
close all

peak_1 = 16e-3;
peak_2 = 12e-3;
damped_period = 0.08;
q_0 = -10e-3;

%part a
log_decrement = log(peak_1 / peak_2);
damping_ratio = log_decrement / (4 * pi^2 + log_decrement^2)^(1/2);
damped_frequency = 2 * pi / damped_period;
damped_frequency_hz = (2 * pi / damped_frequency)^-1;
natural_frequency = damped_frequency / sqrt(1 - damping_ratio^2);
natural_frequency_hz = (2 * pi / natural_frequency)^-1;

%part b
number_to_b = ceil((1 / log_decrement) * log(peak_1 / 0.01e-3));
time_to_b = number_to_b * damped_period;

%part c
%doubling k and halving m would double the natural and damped frequencies,
%decreasing the period and the total time in part b.

%part d
