clear
clc
close all

F0 = 1;
beta = 0.1;
T = 5;
t = linspace(0,10,1000);
exponential = F0 .* heaviside(t) - F0 .* heaviside(t - T) + F0 .* exp(-beta .* (t - T)) .* heaviside(t - T);
combined = t;
combined(t < T) = F0;
combined(t > T) = F0 * exp(-beta * (t(t > T) - T));
step_solution = @(t) (1./(m.*natural_frequency.^2)) .* (1 - exp(-damping_ratio .* natural_frequency .* t) .* cos(damped_frequency .* t) + (damping_ratio .* natural_frequency ./ damped_frequency) .* sin(damped_frequency .* t));

hold on
plot(t,combined,'black')
plot(t,exponential,'red')