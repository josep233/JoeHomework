clear
clc
close all

%unit system: mm g s N

L = 2e-2;
m = 3e-3;
kappa = 20;
damping_ratio = 0.02;
angle_deg = 30;
F0 = kappa * deg2rad(angle_deg) / L;
q_0 = deg2rad(angle_deg);
q_dot_0 = 0;
natural_frequency = sqrt(kappa / m);
damped_frequency = natural_frequency * sqrt(1 - damping_ratio^2);
damped_period = 2 * pi / damped_frequency;
T = 3 * damped_period;
c = damping_ratio * 2 * sqrt(kappa*m);

ts = linspace(0,10*T,725);
% ts = linspace(0,0.1,100);

step_solution = @(t) (1./(m.*natural_frequency.^2)) .* (1 - exp(-damping_ratio .* natural_frequency .* t) .* cos(damped_frequency .* t) + (damping_ratio .* natural_frequency ./ damped_frequency) .* sin(damped_frequency .* t));

q_analytic = exp(-damping_ratio .* natural_frequency .* ts) .* (q_0 .* cos(damped_frequency .* ts) + (q_dot_0 + damping_ratio .* natural_frequency .* q_0 ./ damped_frequency) .* sin(damped_frequency .* ts));
q_f_analytic = F0 .* (step_solution(ts) .* heaviside(ts) - step_solution(ts - T) .* heaviside(ts - T));
hold on
plot(ts,q_f_analytic + q_analytic,'black')

[tout,yout] = ode45(@(t, theta) eom_2_12(t, theta, L, m, kappa, damping_ratio, F0, natural_frequency, damped_frequency, T, c), [ts(1),ts(end)],[q_0; q_dot_0]);

plot(tout,yout(:,1),'red')
ylabel('${\theta}$ (rad)','interpreter','latex')
xlabel('time (s)','interpreter','latex')

function [theta_dot] = eom_2_12(t, theta, L, m, kappa, damping_ratio, F0, natural_frequency, damped_frequency, T, c)    
    % Forcing - sum of step and ramp
    F = F0 * (heaviside(t) - heaviside(t - T));
    % F = 0;
    
    % Equations of Motion
    theta_dot(1,1) = theta(2);
    % theta_dot(2,1) = ((1/3) * m * L^2)^(-1) * (F - c * theta(2) - kappa * theta(1));
    theta_dot(2,1) = (1/m) * (F - c * theta(2) - kappa * theta(1));
end