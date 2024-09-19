clear
clc
close all

L = 2e-2;
m = 3e-3;
kappa = 20;
damping_ratio = 0.02;
angle_deg = 30;
F0 = kappa * deg2rad(angle_deg) / L;
q_0 = deg2rad(angle_deg);
q_dot_0 = 0;
natural_frequency = sqrt(kappa / (1/3 * m * L^2));
damped_frequency = natural_frequency * sqrt(1 - damping_ratio^2);
damped_period = 2 * pi / damped_frequency;
T = 2.5 * damped_period;
c = damping_ratio * 2 * sqrt(kappa*m);

ts = linspace(0,10*T,100);

[tout,yout] = ode45(@(t, theta) eom_2_12(t, theta, L, m, kappa, damping_ratio, F0, natural_frequency, damped_frequency, T, c), [ts(1),ts(end)],[q_0; q_dot_0]);

plot(tout,yout(:,1))
ylabel('${\theta}$ (rad)','interpreter','latex')
xlabel('time (s)','interpreter','latex')

function [theta_dot] = eom_2_12(t, theta, L, m, kappa, damping_ratio, F0, natural_frequency, damped_frequency, T, c)    
    % Forcing - sum of step and ramp
    F = F0 * (t < T && t > 0);
    % F = 0;
    
    % Equations of Motion
    theta_dot(1,1) = theta(2);
    theta_dot(2,1) = (3/(m*L^2)) * (F - c * theta(2) - kappa * theta(1));
end