clear
clc
close all

d = 1e-2;
A = pi * (d / 2)^2;
L = 1;
rho = 7850;
m = rho * A * L;
g = 9.8;
E = 210e9;
k = 2 * E * A / L;
h = 0.5;
delta_s = m * g / k;
natural_frequency = sqrt(k/m);
period = 2 * pi / natural_frequency;

t = linspace(0,3*period,100);

z_0 = 0;
z_dot_0 = sqrt(2*g*h);

[tout,yout] = ode45(@(t, z) eom(z, m, k), [t(1),t(end)],[z_0; z_dot_0]);

plot(tout,yout(:,1))

amplitude_c = abs((L * m * g / (2 * E * A)) - (sqrt(2 * 9.8 * 0.5)/natural_frequency)*1i - (m * natural_frequency^2)^(-1));
max_deflection_c = amplitude_c + (m * natural_frequency^2)^(-1);
max_stress_c = k * max_deflection_c / A;

amplitude = max(yout(:,1));

function [z_dot] = eom(z, m, k)    
    % Forcing - sum of step and ramp
    % F = 0;
    F = m * 9.8;
    
    % Equations of Motion
    z_dot(1,1) = z(2);
    z_dot(2,1) = (1/m) * (F - k * z(1));
end