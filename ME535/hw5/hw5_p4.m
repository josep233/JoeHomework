clear
clc
close all

d = 1e-2;
A = pi * (d / 2)^2;
L = 1/3;
rho = 7850;
m = rho * A * L;
g = 9.8;
E = 210e9;
k = 2 * E * A / L;
ke = k / 2;
h = 0.5;
delta_s = m * g / k;
natural_frequency = sqrt(k/m);
period = 2 * pi / natural_frequency;

M = eye(3) .* m;
K = [3*k, -k, 0; -k, 2*k, -k; 0, -k, k];
F0 = [m * g; m * g; m * g];

% Initial conditions
initial_displacements = [0; 0; 0]; % Initial displacements
initial_velocities = [-sqrt(2*g*h);-sqrt(2*g*h);-sqrt(2*g*h)]; % Initial velocities
y0 = [initial_displacements; initial_velocities];

% Time span
tspan = [0, 10*period]; % From 0 to 10 seconds

% Solve ODE
[t, y] = ode45(@(t, y) dynamics(t, y, M, K, F0), tspan, y0);

disp1 = y(:,1);
disp2 = y(:,2) - y(:,1);
disp3 = y(:,3) - y(:,2);

hold on
plot(t(:,1),disp1)
plot(t(:,1),disp2)
plot(t(:,1),disp3)
legend()

maxdisp1 = min((disp1(t < 0.35e-3)));
maxdisp2 = max(abs(disp2(t < 0.35e-3)));
maxdisp3 = max(abs(disp3(t < 0.35e-3)));

maxstress1 = 2*k * maxdisp1 / A;
maxstress2 = k * maxdisp2 / A;
maxstress3 = k * maxdisp3 / A;

amplitude_c = abs((L * m * g / (2 * E * A)) - (sqrt(2 * 9.8 * 0.5)/natural_frequency)*1i - (m * natural_frequency^2)^(-1));
max_deflection_c = amplitude_c + (m * natural_frequency^2)^(-1);
max_stress_c = k * max_deflection_c / A;

function dydt = dynamics(t, y, M, K, f)
    n = size(M, 1);
    q = y(1:n);       % Displacement
    q_dot = y(n+1:end); % Velocity
    acc = M \ (f - K * q); % Acceleration
    
    dydt = [q_dot; acc]; % Combine velocity and acceleration
end

