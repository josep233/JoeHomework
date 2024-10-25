clear
clc
close all

m = 1;
k = 1;

M = [m, 0; 0, m];
K = [2*k, -k; -k, k];

[vectors, values] = eig(K,M);

natural_frequencies = [sqrt(values(1,1)), sqrt(values(2,2))];
periods = 2 .* pi ./ natural_frequencies;
eigenvectors1 = vectors(:,1) / vectors(1,1);
eigenvectors2 = vectors(:,2) / vectors(1,2);

y0_a = [1,0,0,0];
y0_b = [1,1,0,0];
y0_c = [eigenvectors1(1),eigenvectors1(2),0,0];
y0_d = [0,0,eigenvectors2(1),eigenvectors2(2)];

% Time span
tspan = [0, 30]; % From 0 to 10 seconds

% Solve ODE
[ta, ya] = ode45(@(t, y) eom(t, y, M, K), tspan, y0_a);
[tb, yb] = ode45(@(t, y) eom(t, y, M, K), tspan, y0_b);
[tc, yc] = ode45(@(t, y) eom(t, y, M, K), tspan, y0_c);
[td, yd] = ode45(@(t, y) eom(t, y, M, K), tspan, y0_d);

figure(1)
title('part a')
hold on
plot(ta,ya(:,1))
plot(ta,ya(:,2))
figure(2)
title('part b')
hold on
plot(tb,yb(:,1))
plot(tb,yb(:,2))
figure(3)
title('part c')
hold on
plot(tc,yc(:,1))
plot(tc,yc(:,2))
figure(4)
title('part d')
hold on
plot(td,yd(:,1))
plot(td,yd(:,2))

function dydt = eom(t, y, M, K)
    n = size(M, 1);
    q = y(1:n);       % Displacement
    q_dot = y(n+1:end); % Velocity
    acc = M \ (- K * q); % Acceleration
    
    dydt = [q_dot; acc]; % Combine velocity and acceleration
end