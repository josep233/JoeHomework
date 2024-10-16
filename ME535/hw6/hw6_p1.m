clear
clc
close all

l = 1;
k = 1000;
m = 1;
theta = deg2rad(45);
xs = linspace(-1,1,100);

syms X L THETA M K

L_i = L / cos(THETA);
h = sqrt((L^2 / cos(THETA)^2) - L^2);
L_f = sqrt(h^2 + (L + X)^2);
delta = L_f - L_i;
delta_squared = simplify(delta^2);

potential_energy = (1/2) * K * delta_squared;
K_matrix = simplify(diff(diff(potential_energy,X),X));

K_full = subs(K_matrix,{L,K,THETA}, {l, k, theta});

K_taylor = taylor(K_matrix,X);

% pretty(K_taylor)

K_of_X = subs(K_taylor,{L,K,THETA}, {l, k, theta});

for i = 1:length(xs)
    k_of_x(i) = subs(K_of_X, X, xs(i));
end

t_range = [0 1];
x_0_1 = 0.1;
x_0_2 = 0.5;

% [tout1,yout1] = ode45(@(t, x) eom(t, x, K_full, m), t_range, [x_0_1, 0]);
% [tout2,yout2] = ode45(@(t, x) eom(t, x, K_full, m), t_range, [x_0_2, 0]);

hold on
% plot(xs,k_of_x,'blue')
plot(xs,sqrt(k_of_x/m),'green')
% yline(k * cos(theta)^2,'red')

% hold on
% plot(tout1,yout1(:,2))
% plot(tout2,yout2(:,2))

function [theta_dot] = eom(t, x, k_func, m)        
    % Equations of Motion
    theta_dot(1,1) = x(2);

    syms X

    K_val = subs(k_func, X, x(1));

    theta_dot(2,1) = (1/m) * (-K_val * x(1));
end