clear
clc
close all

dr = linspace(0,1,100);

for i = 1:length(dr)
    m = 1;
    natural_frequency = 5 * 2 * pi;
<<<<<<< Updated upstream
    k = (natural_frequency^2);
    v = 4;
    g = 9.81;
    % g = 0;
=======
    k = (natural_frequency^2) * m;
    v = -4;
    g = 9.81;
>>>>>>> Stashed changes
    damping_ratio = dr(i);
    c = damping_ratio * 2 * sqrt(k*m);
    t = linspace(0,2,100);
    z_0 = 0;
    z_dot_0 = v;
    
<<<<<<< Updated upstream
    % z_func = @(t) exp(-damping_ratio * natural_frequency * t) + (v/(natural_frequency * sqrt(1 - damping_ratio^2))) * sin(natural_frequency * sqrt(1 - damping_ratio^2) * t);
    
    z_ana = exp(-damping_ratio .* natural_frequency .* t) .* (v./(natural_frequency .* sqrt(1 - damping_ratio.^2))) .* sin(natural_frequency .* sqrt(1 - damping_ratio.^2) .* t);
    
    % hold on
    % plot(t,z_ana,'black')
    
    [tout,yout] = ode45(@(t, z) eom_2_12(z, m, k, c, g), [t(1),t(end)],[z_0; z_dot_0]);
    
    % plot(tout,yout(:,1),'red')
    % ylabel('z (m)','interpreter','latex')
    % xlabel('time (s)','interpreter','latex')
=======
    z_func = @(t) exp(-damping_ratio * natural_frequency * t) .* (v/(natural_frequency * sqrt(1 - damping_ratio^2))) * sin(natural_frequency * sqrt(1 - damping_ratio^2) * t);
    
    z = exp(-damping_ratio .* natural_frequency .* t) .* (v./(natural_frequency .* sqrt(1 - damping_ratio.^2))) .* sin(natural_frequency .* sqrt(1 - damping_ratio.^2) .* t);

    hold on
    % plot(t,z)
    
    [tout,yout] = ode45(@(t, z) eom_2_12(z, m, k, c, g), [t(1),t(end)],[z_0; z_dot_0]);
    
    % plot(tout,yout(:,1))
    ylabel('z (m)','interpreter','latex')
    xlabel('time (s)','interpreter','latex')
>>>>>>> Stashed changes
    
    sum_force = (yout(:,2) * c + k * yout(:,1));
    sum_force_abs = abs(sum_force);
    max_abs_force(i) = max(sum_force_abs);

    % hold on
    % plot(tout,sum_force)
    % pause;
<<<<<<< Updated upstream
    % close all
=======
>>>>>>> Stashed changes
end

[minimum, idx] = min(max_abs_force);
best_damping_ratio = dr(idx);

hold on
plot(dr,max_abs_force)

function [z_dot] = eom_2_12(z, m, k, c, g)    
    % Forcing - sum of step and ramp
    % F = 0;
    F = m * g;
    
    % Equations of Motion
    z_dot(1,1) = z(2);
    z_dot(2,1) = (1/m) * (F - c * z(2) - k * z(1));
end