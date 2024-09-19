function [xdot] = eom_2_12(t,x)
    % Equations of motion for Ginsberg example 2.12
    
    % Note: x1 = x(1), x2 = x(2), ...
    
    % Parameter values
    M = 1; K = 1; wn = sqrt(K/M); F0 = 1; T = 10; C = 0;
    
    % Forcing - sum of step and ramp
    F = -F0*(t>0)+(F0/T).*(t-T).*(t-T>0);
    
    % Equations of Motion
    xdot(1,1) = x(2);
    xdot(2,1) = F/M-(K/M)*x(1)-(C/M)*x(2);