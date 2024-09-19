% Example 2-12 using Matlab - Analytical and Numerical
clear all; close all

% Parameter values
M = 1; K = 1; wn = sqrt(K/M); F0 = 1;
T = 10;
q_0 = 0; q_dot_0 = 0; % zero initial conditions! qIC=0

% Time Vector
ts = [0:0.5:50];

% Forcing - sum of step and ramp
% note h(t) written as (t>0) in Matlab
F = -F0*(ts>0) + (F0/T).*(ts-T).*(ts-T>0);

%% Analytical Solution

% Unit step and ramp responses from Ginsberg - includes particular and 
% complimentary solutions
qs = @(t) (1/(M*wn^2))*(1-cos(wn*t))*(t>0)
qr = @(t) (1/(M*wn^3))*(wn*t-sin(wn*t))*(t>0)
    % Note, when creating an anonymous function like this, it sees any
    % variables that are already defined (like M and wn) and uses the
    % values in the Matlab workspace WHEN IT IS CREATED for those.

% Response is a sum of step and ramp responses
for k = 1:length(ts)
    q(k) = -F0*qs(ts(k))+(F0/T)*qr(ts(k)-T);
end

figure(1)
subplot(2,1,1)
plot(ts,F); grid on;
title('Forcing F(t)');
subplot(2,1,2);
plot(ts,q); grid on;
title('Response q(t)');

%% Solution using ODE45

% Define equations of motion in eom_2_12.m
% Note - ode45 requires only the time span, not the whole time vector
tic
[tout,yout] = ode45('eom_2_12',[ts(1),ts(end)],[q_0; q_dot_0]);
t_ode = toc

q_ode = yout(:,1); % the first of the y variables is q(t), the second is q_dot(t)

% Add red dots to plot above
hold on; plot(tout,q_ode,'r.'); hold off;

%% Another option - lsim - requires control system toolbox
% Define Linear State Space System Matrices
 A = [0, 1;
     -M/K, 0];
 B = [0; 1];
 C = [1, 0]; 
 D = 0;
 % Collect into a matlab State Space (SS) model
%  sys = ss(A,B,C,D);
% 
%  tic
%  q_lsim = lsim(sys,F,ts,[q_0; q_dot_0]);
%  t_lsim = toc
% 
%  % Add green circles to plot above for lsim result
% hold on; plot(ts,q_lsim,'go'); hold off;
% legend('Analytical','ode45','lsim','Location','NorthWest');