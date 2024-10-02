clear
clc
close all

m = 10;
natural_frequency = (2 * pi) * 1e3;
driving_frequency_original = (2 * pi) * 1e3;
driving_frequency_95 = (2 * pi) * 0.95e3;
driving_frequency_105 = (2 * pi) * 1.05e3;
k = natural_frequency^2 * m;
amp_original = 2.4e-3;
F = -1200 * 1i;
damping_ratio = abs((F / k) * (2 * 1i * amp_original)^-1);

c = damping_ratio * 2 * sqrt(k * m);

r_95 = driving_frequency_95 / natural_frequency;
r_105 = driving_frequency_105 / natural_frequency;
r_original = driving_frequency_original / natural_frequency;
D_95 = (1 + 2 * 1i * damping_ratio * r_95 - r_95^2)^(-1);
D_105 = (1 + 2 * 1i * damping_ratio * r_105 - r_105^2)^(-1);
D_original = (1 + 2 * 1i * damping_ratio * r_original - r_original^2)^(-1);
X_95 = D_95 * F / k;
X_105 = D_105 * F / k;
X_original = D_original * F / k;

amp_95 = abs(X_95);
amp_105 = abs(X_105);
amp_original = abs(X_original);


t = linspace(0,.01,1000);
q_original = real(X_original .* exp(1i .* natural_frequency .* t));
q_95 = real(X_95 .* exp(1i .* driving_frequency_95 .* t));
q_95_book = 3.1173.*(10.^-5) .* sin(1900.*pi.*t - 0.01234);
q_105 = real(X_105 .* exp(1i .* driving_frequency_105 .* t));
q_105_book = 2.965 .* (10^-5) .* sin(2100 .* pi .* t - 3.129);

hold on
plot(t, q_original,'b-')
plot(t, q_95,'r-')
% plot(t,q_95_book,'r--')
plot(t, q_105,'k-')
% plot(t,q_105_book,'k--')
xlabel('time (s)')
ylabel('response (m)')
legend(["${\omega = 1 kHz}$","${\omega = 0.95 kHz}$","${\omega = 1.05 kHz}$"],'interpreter','latex','location','eastoutside')

disp("amplitude at 0.95 kHz: "+amp_95+" m")
disp("amplitude at 1.05 kHz: "+amp_105+" m")

% amplitude at 0.95 kHz: 3.1173e-05 m
% amplitude at 1.05 kHz: 2.9652e-05 m

