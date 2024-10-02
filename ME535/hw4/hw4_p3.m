clear
clc
close all

F = 1;
r = linspace(0,2,100);
damping_ratio = 0.1; %can also change this to 0.02
X = F ./ (1 + 2 .* 1i .* damping_ratio .* r - r.^2);

amplitude = abs(X);
phase = atan(imag(X)./real(X));

figure(1);
plot(r,amplitude)
xlabel("${r = \omega / \omega_{n}}$",'interpreter','latex')
ylabel("${|\hat{A}| (m/s^{2})}$",'interpreter','latex')

figure(2);
plot(r,phase);
xlabel("${r = \omega / \omega_{n}}$",'interpreter','latex')
ylabel("${\phi (rad)}$",'interpreter','latex')
