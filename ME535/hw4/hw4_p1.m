clear
clc
close all

q0 = 1;
v0 = 1;
omega = 1;
F = 1;
M = 1;
T = linspace(0,100,1000);

Q = q0 - (v0 ./ omega) .* 1i - (F ./ (T .* M .* omega.^3)) .* 1i + exp(-1i .* omega .* T) .* (F ./ (T .* M .* omega.^3)) .* 1i - exp(-1i .* omega .* T) .* (F ./ (M .* omega.^2));

amplitude = abs(Q);

figure(1)
plot(T,amplitude)


figure(2)
other = atan(imag(Q) ./ real(Q));
plot(T,other)
xline(0)

abs(log(q0 / (F/(M*omega^2))) / (-1i * omega)) * omega
