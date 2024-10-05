clear
clc
close all

m = 450;
damping_ratio = 0.03;
force_frequency = 60 * pi;
force_amplitude = 20e3;

q = @(k) abs(force_amplitude / (-force_frequency^2 * m + k + 1i * force_frequency * 2 * damping_ratio * sqrt(k/m) * m));
all = @(k) abs(k * q(k) + 2 * damping_ratio * sqrt(k/m) * m * force_frequency * q(k)) - 2e3;

an = fzero(all,1e6);