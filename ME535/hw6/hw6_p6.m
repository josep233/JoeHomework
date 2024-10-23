clear
clc
close all

M = 80;
g = 9.81;
r = 50e-3;
force_frequency = 145 * 2 * pi / 60;
static_displacement = 40e-3;

k = M * g / static_displacement;
natural_frequency = sqrt(k/M);

damping_ratio = (2 * force_frequency * natural_frequency)^(-1) * tan(deg2rad(75)) * (natural_frequency^2 - force_frequency^2);
damping_ratio2 = (2 * force_frequency * natural_frequency)^(-1) * (tan(deg2rad(15 - 90))) * (natural_frequency^2 - force_frequency^2);

fraction = (-1i * force_frequency^2 / M) / (-force_frequency^2 + 2 * 1i * natural_frequency * force_frequency * damping_ratio + natural_frequency^2);
mult = abs(fraction);

m = 10e-3 / (mult * r);