clear
clc
close all

m1 = 0.5;
m2 = 1;
k = 3.2e3;
mu = 40;
A = 20e-3;
omega = 75;

X = (A * omega * mu / (k + 1i * omega * mu - m1 * omega^2));
F = (-A * 1i - mu * 1i * omega * X / (1i * omega * mu - m2 * omega^2)) * (1i * omega * mu - m2 * omega^2);

disp("magnitude of F for frequency = 75 rad/s: "+abs(F)+" N")
disp("phase of F for frequency = 75 rad/s: "+(atan(imag(F)/real(F))+(pi/2))+" rad")

omega = 85;

X = (A * omega * mu / (k + 1i * omega * mu - m1 * omega^2));
F = (-A * 1i - mu * 1i * omega * X / (1i * omega * mu - m2 * omega^2)) * (1i * omega * mu - m2 * omega^2);

disp("magnitude of F for frequency = 85 rad/s: "+abs(F)+" N")
disp("phase of F for frequency = 85 rad/s: "+(atan(imag(F)/real(F))+(pi/2))+" rad")

% magnitude of F for frequency = 75 rad/s: 104.8818 N
% phase of F for frequency = 75 rad/s: 3.1322 rad
% magnitude of F for frequency = 85 rad/s: 152.6335 N
% phase of F for frequency = 85 rad/s: 3.1351 rad