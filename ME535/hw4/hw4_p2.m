clear
clc
close all

L = 50;
h = 0.4572;
w = 1.2192;
rho = 7800;
damping_ratio = 0.01;
E = 210e9;
I = (w * h^3) / 12;
k_beam = 48 * E * I / L^3;
F0 = 1000;
m = (1/3) * L * w * h * rho;
c = damping_ratio * 2 * sqrt(k_beam*m);
natural_frequency = sqrt(k_beam / m);

force_frequency = linspace(0,8*pi,1000);

X = F0 ./ (k_beam + 1i .* force_frequency .* c - m .* force_frequency.^2);

t = linspace(0,100,1000);

for i = 1:length(force_frequency)
    [max_response(i), indices(i)] = max(real(X(i) .* exp(1i .* natural_frequency .* t)));
end

[maxresponse, maxindex] = max(max_response);

best_frequency = force_frequency(maxindex);

plot(force_frequency,max_response)

disp("frequency to cause maximum displacement: "+string(best_frequency)+" rad/s"+"("+(2*pi/best_frequency)+" Hz)")
disp("displacement at this frequency: "+100*maxresponse+" cm")

% frequency to cause maximum displacement: 3.2957 rad/s(1.9065 Hz)
% displacement at this frequency: 6.165 cm