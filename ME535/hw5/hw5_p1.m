clear
clc
close all

m = 450;
damping_ratio = 0.03;
force_frequency = 60 * pi;
force_amplitude = 20e3;
min_force = 2e3;

q = @(k) abs(force_amplitude / (-force_frequency^2 * m + k + 1i * force_frequency * 2 * damping_ratio * sqrt(k*m)));
all = @(k) abs(k * q(k) + 2 * 1i * damping_ratio * sqrt(k*m) * force_frequency * q(k)) - min_force;

an = fzero(all,100);

k = 1.428e6;

natural_frequency = sqrt(k/m);

new_frequency = linspace(0,60*pi,1000);

X = force_amplitude ./ (-new_frequency.^2 .* m + k + 1i .* new_frequency .* 2 .* damping_ratio .* sqrt(k.*m));

amplitude = abs(X);
max(amplitude)

mnew = m + linspace(0,10000,1000);
knew = natural_frequency^2 .* mnew;

for i = 1:length(mnew)

Xnew(i,:) = force_amplitude ./ (-new_frequency.^2 .* mnew(i) + knew(i) + 1i .* new_frequency .* 2 .* damping_ratio .* sqrt(knew(i).*mnew(i)));
max_amplitude(i) = max(abs(Xnew(i,:)));
end

hold on
plot(mnew,max_amplitude)
yline(0.01)

knew(end)
mnew(end)