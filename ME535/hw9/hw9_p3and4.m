clear
clc
close all

a = 100e-6;
t = 7e-6;
E = 170E9;
rho = 2300;
w = 10e-6;
h = 7e-6;
omega1 = 11.6e3 * 2 * pi;
omega2 = 18.8e3 * 2 * pi;
natural_frequencies = [omega1, omega2];
zeta = [0.001; 0.001];
F = [0;1];

m = a^2 * t * rho;
kappa = sqrt(1/12) * t;
Ig = m * kappa^2;

M = [(Ig / a^2) + (1 / 4) * m, (1 / 4) * m - (Ig / a^2); (1 / 4) * m - (Ig / a^2), (Ig / a^2) + (1 / 4) * m];
K = [(E * w * h^3) / (2), 0; 0, (E * w * h^3) / (2)];

[vectors, values] = eig(K,M);

lambda = [values(1,1); values(2,2)];

norm_vector1 = vectors(:,1) / sqrt(transpose(vectors(:,1)) * M * vectors(:,1));
norm_vector2 = vectors(:,2) / sqrt(transpose(vectors(:,2)) * M * vectors(:,2));

Phi = [norm_vector1, norm_vector2];

L = (lambda./natural_frequencies.^2').^(1/3);

list_frequencies = linspace(0,natural_frequencies(end) * 1.2, 1000);

for i = 1:height(K)
    for j = 1:length(list_frequencies)
        force_frequency = list_frequencies(j);
        Y(i,j) = sum(Phi(:,i) * ((transpose(Phi(:,i)) * F)/(natural_frequencies(i)^2 + 2 * 1i * zeta(i) * natural_frequencies(i) * force_frequency - force_frequency^2)));
    end
end

X = Phi * Y;

figure(1)
plot(list_frequencies,abs(Y))
xline(natural_frequencies)
xlabel('frequency')
ylabel('amplitude')
title('amplitude')

figure(2)
plot(list_frequencies,atan2(imag(X),real(X)))
xline(natural_frequencies)
xlabel('frequency')
ylabel('phase')
title('phase')


% norm_vector1 = vectors(:,1) / sqrt(transpose(vectors(:,1)) * M * vectors(:,1));
% norm_vector2 = vectors(:,2) / sqrt(transpose(vectors(:,2)) * M * vectors(:,2));
% 
% Phi = [norm_vector1, norm_vector2];

