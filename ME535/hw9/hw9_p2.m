clear
clc
close all

m1 = 1;
m2 = 1;
k = 8;
L = 1;
g = 9.81;
c = 0.1;

omegacart = sqrt(k/m2);

M = [(1/3) * m1 * L, 0; 0, m2];
K = [(L^2 / 4) * k + m1 * g * L / 2, (-L/2) * k; (-L/2) * k, k];
C = [0, 0; 0, c];
f = 1;

F = [f * L; 0];

[vectors, values] = eig(K,M);

natural_frequencies = [sqrt(values(1,1)); sqrt(values(2,2))];

norm_vector1 = vectors(:,1) / sqrt(transpose(vectors(:,1)) * M * vectors(:,1));
norm_vector2 = vectors(:,2) / sqrt(transpose(vectors(:,2)) * M * vectors(:,2));

Phi = [norm_vector1, norm_vector2];

force_frequencies = linspace(0,natural_frequencies(end)*1.5,1000);

for i = 1:length(force_frequencies)
    X(i,:) = (K + 1i * force_frequencies(i) * C - force_frequencies(i)^2 * M)^(-1) * F;
end

other = sqrt(k/m2);

hold on
plot(force_frequencies,abs(X))
% plot(force_frequencies,atan2(imag(X),real(X)))
xline(omegacart)
xlabel('force frequency')
ylabel('transfer function')
legend(["\theta/F","x/F",""])
title('P2')
