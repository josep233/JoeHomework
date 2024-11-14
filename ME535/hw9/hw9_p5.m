clear
clc
close all

E = 210e9;
L = 1;
A = 0.0004;
rho = 7800;
N = 5;
k1 = 2 * N * E * A / L;
k = N * E * A / L;
m = rho * A * L / N;
F = [0; 0; 0; 0; 1];
gamma = 0.04;

M = eye(5) * m;
K = [k1 + k, -k, 0, 0, 0; -k, 2*k, -k, 0, 0; 0, -k, 2*k, -k, 0; 0, 0, -k, 2*k, -k; 0, 0, 0, -k, k];

[vectors, values] = eig(K,M);

natural_frequencies = [sqrt(values(1,1)); sqrt(values(2,2)); sqrt(values(3,3)); sqrt(values(4,4)); sqrt(values(5,5))];

alpha = 0.02 * 2 * natural_frequencies(3);
zeta_a = alpha ./ (2 .* natural_frequencies);
zeta_b = 0.02 .* ones(size(natural_frequencies));
force_frequency_list = [0; 0; 0; 0; 1];
% force_frequency_list = [1; 1; 1; 1; 1];

norm_vector1 = vectors(:,1) / sqrt(transpose(vectors(:,1)) * M * vectors(:,1));
norm_vector2 = vectors(:,2) / sqrt(transpose(vectors(:,2)) * M * vectors(:,2));
norm_vector3 = vectors(:,3) / sqrt(transpose(vectors(:,3)) * M * vectors(:,3));
norm_vector4 = vectors(:,4) / sqrt(transpose(vectors(:,4)) * M * vectors(:,4));
norm_vector5 = vectors(:,5) / sqrt(transpose(vectors(:,5)) * M * vectors(:,5));

Phi = [norm_vector1, norm_vector2, norm_vector3, norm_vector4, norm_vector5];

list_frequencies = linspace(0,natural_frequencies(5) * 1.2, 1000);

for i = 1:height(K)
    for j = 1:length(list_frequencies)
        force_frequency = list_frequencies(j);
        zeta_c = (1./2) .* gamma .* natural_frequencies ./ force_frequency;
        Y_a(i,j) = sum(Phi(:,i) * ((transpose(Phi(:,i)) * F)/(natural_frequencies(i)^2 + 2 * 1i * zeta_a(i) * natural_frequencies(i) * force_frequency - force_frequency^2)));
        Y_b(i,j) = sum(Phi(:,i) * ((transpose(Phi(:,i)) * F)/(natural_frequencies(i)^2 + 2 * 1i * zeta_b(i) * natural_frequencies(i) * force_frequency - force_frequency^2)));
        Y_c(i,j) = sum(Phi(:,i) * ((transpose(Phi(:,i)) * F)/(natural_frequencies(i)^2 + 2 * 1i * zeta_c(i) * natural_frequencies(i) * force_frequency - force_frequency^2)));
    end
end

X_a = Phi * Y_a;
X_b = Phi * Y_b;
X_c = Phi * Y_c;

figure(1)
hold on
plot(list_frequencies, abs(X_a(5,:)),'red')
plot(list_frequencies, abs(X_b(5,:)),'blue')
plot(list_frequencies, abs(X_c(5,:)),'black')
xline(natural_frequencies)
xlabel('frequency')
ylabel('amplitude')
legend(["a","b","c"])
title('amplitude')

% figure(2)
% hold on
% plot(list_frequencies, atan2(imag(X_a(5,:)),real(X_a(5,:))),'red')
% plot(list_frequencies, atan2(imag(X_b(5,:)),real(X_b(5,:))),'blue')
% plot(list_frequencies, atan2(imag(X_c(5,:)),real(X_c(5,:))),'black')
% xline(natural_frequencies)
% xlabel('frequency')
% ylabel('phase')
% legend(["a","b","c"])
% title('phase')


