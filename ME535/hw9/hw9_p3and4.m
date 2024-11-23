clear
clc
close all

a = 100e-6;
t = 7e-6;
E = 170E9;
rho = 2300;
w = 10e-6;
h = 2e-6;
omega1 = 11.6e3 * 2 * pi;
omega2 = 18.8e3 * 2 * pi;
natural_frequencies = [omega1, omega2];
zeta = [0.001; 0.001];
F = [0;1];
L = 251e-6;
k=3*E*(w*h^3/12)/L^3;

m = a^2 * t * rho;
% kappa = sqrt(1/12) * t;
Ig = (1/12) * m * a^2;

M = [(Ig / a^2) + (1 / 4) * m, (1 / 4) * m - (Ig / a^2); (1 / 4) * m - (Ig / a^2), (Ig / a^2) + (1 / 4) * m];
K = [(E * w * h^3) / (2), 0; 0, (E * w * h^3) / (2)];
K = [2*k, 0; 0, 2*k];

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

fn=sqrt(diag(lambda))/2/pi;
C=M*Phi*diag(2*(fn*2*pi).*zeta).*Phi.'*M;
for k=1:length(list_frequencies)
H(:,k)=(K+1i*list_frequencies(k)*C-list_frequencies(k)^2*M)\F;
end
X = Phi * Y;

semilogy(list_frequencies/2/pi,abs(H)); xlabel('Frequency (Hz)');
legend('H_1(\omega)','H_2(\omega)');

% figure(1)
% semilogy(list_frequencies./(2*pi),abs(H))
% xline(natural_frequencies)
% xlabel('frequency')
% ylabel('amplitude')
% title('amplitude')

% figure(2)
% semilogy(list_frequencies,atan2(imag(X),real(X)))
% xline(natural_frequencies)
% xlabel('frequency')
% ylabel('phase')
% title('phase')


% norm_vector1 = vectors(:,1) / sqrt(transpose(vectors(:,1)) * M * vectors(:,1));
% norm_vector2 = vectors(:,2) / sqrt(transpose(vectors(:,2)) * M * vectors(:,2));
% 
% Phi = [norm_vector1, norm_vector2];

