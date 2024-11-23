clear
clc
close all

M = [600, 400, 200; 400, 1200, 0; 200, 0 800];
K = [300, 0, -200; 0, 500, 300; -200, 300, 700] .* 1e3;
C = [500, 300, -400; 300, 900, 600; -400, 600, 1300];
C_light = diag(C) .* eye(3);

amps = [200; 0; 0];
% force_frequencies = [16; 0; 0];
force_frequencies = [16; 16; 16];

q_0 = [0; 0; 0];
qdot_0 = [0; 0; 0];

[vectors, values] = eig(K,M);

natural_frequencies = [sqrt(values(1,1)); sqrt(values(2,2)); sqrt(values(3,3))];

norm_vector1 = vectors(:,1) / sqrt(transpose(vectors(:,1)) * M * vectors(:,1));
norm_vector2 = vectors(:,2) / sqrt(transpose(vectors(:,2)) * M * vectors(:,2));
norm_vector3 = vectors(:,3) / sqrt(transpose(vectors(:,3)) * M * vectors(:,3));

Phi = [norm_vector1, norm_vector2, norm_vector3];

zeta = diag((transpose(Phi) * C * Phi)) .* (2 .* natural_frequencies).^(-1);
damped_frequencies = natural_frequencies .* sqrt(1 - zeta.^2);

eta_0 = transpose(Phi) * M * q_0;
etadot_0 = transpose(Phi) * M * qdot_0;

t = linspace(0,140,10000);

eta_homogeneous = exp(-zeta .* natural_frequencies .* t) .* (eta_0 .* cos(damped_frequencies .* t)...
    + ((etadot_0  + zeta .* natural_frequencies .* eta_0) ./ damped_frequencies) .* sin(damped_frequencies .* t));

% eta1_homogeneous = exp(-zeta(1) .* natural_frequencies(1) .* t) .* (eta_0(1) .* cos(damped_frequencies(1) .* t) + ((etadot_0(1)  + zeta(1) .* natural_frequencies(1) .* eta_0(1)) / damped_frequencies(1)) .* sin(damped_frequencies(1) .* t));
% eta2_homogeneous = exp(-zeta(2) .* natural_frequencies(2) .* t) .* (eta_0(2) .* cos(damped_frequencies(2) .* t) + ((etadot_0(2)  + zeta(2) .* natural_frequencies(2) .* eta_0(2)) / damped_frequencies(2)) .* sin(damped_frequencies(2) .* t));
% eta3_homogeneous = exp(-zeta(3) .* natural_frequencies(3) .* t) .* (eta_0(3) .* cos(damped_frequencies(3) .* t) + ((etadot_0(3)  + zeta(3) .* natural_frequencies(3) .* eta_0(3)) / damped_frequencies(3)) .* sin(damped_frequencies(3) .* t));

% for i = 1:length(t)
%     for j = 1:3
%     eta_forced(j,i) = Phi(:,j) * amps(j) * ((natural_frequencies(j)^2 - force_frequencies(j)^2)^2 + 4 * zeta(j)^2 ...
%         * natural_frequencies(j)^2 * force_frequencies(j)^2)^(-1);
%     end
% end

for i = 1:3
ETA(i,:) = Phi(i,1) * amps(1) * ((natural_frequencies(i)^2 - force_frequencies(i)^2)^2 + 4 .* zeta(i)^2 .* natural_frequencies(i)^2 .* ...
    force_frequencies(i)^2)^(-1) .* ((natural_frequencies(i)^2 - force_frequencies(i)^2) .* cos(force_frequencies(i) .* t) ...
    + 2 .* zeta(i) * natural_frequencies(i) * force_frequencies(i) .* sin(force_frequencies(i) .* t) - exp(-zeta(i) ...
    .* natural_frequencies(i) .* t) .* ((natural_frequencies(i)^2 - force_frequencies(i)^2) .* cos(damped_frequencies(i) .* t) ...
    + ((zeta(i) .* natural_frequencies(i) .* (natural_frequencies(i) + force_frequencies(i)^2)./damped_frequencies(i)).*sin(damped_frequencies(i) .* t))));
end

% ETA = eta_forced + eta_homogeneous;

% ETA = Phi * ETA;

Q = Phi * ETA;

plot(t,ETA(3,:))

tau1 = 110;
tau2 = 32;
tau3 = 7;

tau = ((4)./(zeta .* natural_frequencies));