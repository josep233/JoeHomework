clear
clc
close all

const = 0.05;

M = (1/3) .* [1, 0, 0; 0, 1, 0; 0, 0, 1];
K = [1, -1, 0; -1, 2, -1; 0, -1, 1] .* const + eye(3) * (1/2);

q_0 = [0; 0; 0];
qdot_0 = [0; 2; 0];

[vectors, frequencies] = eig(K, M);

natural_frequencies = real([sqrt(frequencies(1,1)), sqrt(frequencies(2,2)), sqrt(frequencies(3,3))]);

natural_frequencies((natural_frequencies) == 0) = 1e-18;

norm_vector1 = vectors(:,1) / sqrt(transpose(vectors(:,1)) * M * vectors(:,1));
norm_vector2 = vectors(:,2) / sqrt(transpose(vectors(:,2)) * M * vectors(:,2));
norm_vector3 = vectors(:,3) / sqrt(transpose(vectors(:,3)) * M * vectors(:,3));

Phi = [norm_vector1,norm_vector2,norm_vector3];

eta_0 = transpose(Phi) * M * q_0;
etadot_0 = transpose(Phi) * M * qdot_0;

t = linspace(0,80,1000);

response1 = (1/natural_frequencies(1)) * etadot_0(1) * sin(natural_frequencies(1) * t);
response2 = (1/natural_frequencies(2)) * etadot_0(2) * sin(natural_frequencies(2) * t);
response3 = (1/natural_frequencies(3)) * etadot_0(3) * sin(natural_frequencies(3) * t);

ETA = [response1; response2; response3];

Q = Phi * ETA;

hold on
plot(t,Q)
title('problem 5')
xlabel("t'")
ylabel("y(t')")



