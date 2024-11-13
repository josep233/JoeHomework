clear
clc
close all

syms k m t

M = [0.41, 0.09; 0.09, 0.41];
K = [1.5, 0; 0, 1];
q_0 = [1; 0];
qdot_0 = [0; 0];

[vectors, values] = eig(K,M);

natural_frequencies = [sqrt(values(1,1)),sqrt(values(2,2))];
vector1 = vectors(:,1) / vectors(1,1);
vector2 = vectors(:,2) / vectors(1,2);

Phi1 = vector1 ./ (sqrt(transpose(vector1) * M * vector1));
Phi2 = vector2 ./ (sqrt(transpose(vector2) * M * vector2));

Phi = [Phi2,Phi1];

eta_0 = transpose(Phi) * M * q_0;
etadot_0 = transpose(Phi) * M * qdot_0;

q = Phi * eta_0;

t = linspace(0,50,1000);
y1 = 0.838 * cos(2 * t) + 0.159 * cos(1.5 * t);
y2 = -0.446 * cos(2 * t) + 0.448 * cos(1.5 * t);

hold on
plot(t, y1, 'blue')
plot(t,y2,'red')
title('problem 4')
xlabel("t'")
ylabel("y(t')k/mg")
