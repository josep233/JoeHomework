clear
clc
close all

syms x L rho A E I

% m = (1/2) * rho * A * L;
% k = 10 * E * I / L^3;

m = 0;
k = 0;

alpha(1) = 1.8751;
alpha(2) = 4.6941;
alpha(3) = (2 * 3 - 1) * pi / 2;

for i = 1:3
    R(i) = -((sin(alpha(i)) + sinh(alpha(i))) / (cos(alpha(i)) + cosh(alpha(i))));
    psi(i) = sin(alpha(i) .* x/L) - sinh(alpha(i) .* x/L) + R(i) .* (cos(alpha(i) .* x/L) - cosh(alpha(i) .* x/L));
end

% for i = 1:3
%     psi(i) = (x/L)^(i);
% end

for i = 1:3
    for j = 1:3
        M(i,j) = int(psi(i) * psi(j) * rho * A, x, 0, L) + sum(m * subs(psi(i),x,L) * subs(psi(j),x,L));
        K(i,j) = int(E * I * diff(psi(i),x,2) * diff(psi(j),x,2),x,0,L) + sum(k * subs(psi(i),x,L) * subs(psi(j),x,L));
    end
end

M
K

M_num = double(subs(M,[L rho A E I],[1, 1, 1, 1, 1]));
K_num = double(subs(K,[L rho A E I],[1, 1, 1, 1, 1]));

[vectors, values] = eig(K_num,M_num);

natural_frequencies = [sqrt(values(1,1)); sqrt(values(2,2)); sqrt(values(3,3))];

% norm_vector1 = vectors(:,1) / max(vectors(:,1));
% norm_vector2 = vectors(:,2) / max(vectors(:,2));
% norm_vector3 = vectors(:,3) / max(vectors(:,3));

norm_vector1 = vectors(:,1) / sqrt(transpose(vectors(:,1)) * M_num * vectors(:,1));
norm_vector2 = vectors(:,2) / sqrt(transpose(vectors(:,2)) * M_num * vectors(:,2));
norm_vector3 = vectors(:,3) / sqrt(transpose(vectors(:,3)) * M_num * vectors(:,3));

Phi = [norm_vector1, norm_vector2, norm_vector3];

xs = linspace(0,1,100);

alpha(1) = 1.8751;
alpha(2) = 4.6941;
alpha(3) = (2 * 3 - 1) * pi / 2;

for i = 1:3
    R(i) = -((sin(alpha(i)) + sinh(alpha(i))) / (cos(alpha(i)) + cosh(alpha(i))));
    psi_analytic(i,:) = sin(alpha(i) .* xs) - sinh(alpha(i) .* xs) + R(i) .* (cos(alpha(i) .* xs) - cosh(alpha(i) .* xs));
end

for i = 1:3
    for j = 1:length(xs)
        modefunction(i,j) = sum(double(subs(Phi(:,i) * psi(i),[x,L],[xs(j),1])));
    end
end

hold on
plot(xs,psi_analytic)
plot(xs,modefunction,'black')
