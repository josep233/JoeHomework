clear
clc
close all

syms x L A0 rho E

for i = 1:4
    psi(i) = sin((i * pi * x)/L);
    % psi(i) = (x/L) * (1 - x/L) * sin((i * pi * x)/L);
end

A = A0 * (1 + x/L);

for i = 1:4
    for j = 1:4
        M(i,j) = int(psi(i) * psi(j) * rho * A, x, 0, L);
        K(i,j) = int(diff(psi(i),x,1) * diff(psi(j),x,1) * E * A, x, 0, L);
    end
end

M_num = double(subs(M,[L,rho,A0,E],[1,1,1,1]));
K_num = double(subs(K,[L,rho,A0,E],[1,1,1,1]));

