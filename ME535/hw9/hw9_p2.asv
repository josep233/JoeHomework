clear
clc
close all

syms m1 m2 g L k c omegas

M = [(1/3) * m1 * L, 0; 0, m2];
K = [(L^2 / 4) * k + m1 * g * L / 2, (-L/2) * k; (-L/2) * k, k];
C = [0, 0; 0, c];

foil = simplify(expand((K(1,1) - M(1,1) * omegas) * (K(2,2) - M(2,2) * omegas) - (K(1,2) - M(1,2) * omegas)^2));

simplify(expand(solve(foil == 0, omegas)))

