clear
clc
close all

%given information from Goland / Dr Allen
rho = 0.002378; %density - lbm/ft^3
L = 20; %length of wing - ft
c = 6; %chord - ft
b = c/2; %half chord -ft
elastic_center = 0.33 * c; %elastic center - ft
center_of_gravity = 0.43 * c; %center of gravity - ft
aerodynamic_center = 0.010 * c + c / 4; %aerodynamic center - ft
a = (elastic_center - b) / b; %parameter 'a' from slides - frac
mw = 0.746; %mass per length - slugs / ft
Sy = 0.447; % slugs-ft / ft (center of gravity - elastic_center) * mw
Iy = 1.943; %Inertia of wing about elastic axis - slugs-ft^2/ft
E = 1.4976e9; %elastic modulus - lb-ft^2
G = 5.616e8; %shear modulus - lb-ft^2
Iea = 1.58e-2; %area moment of inertia of wing in bending - ft^4
J = 4.25e-3; %area moment of inertia of wing in torsion - ft^4
x = (center_of_gravity - elastic_center);
A = mw / rho; %this is inferred, really.

%these mode functions came from the book
syms y real
alpha1 = 1.8571;
R = -(sin(alpha1) + sinh(alpha1))/(cos(alpha1) + cosh(alpha1));
psiB = sin(alpha1*y) - sinh(alpha1*y) + R*(cos(alpha1*y) - cosh(alpha1*y));
psiT = sin(pi*y/2);

%from hand deriviation, M
Mb = rho * A * double(int(psiB^2,0,1));
Mbt = rho * A * x * double(int(psiB*psiT,0,1));
Mt = rho * A * x^2 * double(int(psiT^2,0,1)) + Iy * double(int(psiT^2,0,1));
M = [Mb,Mbt;Mbt,Mt];

%from hand deriviation, K
Kb = (E * Iea / L^4) * double(int(diff(psiB,y,2)^2,0,1));
Kbt = x * (E * Iea / L^4) * double(int(diff(psiB,y,2) * diff(psiT,y,2),0,1));
Kt = x^2 * (E * Iea / L^4) * double(int(diff(psiT,y,2)^2,0,1)) + (G * J / L^2) * double(int(diff(psiT,y,1)^2,0,1));
K = [Kb,Kbt;Kbt,Kt];

%estimate frequencies at v = 0
actual_uncoupled_frequencies = [1.8755^2*sqrt(E*Iea/(mw*L^4)); (pi/(2*L))*sqrt(G*J/Iy)];
calculated_uncoupled_frequencies = sqrt(eig(K,M));

% loop over ks
ks = [0.005:0.1:1,2:2:20]; lam = zeros(2,length(ks)); U = lam; g = lam; omega = lam;
for ki = 1:length(ks)
    k = ks(ki);

    %theodorsen stuff, from slides
    C = besselh(1,2,k)./(besselh(1,2,k)+1i*besselh(0,2,k)); % F+i*G
    L_h = -((- k + C*2i))/k;
    L_alpha = -(b*(2*C + k*1i + C*k*1i + a*k^2 - C*a*k*2i))/k^2;
    M_h = -(b*(a*k + C*a*2i - C*1i))/k ;
    M_alpha = (b^2*(8*C - k*4i - 16*C*a + C*k*4i + a*k*8i + k^2 - 8*a^2*k^2 - C*a*k*16i + C*a^2*k*16i))/(8*k^2);
    
    %from hand deriviation, Q
    Q_h = double(int(pi * rho * b^2 * L_h * psiB^2,0,1));
    Q_halpha = double(int(x * pi * rho * b^2 * L_alpha * psiB * psiT,0,1));
    Q_alphah = double(int(x * pi * rho * b^2 * M_h * psiT * psiB,0,1));
    Q_alpha = double(int(x^2 * pi * rho * b^2 * M_alpha * psiT^2,0,1));
    Q = [Q_h, Q_halpha; Q_alphah, Q_alpha];

    %EVP
    [Phi,Lam] = eig(K,(M+Q));
    lam(:,ki) = diag(Lam);
    omega(:,ki) = sqrt(abs(lam(:,ki))); % in rad/s
    zts(:,ki) = imag(sqrt(lam(:,ki)))./sqrt(abs(lam(:,ki)));
    U(:,ki) = omega(:,ki)*b/k;
end

%I got all of this stuff from Dr Allen's code
search_ind = find(U(2,:) < 600);
fl_mode_ind = 2;
fl_U_ind = find(zts(fl_mode_ind,:) < -0.03,1,'last');
U_flutter = U(fl_mode_ind,fl_U_ind)
omega_flutter = omega(fl_mode_ind,fl_U_ind);
omega_c = omega; zts_c = zts; U_c = U;
for k = size(omega_c,2):-1:2
    if abs(omega_c(1,k) - omega_c(1,k-1)) > 5 && abs(omega_c(2,k) - omega_c(1,k-1)) <= 5
        omega_c([1,2],k-1:-1:1) = omega_c([2,1],k-1:-1:1);
        zts_c([1,2],k-1:-1:1) = zts_c([2,1],k-1:-1:1);
        U_c([1,2],k-1:-1:1) = U_c([2,1],k-1:-1:1);
    end
end
figure(1)
subplot(2,1,1);
plot(U_c.',omega_c.',U(fl_mode_ind,fl_U_ind),omega_c(fl_mode_ind,fl_U_ind),'*'); grid on;
ylabel('\bfFrequencies (rad/s)');
set(get(gca,'Children'),'LineWidth',2)
subplot(2,1,2);
plot(U_c.',zts_c.',U_c(fl_mode_ind,fl_U_ind),zts_c(fl_mode_ind,fl_U_ind),'*'); grid on;
xlabel('\bfSpeed (ft/s)'); ylabel('\bfDamping Ratio');

