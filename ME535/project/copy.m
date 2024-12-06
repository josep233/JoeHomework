clear
clc
close all

Mbactual = 1.3843;
Mtactual = 0.9715;
Mbtactual = -0.4128;
Kbactual = 3.3926e3;
Ktactual = 7.3615e3;

rho = 0.002378;
% rho = 0.0025;

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

A = mw / rho;

%Ritz
N=1;
Ndof = 2*N;

syms y real
% Find alpha values for mode functions
for jj = 1:1:N
    r = (2*jj-1)*pi/2;
    al_v(jj) = fzero(@(y)cos(y)+1/cosh(y),r);
end

for jj=1:N
    R = (sinh(al_v(jj))+sin(al_v(jj)))/(cosh(al_v(jj))+cos(al_v(jj)));
    PsiB(jj)=sin(al_v(jj)*y)-sinh(al_v(jj)*y)-R*(cos(al_v(jj)*y)-cosh(al_v(jj)*y));
    PsiT(jj)=sin((2*jj-1)*pi*y/2);
    % Psi_pw(jj)=al_v(jj)*(cos(al_v(jj)*y)-cosh(al_v(jj)*y)+R*(sin(al_v(jj)*y)+sinh(al_v(jj)*y)));
end

Mb = rho * A * double(int(PsiB^2,0,1));
Mbt = -rho * A * x * double(int(PsiB*PsiT,0,1));
Mt = rho * A * x^2 * double(int(PsiT^2,0,1)) + Iy * double(int(PsiT^2,0,1));
% Mt = Iy * double(int(PsiT^2,0,1)) + rho * J * double(int(PsiT^2,0,1));

Kb = (E * Iea / L^4) * double(int(diff(PsiB,y,2)^2,0,1));
Kbt = -x * (E * Iea / L^4) * double(int(diff(PsiB,y,2) * diff(PsiT,y,2),0,1));
Kt = x^2 * (E * Iea / L^4) * double(int(diff(PsiT,y,2)^2,0,1)) + (G * J / L^2) * double(int(diff(PsiT,y,1)^2,0,1));

M = [Mb,Mbt;Mbt,Mt];
K = [Kb,Kbt;Kbt,Kt];

% Zero stores
nstores = 3;
d = 0.2*b*[1,1,1]; % 0.2 is x_alpha or how far the CG is behind the elastic axis
    % See p. 533 of BAH
mstores = 0*mw*[1,1,1];  % mass of store (assumes 2 ft width)
Istores = 0*Iy + mstores.*d.^2; % Inertia of store
Sstores = 0*mstores.*d;     % slugs-ft/ft - Store static imbalance
astores = 0*[1,1,1];     % ft - a_s*b_s = store CG Location relative to elastic axis
cstores = 0*[1,1,1];    % ft - store chord (see above, length of store in wing chord direction)
bstores = 0*cstores/2;  % ft - store chord (length)
ystores = [20,10,15]./L;    % nondim - store position.
swidth = 1; % Store width??

% Compute Frequencies for V = 0
disp('Frequencies for Uncoupled Bending and Torsion: [wb1,wt1]')
[1.8755^2*sqrt(E*Iea/(mw*L^4)); (pi/(2*L))*sqrt(G*J/Iy)]
disp('Frequencies for U = 0 (rad/s)');
wnsq_U0 = sqrt(eig(K,M))

% K-method - loop over various k's
ks = [0.005:0.005:1,2:2:20]; lam = zeros(Ndof,length(ks)); U = lam; g = lam; omega = lam;
for ki = 1:length(ks)
    k = ks(ki);

% Aerodyamics - Theodorsen's function Lift and Moment coefficients
% k = 1; (above)
C = besselh(1,2,k)./(besselh(1,2,k)+i*besselh(0,2,k)); % F+i*G
M_inf = 0;
M_corr = 1/sqrt(1+M_inf^2);

L_alpha = M_corr*(0.5 - i*(1/k)*(1+2*C)-2*(1/k)^2*C);
L_h = M_corr*(1-2*i*(1/k)*C);
M_alpha = M_corr*(3/8 - i*(1/k));
M_h = M_corr*(1/2);

Q_w = zeros(N,N); Q_wo = zeros(N,N); Q_ow = zeros(N,N); Q_o = zeros(N,N);
for m = 1:N
    for n = 1:m
          q_integrand = inline(strrep(strrep(strrep(char(b^2*L_h*PsiB(m)*PsiB(n)),...
              '*','.*'),'/','./'),'^','.^'));
          q_stores = inline(strrep(strrep(strrep(char(PsiB(m)*PsiB(n)),...
              '*','.*'),'/','./'),'^','.^'));
        Q_w(m,n) = pi*rho*(quadl(q_integrand,0,1) + (swidth/L)*sum(L_h*bstores.^2.*q_stores(ystores)));
        % Note, +omega^2 multiplies each term, but the EVP below takes care of
        % this.
        
          q_integrand = inline(strrep(strrep(strrep(char(b^3*(L_alpha-L_h*(0.5+a))*PsiB(m)*PsiT(n)),...
              '*','.*'),'/','./'),'^','.^'));
          q_stores = inline(strrep(strrep(strrep(char(PsiB(m)*PsiT(n)),...
              '*','.*'),'/','./'),'^','.^'));
        Q_wo(m,n) = -pi*rho*(quadl(q_integrand,0,1) + ...
            (swidth/L)*sum(bstores.^3.*(L_alpha-L_h*(0.5+astores)).*q_stores(ystores)));
        
          q_integrand = inline(strrep(strrep(strrep(char(b^3*(M_h-L_h*(0.5+a))*PsiT(m)*PsiB(n)),...
              '*','.*'),'/','./'),'^','.^'));
          q_stores = inline(strrep(strrep(strrep(char(PsiT(m)*PsiB(n)),...
              '*','.*'),'/','./'),'^','.^'));
        Q_ow(m,n) = -pi*rho*(quadl(q_integrand,0,1) + ...
            (swidth/L)*sum(bstores.^3.*(M_h-L_h*(0.5+astores)).*q_stores(ystores)));
        
          q_integrand = inline(strrep(strrep(strrep(char(b^4*(M_alpha-(M_h+L_alpha)*(0.5+a)+L_h*(0.5+a)^2)...
              *PsiT(m)*PsiT(n)),'*','.*'),'/','./'),'^','.^'));
          q_stores = inline(strrep(strrep(strrep(char(PsiT(m)*PsiT(n)),...
              '*','.*'),'/','./'),'^','.^'));
        Q_o(m,n) = pi*rho*(quadl(q_integrand,0,1) + ...
            (swidth/L)*sum(bstores.^4.*(M_alpha-(M_h+L_alpha)*(0.5+astores)+L_h*(0.5+astores).^2).*q_stores(ystores)));
        
        if n ~= m;
            Q_w(n,m) = Q_w(m,n);
            Q_wo(n,m) = Q_wo(m,n);
            Q_ow(n,m) = Q_ow(m,n);
            Q_o(n,m) = Q_o(m,n);
        end
    end
end

% Assemble
Q = [Q_w, Q_wo;
    Q_ow, Q_o];

[Phi,Lam] = eig(K,(M+Q));
% V-g method
    lam(:,ki) = diag(Lam);
    omega(:,ki) = sqrt(abs(lam(:,ki))); % in rad/s
    zts(:,ki) = imag(sqrt(lam(:,ki)))./sqrt(abs(lam(:,ki)));
    U(:,ki) = omega(:,ki)*b/k;
end


% Plots

% Find Flutter Point:
    % assume that flutter happens in 2nd mode, below 600 fps
search_ind = find(U(2,:) < 600);
fl_mode_ind = 2;
fl_U_ind = find(zts(fl_mode_ind,:) < -0.03,1,'last')
U_flutter = U(fl_mode_ind,fl_U_ind)
omega_flutter = omega(fl_mode_ind,fl_U_ind)

% Eliminate jumps in v-g diagram - for one mode only
omega_c = omega; zts_c = zts; U_c = U;
for k = size(omega_c,2):-1:2;
    if abs(omega_c(1,k) - omega_c(1,k-1)) > 5 && abs(omega_c(2,k) - omega_c(1,k-1)) <= 5
        omega_c([1,2],k-1:-1:1) = omega_c([2,1],k-1:-1:1);
        zts_c([1,2],k-1:-1:1) = zts_c([2,1],k-1:-1:1);
        U_c([1,2],k-1:-1:1) = U_c([2,1],k-1:-1:1);
        k
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
sdaxsetb([0,800])
set(get(gca,'Children'),'LineWidth',2)
% ylim([-0.5,1])

