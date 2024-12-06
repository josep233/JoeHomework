clear% Solve equations of motion for an aeroelastic Goland wing based on strip
% theory for the aerodynamics.
%
% Based on Bisplinghoff, Ashley and Halfman 1996 textbook (orig 1955) solution and
% (?unpublished?) paper by Sanders, Eastep & Huttsell.  Those works were
% based on Goland's paper:
% Goland, M. "The Flutter of a Uniform Cantilever Wing." ASME. J. Appl.
% Mech. December 1945; 12(4): A197–A208. https://doi.org/10.1115/1.4009489  
%
% See also: M. S. Allen and J. A. Camberos, "Comparison of Uncertainty
% Propagation / Response Surface Techniques for Two Aeroelastic Systems,"
% in 50th AIAA Structural Dynamics and Materials Conference Palm Springs,
% CA, 2009.

clear all; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define values for material properties and system geometry
% E:	Young's modulus (Pa)
% rho: Density (kg/m3)
% nu:	Poisson's ratio
% G:	Shear modulus (Pa)
% A:	Cross sectional area (m2)
% I:	Moment of Inertia (m4)
% J:	Polar moment of inertia (m4)
% L1:	Length of beam 1 (m)
% L2:	Length of beam 2 (m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rho_inf = 0.002378; % slugs/ft^3 - air density at sea level

% Light Goland wing (Eastep GolandV5a.nb)
Lw = 20;    % ft - wing length
cw = 6;     % ft - wing chord
bw = cw/2;  % ft - wing half chord - could be variable, see below
aw = -0.34;  % nondim - wing elastic axis location (0.33c?) aw*bw = distance from 
            %   wing centerline to elastic axis
xthw = 0.43;% nondim - distance from wing elastic axis to wing center of mass.
            %   (0.43 c ??) NOT USED - included in Sy
mw = 0.746; % slugs/ft - wing mass per unit length (can be variable)
Sy = 0.447; % slugs-ft/ft - 
Iy = 1.943;% slugs-ft^2/ft - Inertia of wing about elastic axis
E = 1.4976e9; % lb-ft^2 (does he mean lb/ft^2) - modulus
G = 5.616e8;% lb-ft^2 - (does he mean lb/ft^2)  - shear modulus
Iea = 1.58e-2;% ft^4 - area moment of inertia of wing in bending
J = 4.25e-3; % ft^4 - area moment of inertia of wing in torsion

% 
% Zero stores
nstores = 3;
d = 0.2*bw*[1,1,1]; % 0.2 is x_alpha or how far the CG is behind the elastic axis
    % See p. 533 of BAH
mstores = 0*mw*[1,1,1];  % mass of store (assumes 2 ft width)
Istores = 0*Iy + mstores.*d.^2; % Inertia of store
Sstores = 0*mstores.*d;     % slugs-ft/ft - Store static imbalance
astores = 0*[1,1,1];     % ft - a_s*b_s = store CG Location relative to elastic axis
cstores = 0*[1,1,1];    % ft - store chord (see above, length of store in wing chord direction)
bstores = 0*cstores/2;  % ft - store chord (length)
ystores = [20,10,15]./Lw;    % nondim - store position.
swidth = 1; % Store width??
    % Doesn't mesh with paper, there the stores have a length along the wing
    % (span direction).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ritz Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define Psi vectors using Matlab's Symbolic Toolbox
% N:	Series length
% Psi_w:	Basis function vector for displacement
% Psi_o:	Basis function vector for rotation
% For all Psi vectors, index = basis function number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
N=1;
Ndof = 2*N;

syms y real
% Find alpha values for mode functions
for jj = 1:1:N;
    r = (2*jj-1)*pi/2;
    al_v(jj) = fzero(@(y)cos(y)+1/cosh(y),r);
end

for jj=1:N
    R = (sinh(al_v(jj))+sin(al_v(jj)))/(cosh(al_v(jj))+cos(al_v(jj)));
    Psi_w(jj)=sin(al_v(jj)*y)-sinh(al_v(jj)*y)-R*(cos(al_v(jj)*y)-cosh(al_v(jj)*y));
    Psi_o(jj)=sin((2*jj-1)*pi*y/2);
    % Psi_pw(jj)=al_v(jj)*(cos(al_v(jj)*y)-cosh(al_v(jj)*y)+R*(sin(al_v(jj)*y)+sinh(al_v(jj)*y)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Form the "small" partitions of the M,K,C matrices for system components
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M_w = zeros(N,N); K_w = zeros(N,N); M_o = zeros(N,N); K_o = zeros(N,N); M_wo = zeros(N,N);
for m = 1:N
    for n = 1:m
        % Bending
          m_integrand = inline(strrep(strrep(strrep(char(1*Psi_w(m)*Psi_w(n)),'*','.*'),'/','./'),'^','.^'));
            % replace 1* with mw(y)/mw_ref to make wing mass/unit length
            % variable.  Need separate function for stores if this is done.
        M_w(m,n) = mw*quadl(m_integrand,0,1) + sum((mstores/Lw).*m_integrand(ystores));
            % Lw divided out here and in Q
            
          m_integrand = inline(strrep(strrep(strrep(char(1*Psi_w(m)*Psi_o(n)),'*','.*'),'/','./'),'^','.^'));
            % replace 1* with Sy(y)/sy_ref to make wing mass/unit length
            % variable.  Need separate function for stores if this is done.
        M_wo(m,n) = -Sy*quadl(m_integrand,0,1) - sum((Sstores/Lw).*m_integrand(ystores));
            % Lw divided out here and in Q
            
          m_integrand = inline(strrep(strrep(strrep(char(1*Psi_o(m)*Psi_o(n)),'*','.*'),'/','./'),'^','.^'));
            % replace 1* with Sy(y)/sy_ref to make wing mass/unit length
            % variable.  Need separate function for stores if this is done.
        M_o(m,n) = Iy*quadl(m_integrand,0,1) + sum((Istores/Lw).*m_integrand(ystores));
            % Lw divided out here and in Q
        
        % K matrix - stores don't affect stiffness for this model.
          k_integrand = inline(strrep(strrep(strrep(char(1*diff(Psi_w(m),y,2)*diff(Psi_w(n),y,2)),'*','.*'),'/','./'),'^','.^'));
        K_w(m,n) = (E*Iea/Lw^4)*quadl(k_integrand,0,1); 
        
          k_integrand = inline(strrep(strrep(strrep(char(1*diff(Psi_o(m),y,1)*diff(Psi_o(n),y,1)),'*','.*'),'/','./'),'^','.^'));
        K_o(m,n) = (G*J/Lw^2)*quadl(k_integrand,0,1);
        
        if n ~= m;
            M_w(n,m) = M_w(m,n);
            K_w(n,m) = K_w(m,n);
            M_o(n,m) = M_o(m,n);
            K_o(n,m) = K_o(m,n);
            M_wo(n,m) = M_wo(m,n);
        end
    end
end

% Assemble Mass and Stiffness Matrices
M = [M_w, M_wo;
    M_wo.', M_o];
K = [K_w, zeros(N,N);
    zeros(N,N), K_o];

% Compute Frequencies for V = 0
disp('Frequencies for Uncoupled Bending and Torsion: [wb1,wt1]')
[1.8755^2*sqrt(E*Iea/(mw*Lw^4)); (pi/(2*Lw))*sqrt(G*J/Iy)]
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
          q_integrand = inline(strrep(strrep(strrep(char(bw^2*L_h*Psi_w(m)*Psi_w(n)),...
              '*','.*'),'/','./'),'^','.^'));
          q_stores = inline(strrep(strrep(strrep(char(Psi_w(m)*Psi_w(n)),...
              '*','.*'),'/','./'),'^','.^'));
        Q_w(m,n) = pi*rho_inf*(quadl(q_integrand,0,1) + (swidth/Lw)*sum(L_h*bstores.^2.*q_stores(ystores)));
        % Note, +omega^2 multiplies each term, but the EVP below takes care of
        % this.
        
          q_integrand = inline(strrep(strrep(strrep(char(bw^3*(L_alpha-L_h*(0.5+aw))*Psi_w(m)*Psi_o(n)),...
              '*','.*'),'/','./'),'^','.^'));
          q_stores = inline(strrep(strrep(strrep(char(Psi_w(m)*Psi_o(n)),...
              '*','.*'),'/','./'),'^','.^'));
        Q_wo(m,n) = -pi*rho_inf*(quadl(q_integrand,0,1) + ...
            (swidth/Lw)*sum(bstores.^3.*(L_alpha-L_h*(0.5+astores)).*q_stores(ystores)));
        
          q_integrand = inline(strrep(strrep(strrep(char(bw^3*(M_h-L_h*(0.5+aw))*Psi_o(m)*Psi_w(n)),...
              '*','.*'),'/','./'),'^','.^'));
          q_stores = inline(strrep(strrep(strrep(char(Psi_o(m)*Psi_w(n)),...
              '*','.*'),'/','./'),'^','.^'));
        Q_ow(m,n) = -pi*rho_inf*(quadl(q_integrand,0,1) + ...
            (swidth/Lw)*sum(bstores.^3.*(M_h-L_h*(0.5+astores)).*q_stores(ystores)));
        
          q_integrand = inline(strrep(strrep(strrep(char(bw^4*(M_alpha-(M_h+L_alpha)*(0.5+aw)+L_h*(0.5+aw)^2)...
              *Psi_o(m)*Psi_o(n)),'*','.*'),'/','./'),'^','.^'));
          q_stores = inline(strrep(strrep(strrep(char(Psi_o(m)*Psi_o(n)),...
              '*','.*'),'/','./'),'^','.^'));
        Q_o(m,n) = pi*rho_inf*(quadl(q_integrand,0,1) + ...
            (swidth/Lw)*sum(bstores.^4.*(M_alpha-(M_h+L_alpha)*(0.5+astores)+L_h*(0.5+astores).^2).*q_stores(ystores)));
        
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
    U(:,ki) = omega(:,ki)*bw/k;
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
% 
% figure(1)
% subplot(2,1,1);
% plot(U_c.',omega_c.',U(fl_mode_ind,fl_U_ind),omega_c(fl_mode_ind,fl_U_ind),'*'); grid on;
% ylabel('\bfFrequencies (rad/s)');
% set(get(gca,'Children'),'LineWidth',2)
% subplot(2,1,2);
% plot(U_c.',zts_c.',U_c(fl_mode_ind,fl_U_ind),zts_c(fl_mode_ind,fl_U_ind),'*'); grid on;
% xlabel('\bfSpeed (ft/s)'); ylabel('\bfDamping Ratio');
% sdaxsetb([0,800])
% set(get(gca,'Children'),'LineWidth',2)
% % ylim([-0.5,1])



