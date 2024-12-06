clear
clc
close all

%cross section of airfoil information, approximated as rectangle. Using
%NACA0009, which means the thickness is 9% of the chord.
chord = 6;
thickness = chord * .09;
L = 20;
% I = (chord * thickness^3) / 12; %this is just bh^3 / 12
I = 1.943;
J = (chord * thickness^3 / 16) * (16 / 3 - 3.36 * thickness / chord * (1 - thickness^4 / (12 * chord^4))); %from wikipedia
A = thickness * chord;


% m = (rho * A * L);
m = 0.746;

rho = 32; %kg/m^3
% E = 6.5e6; %Pa
E = 31.7E6 * m / I;
nu = 0.33;
G = E / (2 * (1 + nu));
J = 1.23E6 * I / G;

aerodynamic_center = 0.010 * chord + chord / 4; %from Figure 10
center_of_gravity = chord * 0.43; %for a symmetric airfoil, the center of gravity is at the center.
elastic_center = chord * 0.33;
e = abs(elastic_center - aerodynamic_center);
x = abs(elastic_center - center_of_gravity);

syms y

beta = 0.6 * pi / L;
Cb = 1;
Ct = 1;

% phiB = Cb * (cosh(beta * y) - 0.734 * sinh(beta * y) - cos(beta * y) + 0.734 * sin(beta * y));
% phiT = Ct * sin((pi * y)/(2 * L));

N = 1;

% Find alpha values for mode functions
for jj = 1:1:N;
    r = (2*jj-1)*pi/2;
    al_v(jj) = fzero(@(y)cos(y)+1/cosh(y),r);
end

for jj=1:N
    R = (sinh(al_v(jj))+sin(al_v(jj)))/(cosh(al_v(jj))+cos(al_v(jj)));
    phiB(jj)=sin(al_v(jj)*y)-sinh(al_v(jj)*y)-R*(cos(al_v(jj)*y)-cosh(al_v(jj)*y));
    phiT(jj)=sin((2*jj-1)*pi*y/2);
    % Psi_pw(jj)=al_v(jj)*(cos(al_v(jj)*y)-cosh(al_v(jj)*y)+R*(sin(al_v(jj)*y)+sinh(al_v(jj)*y)));
end

M = double([int(rho * A * phiB^2,y,0,L), -x * int(rho * A * phiB * phiT,y,0,L); -x * int(rho * A * phiB * phiT,y,0,L), int(rho * A * phiT^2 + rho * J * phiT^2,y,0,L)]);
K = double([int(E * I * diff(phiB,y,2)^2,y,0,L), -x * int(E * I * diff(phiB,y,2) * diff(phiT,y,2),y,0,L); -x * int(E * I * diff(phiB,y,2) * diff(phiT,y,2),y,0,L), int(E * I * diff(phiT,y,2)^2 + G * J * diff(phiT,y,1),y,0,L)]);

omega_b = (3.55 / L^2) * sqrt(E * I / m);
omega_t = (pi / (2 * L)) * sqrt(G * J / I);

% Dynamic Aeroelasticity with Theodorsen Aerodynamics
%
b=chord/2; % half chord
a=e/b; % nondimensional elastic axis 
% Find range of k to study.
fn0 = [omega_b,omega_t] ./ (2*pi*b);
Udes=[1,500]; % m/s
kdes=fn0(1)*2*pi*b./Udes(2:-1:1) % use first mode frequency to find k's
ksD=linspace(kdes(1),kdes(2),100);
fnsD=zeros(2,length(ksD)); lamsD=fnsD; ztsD=fnsD;

for kk=1:length(ksD)
    k=ksD(kk);
    C = besselh(1,2,k)./(besselh(1,2,k)+1i*besselh(0,2,k)); % F+i*G
    Lh = -((- k + C*2i))/k;
    Lal = -(b*(2*C + k*1i + C*k*1i + a*k^2 - C*a*k*2i))/k^2;
    Mh = -(b*(a*k + C*a*2i - C*1i))/k;
    Mal = (b^2*(8*C - k*4i - 16*C*a + C*k*4i + a*k*8i + k^2 - 8*a^2*k^2 - C*a*k*16i + C*a^2*k*16i))/(8*k^2);
    Maero=pi*rho*b^2*[Lh,Lal; Mh,Mal];

    hs = [Lh;Mh];
    as = [Lal,Mal];
    Q = double(-int(phiB,y,0,L)) * hs - double(int(phiT,y,0,L)) * hs + double(int(phiB,y,0,L)) * as;
    
    % Solve for state-space eigenvalues
    [Phic,lamss]=polyeig(K,0*M,M+Q);
    lkp=lamss(imag(lamss)>0); % keep only root with positive imaginary part
        % Note, this will fail if the roots are not in complex conjugate
        % pairs, i.e. they become real.
    lamsD(:,kk)=sort(lkp); 
    fnsD(:,kk)=abs(lamsD(:,kk))/2/pi;
    ztsD(:,kk)=-real(lamsD(:,kk))./abs(lamsD(:,kk));
end

% Find U values for each mode
Um1=(fnsD(1,:)*2*pi*b./ksD)./1.467;
Um2=(fnsD(2,:)*2*pi*b./ksD)./1.467;

figure(2);
subplot(2,1,1);
plot(Um1,fnsD*2*pi,Um2,fnsD*2*pi); grid on;
legend('M1 Flutter M1','M1 Flutter M2','M2 Flutter M1','M2 Flutter M2');
xlabel('Speed (m/s)'); ylabel('Mode Frequencies');
subplot(2,1,2);
plot(Um1,ztsD,Um2,ztsD); grid on;
legend('M1 Flutter M1','M1 Flutter M2','M2 Flutter M1','M2 Flutter M2');
xlabel('Speed (m/s)'); ylabel('Damping Ratios');