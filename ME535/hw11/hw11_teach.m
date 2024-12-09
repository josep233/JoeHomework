clear
clc
close all

%note: a lot of this came from Dr Allen's code.

f0_1 = 4.57;
f0_2 = 5.95;

m = 6.47;
c = 0.36;
b = c/2;
S = 1 * c;
Ialpha = 0.0699;
% Kalpha = 97.7;
aerodynamic_center = 0.227 * c;
elastic_center = 0.374 * c;
e = abs(aerodynamic_center - elastic_center);
xcg = abs(elastic_center - 0.4 * c);
dCL_da = 4.07;
rhoair = 1.2;

% Kh = 755;
Kalpha = (f0_2 * 2 * pi)^2 * Ialpha;
Kh = (f0_1 * 2 * pi)^2 * m;
Us = linspace(0,60,1000);
subtract = zeros(2,2);
x_alpha = xcg / b;
S_alpha = m * x_alpha * b;

actual_airspeed = [22.5:0.25:25.5];
actual_freq1 = [5.9,5.95,5.8];

% Check frequencies of bounce and torsion modes when flow is zero
    K=[Kh,0;
        0, Kalpha];
    M=[m,xcg*m;
        xcg*m, Ialpha];
[Phi0,lam0]=eig(K,M);
disp('Natural Frequencies at U=0');
fn0=sqrt(diag(lam0))/2/pi

% Us=[0:0.1:35]; % Airspeed in m/s
qs=rhoair*Us.^2/2;
fns=zeros(2,length(qs)); lams=fns;

for k=1:length(qs)
    K=[Kh,qs(k)*S*dCL_da;
        0, Kalpha-qs(k)*e*S*dCL_da];
    M=[m,xcg*m;
        xcg*m, Ialpha];
[Phi,lam]=eig(K,M);
fn=sqrt(diag(lam))/2/pi;
fns(:,k)=sort(fn);
% Alternative, solve for state-space eigenvalues instead.
    [Phic,lamss]=polyeig(K,0*M,M);
        lkp=lamss(imag(lamss)>0); % keep only root with positive frequency
    % lams(:,k)=sort(lkp);
    lams = 0;
end

figure(1);
subplot(2,1,1)
hold on
plot(Us,real(fns),'-',Us,imag(lams)/2/pi,'--'); grid on;
set(get(gca,'Children'),'LineWidth',2)
legend('M1 Clas.','M2 Clas.','M1 SS','M2 SS');
xlabel('Speed (m/s)'); ylabel('Mode Frequencies');
subplot(2,1,2)
hold on
plot(Us,-real(1i*fns)./abs(fns),'-',Us,-real(lams)./abs(lams),'--'); grid on;
set(get(gca,'Children'),'LineWidth',2)
xlabel('Speed (m/s)'); ylabel('Damping Ratio');

% Dynamic Aeroelasticity with Theodorsen Aerodynamics
%
b=c/2; % half chord
a=e/b; % nondimensional elastic axis 
% Find range of k to study.
Udes=[1,60]; % m/s
kdes=fn0(1)*2*pi*b./Udes(2:-1:1); % use first mode frequency to find k's
ksD=linspace(kdes(1),kdes(2),1000);
fnsD=zeros(2,length(ksD)); lamsD=fnsD; ztsD=fnsD;
% Structural stiffness matrices are constant
K=[Kh,0;
    0, Kalpha];
M=[m,xcg*m;
    xcg*m, Ialpha];

for kk=1:length(ksD);
    k=ksD(kk);
    C = besselh(1,2,k)./(besselh(1,2,k)+1i*besselh(0,2,k)); % F+i*G
    Lh = -((- k + C*2i))/k;
    Lal = -(b*(2*C + k*1i + C*k*1i + a*k^2 - C*a*k*2i))/k^2;
    Mh = (b*(- a*k + C*1i + C*a*2i))/k; 
    Mal = (b^2*(8*C - k*4i + 16*C*a + C*k*4i + a*k*8i + k^2 - 8*a^2*k^2 - C*a^2*k*16i))/(8*k^2);
    Maero=pi*rhoair*b^2*[Lh,Lal; Mh,Mal];
    
    % Solve for state-space eigenvalues
    [Phic,lamss]=polyeig(K,0*M,M+Maero);
    lkp=lamss(imag(lamss)>0); % keep only root with positive imaginary part
        % Note, this will fail if the roots are not in complex conjugate
        % pairs, i.e. they become real.
    lamsD(:,kk)=sort(lkp); 
    fnsD(:,kk)=abs(lamsD(:,kk))/2/pi;
    ztsD(:,kk)=-real(lamsD(:,kk))./abs(lamsD(:,kk));
end

% Find U values for each mode
Um1=fnsD(1,:)*2*pi*b./ksD;
Um2=fnsD(2,:)*2*pi*b./ksD;

figure(2);
subplot(2,1,1);
hold on
plot(Um1,fnsD,Um2,fnsD); grid on;
legend('M1 Flutter M1','M1 Flutter M2','M2 Flutter M1','M2 Flutter M2');
xlabel('Speed (m/s)'); ylabel('Mode Frequencies');
subplot(2,1,2);
hold on
plot(Um1,ztsD,Um2,ztsD); grid on;
legend('M1 Flutter M1','M1 Flutter M2','M2 Flutter M1','M2 Flutter M2');
xlabel('Speed (m/s)'); ylabel('Damping Ratios');