clear
clc
close all

m = 6.47;
c = 0.36;
S = 1;
I_alpha = 0.0699;
k_alpha = 97.7;
aerodynamic_center = 0.227 * c;
elastic_center = 0.374 * c;
e = abs(aerodynamic_center - elastic_center);
cg = 0.4 * c;
dCL_dalpha = 4.07;
rho = 1.2;
b = c/2;
k_h = 755;
Us = linspace(0,60,1000);
subtract = zeros(2,2);
x_alpha = cg / b;
S_alpha = m * x_alpha * b;

for i = 1:length(Us)
    q = (1/2) * rho * Us(i)^2;
    L = q * S * dCL_dalpha;
    My = q * S * e * dCL_dalpha;
    subtract(:,2) = [-L;My];
    
    M = [m, S_alpha; S_alpha, I_alpha];
    K = [k_h, 0; 0, k_alpha];
    Kp = K - subtract;
    
    [Phi,lam]=eig(K,M);
    fn=sqrt(diag(lam))/2/pi;
    fns(:,i)=sort(fn);

    [Phic,lamss]=polyeig(Kp,0*M,M);
        lkp=lamss(imag(lamss)>0); % keep only root with positive frequency
    % lams(:,i)=sort(lkp);
    lams = 0;
end

figure(1);
subplot(2,1,1)
plot(Us,real(fns),'-',Us,imag(lams)/2/pi,'--'); grid on;
set(get(gca,'Children'),'LineWidth',2)
legend('M1 Clas.','M2 Clas.','M1 SS','M2 SS');
xlabel('Speed (m/s)'); ylabel('Mode Frequencies');
subplot(2,1,2)
plot(Us,-real(1i*fns)./abs(fns),'-',Us,-real(lams)./abs(lams),'--'); grid on;
set(get(gca,'Children'),'LineWidth',2)
xlabel('Speed (m/s)'); ylabel('Damping Ratio');




