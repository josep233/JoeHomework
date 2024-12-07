clear
clc
close all

x_alpha = 0.05;
r_alpha = 0.5;
omegah_omegaalpha = 0.5;
bigratio = 7.5;
e_b = 0.4;
dcl_dalpha = 2 * pi;

Us = linspace(0,1.2,100);

for i = 1:length(Us)

    M = [1, x_alpha; x_alpha, r_alpha^2];
    K = [omegah_omegaalpha^2, 2 * bigratio^(-1) * Us(i)^2; 0, r_alpha^2 - 2 * bigratio^(-1) * e_b * Us(i)^2];

[Phi,lam]=eig(K,M);
fn=sqrt(diag(lam));
fns(:,i)=sort(fn);
% Alternative, solve for state-space eigenvalues instead.
    [Phic,lamss]=polyeig(K,0*M,M);
        lkp=lamss(imag(lamss)>0); % keep only root with positive frequency
    % lams(:,i)=sort(lkp);
    lams = 0;

end

figure(1);
subplot(2,1,1)
hold on
plot(Us,real(fns)/2/pi,'-',Us,imag(lams)/2/pi,'--'); grid on;
set(get(gca,'Children'),'LineWidth',2)
legend('M1 Clas.','M2 Clas.','M1 SS','M2 SS');
xlabel('Speed (m/s)'); ylabel('Mode Frequencies');
subplot(2,1,2)
hold on
plot(Us,-real(1i*fns)./abs(fns),'-',Us,-real(lams)./abs(lams),'--'); grid on;
set(get(gca,'Children'),'LineWidth',2)
xlabel('Speed (m/s)'); ylabel('Damping Ratio');


