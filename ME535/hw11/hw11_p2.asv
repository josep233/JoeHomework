clear
clc
close all

m = 6.47;
c = 0.36;
S = 1;
I_alpha = 0.0699;
k_alpha = 97.7;
ac = 0.227 * c;
e = 0.374 * c;
cg = 0.4 * c;
dCL_dalpha = 4.07;
rho = 1.2;
b = c/2;
k_h = 755;
v = linspace(0,25.5,100);
subtract = zeros(2,2);
x_alpha = cg / b;
S_alpha = m * x_alpha * b;

for i = 1:length(v)
    q = (1/2) * rho * v(i)^2;
    L = q * S * dCL_dalpha;
    My = q * S * e * dCL_dalpha;
    subtract(:,2) = [-L;My];
    
    M = [m, S_alpha; S_alpha, I_alpha];
    K = [k_h, 0; 0, k_alpha];
    Kp = K - subtract;
    
    [Phic,lamss]=polyeig(Kp,0*M,M);
        lkp=lamss(imag(lamss)>0); % keep only root with positive frequency
    lams(:,i)=sort(lkp); 
end

% plot(v,imag(lams)/2/pi)





