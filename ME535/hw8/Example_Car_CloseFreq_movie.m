% Solution of 2DOF Car model for various initial displacements.

roa = 1/2;%
 M = [.25+roa^2 .25-roa^2; 
     .25-roa^2 .25+roa^2];
 K = diag([1    ,1]); % diag([1,1.5]);
 
% % Case with equal natural frequencies, r/a=1/2
% M=[0.5, 0; 0, 0.5];
% K=eye(2);
 
 q_0 = [1; 0];
% q_0=[1; 1];

[phi,lam] = eig(K,M);
wns = sqrt(diag(lam));

 % Sort & Normalize Eigenvectors to unity modal mass and Check Orthogonality
[lam_sort,lam_indx] = sort(diag(lam));
wns = sqrt(lam_sort); % *sqrt(k/m)

phi_sort = (phi(:,lam_indx));
mu = phi_sort.'*M*phi_sort;
PHI = real(phi_sort*sqrt(inv(mu)))

check_orth = norm(PHI.'*M*PHI-eye(size(phi)))

n_0 = PHI'*M*q_0*9.81; % *m/k
% Initial velocities assumed to be zero.
t = [0:0.2:30];
eta_t = [n_0(1)*cos(wns(1)*t);
    n_0(2)*cos(wns(2)*t)];
q = PHI*eta_t;

figure(1)
subplot(2,1,1)
plot(t,eta_t); grid on;
xlabel('time t*(k/m)^0.5'); ylabel('Modal Disp.');
subplot(2,1,2);
plot(t,q); grid on;
xlabel('time t*(k/m)^0.5'); ylabel('Physical Disp.');

figure(2)
for ii = 1:1:length(t);
    plot([-0.5 0.5].', [10 10].','o:',[-0.7 0.7].', [0 0].','k');
    line([-0.5 0.5].',q(:,ii)+10,'LineWidth', 7); grid on;
    xlabel('X-position (*L)'); ylabel('Physical Displacement (m)');
    title(['Time (m/k)^0.5 = ' num2str(t(ii),'%2.1f')])
    axis([-0.7 0.7 -5 25]);
    mov1(ii) = getframe(2);
%     MkGifMovie('CarICResp_M2',ii);
end

movie(mov1,1,10)

% [fname_done] = MkGifMovie('CarICResp_M2',ii,'done');

% movie2avi(mov1,'CarModel.avi','compression', 'indeo3','quality',100)
