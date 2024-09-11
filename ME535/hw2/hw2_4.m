clear
clc
close all

A_hat = 1.2 * (0.8315 - 0.556 * 1i);

omega = (2 * pi) / 16;

t = linspace(0,20,100);

v = real(A_hat * exp(1i * omega * t));

v_dot = real(A_hat * omega * 1i * exp(1i * omega * t));

hold on
plot(t,v,'blue')
plot(t,v_dot,'red')
yline(0)
legend(["${x}$","${\dot{x}}$"],'Interpreter','latex','FontSize',12,'location','southwest')
xlabel('time (ms)')
ylabel('voltage (V)')