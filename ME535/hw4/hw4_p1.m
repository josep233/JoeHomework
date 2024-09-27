clear
clc
close all

F0 = 1;

T = 5;
t = linspace(0,10,100);
f = (-F0./T).*t.*heaviside(t) - (-F0./T).*(t-T).*heaviside(t-T) + F0*heaviside(t-T);

plot(t,f);