clear
clc
close all

k = 392.4;
m2 = 2;

func = @(m) (sqrt(k/m) - 10 * pi) - sqrt(k/(m+m2));

a = fzero(func,1)

