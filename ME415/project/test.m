clear
clc
close all

d = readmatrix('jou.dat');

span = 1.524;

chord = linspace(1e-2,1,100);
% chord = 0.1;

density = 160.185;

for i = 1:length(chord)
    d2 = d .* chord(i);
    area(i) = polyarea(d2(:,1),d2(:,2));
    volume(i) = area(i) * span;
    mass(i) = volume(i) * density;
end

hold on
plot(chord,area)

p = polyfit(chord, area, 3);

plot(chord, p(4) + p(3).*chord + p(2).*chord.^2 + p(1).*chord.^3,'black')

r = corrcoef(chord,area)

% set(gca,'yscale','log')

% plot(d(:,1),d(:,2),'.-')
% axis equal

