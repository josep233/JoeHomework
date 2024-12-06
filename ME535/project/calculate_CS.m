function [cx,cy,Ix,Iy,J,Area] = calculate_CS(naca,chord)

% clear
% clc
% close all
% 
% naca = load('naca0009.txt');
% chord = 1;

% Generate the airfoil coordinates (requires naca4gen function)
naca = chord .* naca;

x = naca(:,1);
y = naca(:,2);

% Find the centroid of the airfoil
cx = mean(x);
cy = mean(y);

% Shift coordinates to the centroidal axis
x_c = x - cx;
y_c = y - cy;

% Initialize variables for area and moments of inertia
Area = 0;  % Area
Ix = 0; % Moment of inertia about x-axis
Iy = 0; % Moment of inertia about y-axis
J = 0;  % Polar moment of inertia

for i = 1:length(x)-1
    % Signed area of the trapezoid
    dA = 0.5 * (x(i) * y(i+1) - x(i+1) * y(i));
    Area = Area + dA; % Accumulate signed area
    
    % Contribution to moments of inertia
    Ix = Ix + abs(dA) * y_c(i)^2;
    Iy = Iy + abs(dA) * x_c(i)^2;
    J = J + abs(dA) * (x_c(i)^2 + y_c(i)^2);
end

% Final area correction
Area = abs(Area); % Ensure positive total area

% Display results
disp(['Area = ', num2str(Area), ' m^2']);
disp(['Moment of inertia about x-axis I_x = ', num2str(Ix), ' m^4']);
disp(['Moment of inertia about y-axis I_y = ', num2str(Iy), ' m^4']);
disp(['Polar moment of inertia J = ', num2str(J), ' m^4']);
