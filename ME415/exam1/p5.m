clear
clc
close all

diameter = 10;
length = 70;
mach_number = 0.5;
speed_of_sound = 300;
span = 60;
kinematic_viscosity = 1e-5;
fuselage_wetted_area_div_wing_reference_area = 8;

%inferred values
temperature = (speed_of_sound^2) / (1.4 * 287.053);
dynamic_viscosity = 1.458e-6 * temperature^(3/2) / (temperature + 110.4);
density = dynamic_viscosity / kinematic_viscosity;   
flight_speed = speed_of_sound * mach_number;
% Re = density * flight_speed * length / dynamic_viscosity;
Re = kinematic_viscosity^(-1) * flight_speed * length;

cf_inc = 0.074 / (Re^0.2);
cf = cf_inc * 1;
% cf = (1 + 0.144*mach_number^2)^(-0.65) * cf_inc;
Swet = pi * diameter * length;
dynamic_pressure = (1/2) * density * flight_speed^2;
drag_skin_friction = cf * dynamic_pressure * Swet;
skin_friction_coefficient = drag_skin_friction / (dynamic_pressure * Swet / fuselage_wetted_area_div_wing_reference_area)
skin_friction_coefficient = cf * Swet / (Swet / 8)