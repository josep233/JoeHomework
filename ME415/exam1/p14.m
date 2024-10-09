clear
clc
close all

span = 2;
aspect_ratio = 8;
thickness_to_chord = 0.12;
flight_speed = 20;
kinematic_viscosity = 1e-5;

S = span^2 / aspect_ratio;
mean_geometric_chord = S / span;

Re = kinematic_viscosity^(-1) * flight_speed * mean_geometric_chord;

cf = .074 / Re^0.2;
k = 1 + 2*thickness_to_chord + 100 * thickness_to_chord^4;
planform_area = (span) * mean_geometric_chord;
Sexposed = span * mean_geometric_chord;
Swet = 2 * (1 + 0.2 * thickness_to_chord) * Sexposed;

coefficient_parasitic_drag = k * cf * Swet / planform_area
