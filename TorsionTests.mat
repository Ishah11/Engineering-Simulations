%This experiment carried four torsion tests on samples made of 1018 steel, 2024 Aluminum, and 6061-T6 Aluminum, with both solid and hollow for each metal type. First, load-extensometer data was collected for all circular rods to determine the angle of twist for each value of torque. Next, the value of the shear modulus and the torque-angle twist relations were calculated. These were used to find the variations of shear stress and shear strain. Experimental results closely matched theoretical values, with minor discrepancies attributed to material inconsistencies and test setup limitations.

clear; clc;
filenames = {'Th_100_T_HA', 'Th_100_T_SA', 'Th_100_T_HS', 'W_0730_G1_SS.csv'};
materials = {'Hollow Aluminum', 'Solid Aluminum', 'Hollow Steel', 'Solid Steel'};
colors = {'b-', 'r-', 'k--', 'k-'}; % Colors for different materials
for i = 1:length(filenames)
   filename = filenames{i};
   opts = detectImportOptions(filename);
   opts.DataLines = [11, Inf];
   % Import the data
   data = readtable(filename, opts);
   % Extract relevant columns
   Torque = data{:, 3};
   Rel_Angle = data{:, 7};
   % Create a new figure for each material
   figure(i);
   plot(Rel_Angle, Torque, colors{i}, 'LineWidth', 1.5);
  
   % Customize plot
   grid on;
   xlabel('Relative Angle of Twist (deg)');
   ylabel('Torque (Lbs-in)');
   title(['Torque vs. Relative Angle of Twist for ', materials{i}]);
   legend(materials{i}, 'Location', 'best');
end
d_solid_Al = 9.61 / 25.4; % Solid Aluminum diameter
d_solid_Steel = 9.42 / 25.4; % Solid Steel diameter
d_hollow_Al_outer = 13.08 / 25.4; % Hollow Aluminum outer diameter
d_hollow_Steel_outer = 12.75 / 25.4; % Hollow Steel outer diameter
t = 65 / 1000; % Thickness in inches
% Compute inner diameters
d_hollow_Al_inner = d_hollow_Al_outer - 2 * t;
d_hollow_Steel_inner = d_hollow_Steel_outer - 2 * t;
% Compute polar moment of inertia (J)
J_solid_Al = (pi / 32) * d_solid_Al^4;
J_solid_Steel = (pi / 32) * d_solid_Steel^4;
J_hollow_Al = (pi / 32) * (d_hollow_Al_outer^4 - d_hollow_Al_inner^4);
J_hollow_Steel = (pi / 32) * (d_hollow_Steel_outer^4 - d_hollow_Steel_inner^4);
% Define known experimental values
L = 11.25;
J_values = [J_hollow_Al, J_solid_Al, J_hollow_Steel, J_solid_Steel];
E_values = [11217755.28, 11217755.28, 30549791.827, 30549791.827];
% Literature values for comparison
G_lit = [3.7e6, 3.7e6, 11.5e6, 11.5e6]; % Literature shear modulus values (psi)
nu_lit = [0.34, 0.33, 0.33, 0.33]; % Literature Poisson’s ratio
% Storage for computed values
G_estimates = zeros(1, length(filenames)); % Estimated Shear Modulus
nu_estimates = zeros(1, length(filenames)); % Estimated Poisson's Ratio
filenames1 = {'MEE324_Th_100_T_HA - MEE324_Th_100_T_HA.csv', 'MEE324_Th_100_T_SA - MEE324_Th_100_T_SA.csv', 'MEE324_Th_100_T_HS - MEE324_Th_100_T_HS.csv', 'AEE325_W_0730_G1_SS - AEE325_W_0730_G1_SS.csv'};
for i = 1:length(filenames1)
   filename = filenames1{i};
   % Read data starting from the row where the table begins
   opts = detectImportOptions(filename);
   opts.DataLines = [11, Inf]; % Skip metadata and start from actual data
   % Import the data
   data = readtable(filename, opts);
   % Extract relevant columns
   Torque = data{:, 3}; % Column 3: Torque (Lbs-in)
   Rel_Angle = data{:, 7} * (pi/180); % Convert deg to radians
   % Perform linear regression: T = m * theta
   p = polyfit(Rel_Angle, Torque, 1); % Fit first-degree polynomial (linear)
   m = p(1); % Slope of regression line
   % Compute shear modulus G
   G_estimates(i) = (m * L) / J_values(i);
   % Compute Poisson’s ratio
   nu_estimates(i) = (E_values(i) / (2 * G_estimates(i))) - 1;
end
% Tabulate results
Results = table(materials', G_estimates', nu_estimates', G_lit', nu_lit', ...
   'VariableNames', {'Material', 'Estimated_G', 'Estimated_nu', 'Literature_G', 'Literature_nu'});
% Display table
disp('Comparison of Experimental and Literature Values:');
disp(Results);
% Define parameters
r_solid = linspace(0, d_solid_Al/2, 100); % Radial positions for solid shaft
r_hollow = linspace(d_hollow_Al_inner/2, d_hollow_Al_outer/2, 100); % Radial positions for hollow shaft
% Normalize stress values
tau_solid = r_solid / max(r_solid);
tau_hollow = (r_hollow - min(r_hollow)) / (max(r_hollow) - min(r_hollow));
% Plot results
figure(5);
subplot(1,2,1);
plot(r_solid, tau_solid, 'b', 'LineWidth', 1.5);
grid on;
xlabel('Radial Position (in)');
ylabel('Normalized Shear Stress');
title('Shear Stress in Solid Shaft');
subplot(1,2,2);
plot(r_hollow, tau_hollow, 'r', 'LineWidth', 1.5);
grid on;
xlabel('Radial Position (in)');
ylabel('Normalized Shear Stress');
title('Shear Stress in Hollow Shaft');
% Given parameters
r_i = 0.5; % Inner radius (in)
r_o = 1.5; % Outer radius (in)
T = 100; % Torque
J = (pi/2) * (r_o^4 - r_i^4); % Approximate polar moment of inertia
G_i = 80e9; % Inner shear modulus (Pa)
G_o = 30e9; % Outer shear modulus (Pa)
% Radial positions
r = linspace(r_i, r_o, 100);
% Shear stress distribution
tau = (T .* r) / J;
% Shear modulus variation
G = G_i + (G_o - G_i) * ((r - r_i) / (r_o - r_i));
% Shear strain distribution
gamma = tau ./ G;
% Plot
figure(6);
plot(tau, gamma, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('Shear Stress (Pa)');
ylabel('Shear Strain');
title('Shear Stress vs. Shear Strain in a Composite Bar');

