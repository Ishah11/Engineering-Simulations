%The objective of this experiment is to measure the pressure distribution for two airfoils, at various angles of attack and study the changes in lift and drag. 
% NACA 0012 Code
% Load the data and extract relevant columns
D = importdata('Lab2_FriA730_0012.txt');

WT = D.data(:,1); % Wind tunnel setting
T = D.data(:,2); % Air temperature C
DP = D.data(:,3); % Dynamic Pressure (Pa)
AA = D.data(:,4); % Angle of Attack (degrees)
SR = D.data(:,5)*248.84; % Scannivalve Reading(in H20) converted to Pa
SP = D.data(:,6); % Scannivalve Port number
V = D.data(:,7); % Velocity (m/s)
AP = D.data(:,8); % Ambient Pressure (kPa)

% Define parameters
alpha = [0, 4, 8, 10, 12]; % Angles of attack
chord = 100; % Chord length
xd = [0, 2, 9, 19, 30, 40, 50, 60, 70, 80, 100]; % x distance along the chord
x = xd / chord; % Normalized chord length
z = (0.12/0.2)*(0.296*(x).^(0.5) - 0.126*(x) - 0.3516*(x).^2 + 0.2843*(x).^3 - 0.1015*(x).^4); % Camber line

% Error calculations for Scannivalve and Pressure Transducer
Scanni_err = 0.005 * 578.351151715723; % Error of scannivalve in Pa
PresTrans_err = 0.01 * 1365.79987240739; % Error of pressure transducer in Pa

% Total Error
Error = sqrt((DP.^(-1) * Scanni_err).^2 + (-(SR ./ (DP.^2)) * PresTrans_err).^2);

% Pressure Coefficient
Cp = SR ./ DP;

% Initialize matrices for upper and lower Cp and errors
numPoints = length(x); % Number of chord points
CpL = zeros(numPoints, 5); % Lower Cp
CpU = zeros(numPoints, 5); % Upper Cp

ErrorLowermin = zeros(numPoints, 5); % Error bounds for lower Cp
ErrorLowermax = zeros(numPoints, 5);
ErrorUppermin = zeros(numPoints, 5); % Error bounds for upper Cp
ErrorUppermax = zeros(numPoints, 5);

CpU(1, :) = 1; % Set upper Cp initial condition
CpL(1, :) = 1; % Set lower Cp initial condition

% Extract Cp values correctly with matching dimensions (9 Cp values for each surface)
for n = 1:5
   % Each angle of attack corresponds to 9 rows of Cp values.
   CpU(2:10, n) = Cp(((n-1)*9+1):(n*9));  % Corrected index for upper surface
   ErrorUppermax(2:10, n) = CpU(2:10, n) + Error(((n-1)*9+1):(n*9));
   ErrorUppermin(2:10, n) = CpU(2:10, n) - Error(((n-1)*9+1):(n*9));
end

% For the lower surface, similarly:
for n = 6:9
   CpL(2:10, n-4) = Cp(((n-1)*9+1):(n*9)); 
   ErrorLowermin(2:10, n-4) = CpL(2:10, n-4) - Error(((n-1)*9+1):(n*9));
   ErrorLowermax(2:10, n-4) = CpL(2:10, n-4) + Error(((n-1)*9+1):(n*9));
end

% Plotting loop for each angle of attack
for n = 1:5
   a = alpha(n);  % Define the current angle of attack
   figure;
  
     % Plot Cp values for lower surface with solid line
   h1 = plot(x, CpL(:,n), 'b-', 'DisplayName', 'Cp Lower Surface'); 
   hold on;
  
   % Add dots for lower surface data points
   plot(x, CpL(:,n), 'o', 'MarkerFaceColor', 'b'); % Blue dots for lower surface
  
   % Plot Cp values for upper surface with dashed line
   h2 = plot(x, CpU(:,n), 'b--', 'DisplayName', 'Cp Upper Surface');
  
   % Add dots for upper surface data points
   plot(x, CpU(:,n), 'o', 'MarkerFaceColor', 'r'); % Red dots for upper surface
  
   % Plot error lines
   h3 = plot(x, ErrorUppermin(:,n), 'r--', 'DisplayName', 'Error Upper Min');
   h4 = plot(x, ErrorUppermax(:,n), 'r', 'DisplayName', 'Error Upper Max');
   h5 = plot(x, ErrorLowermin(:,n), 'k--', 'DisplayName', 'Error Lower Min');
   h6 = plot(x, ErrorLowermax(:,n), 'k', 'DisplayName', 'Error Lower Max');
  
   % Hold off plotting for the current figure
   hold off;
  
   % Set plot attributes
   grid on;
   axis ij;
   xlabel('Dimensionless Chord');
   ylabel('Pressure Coefficient');
   title(['Pressure Distribution of NACA 0012 Airfoil at \alpha = ' num2str(a)]);
  
   % Create the legend to show which lines are errors
   legend([h1, h2, h3, h4, h5, h6], 'Location', 'best');
  
   % Set legend location and grid
   legend show;
   legend('boxoff');
end


% Initialize angles of attack (in radians)
angles = deg2rad(alpha);
numAngles = length(angles);
cl = zeros(1, numAngles);
cd = zeros(1, numAngles);

% Loop to calculate lift and drag coefficients using predefined variables
for i = 1:numAngles
   % For upper surface
   x = x'; % x-coordinates for upper surface
   z = z'; % z-coordinates for upper surface
   c = CpU(((i-1)*9+1):(i*9)); % Corresponding Cp for upper surface

   % Compute cfx and cfz using the trapezoidal rule for upper surface
   cfxU = trapz(z, c);       % Integral in z-direction
   cfzU = trapz(x, -c);      % Integral in x-direction

   % Calculate lift and drag coefficients (cl and cd) for upper surface
   cl(i) = cfzU * cos(angles(i)) - cfxU * sin(angles(i));
   cd(i) = cfzU * sin(angles(i)) + cfxU * cos(angles(i));

   % For lower surface
   x = x'; % x-coordinates for lower surface
   z = z'; % z-coordinates for lower surface
   c = CpL(((i-1)*9+1):(i*9)); % Corresponding Cp for lower surface

   % Compute cfx and cfz using the trapezoidal rule for lower surface
   cfxL = trapz(z, c);       % Integral in z-direction
   cfzL = trapz(x, -c);      % Integral in x-direction

   % Calculate lift and drag coefficients (cl and cd) for lower surface
   cl(i) = cl(i) + (cfzL * cos(angles(i)) - cfxL * sin(angles(i)));
   cd(i) = cd(i) + (cfzL * sin(angles(i)) + cfxL * cos(angles(i)));
end

% Calculate Lift/Drag ratio
L_D = cl ./ cd;

% Plotting results
figure;
hold on;

% Loop through angles to plot error bars for lower surface
for n = 1:numAngles
   % Plot error bars for lower surface
   errorbar(x, CpL(:, n), ErrorLowermin(:, n), ErrorLowermax(:, n), ...
       'b', 'DisplayName', 'Lower Surface', 'LineStyle', '-', 'Marker', 'o');
  
   % Plot error bars for upper surface
   errorbar(x, CpU(:, n), ErrorUppermin(:, n), ErrorUppermax(:, n), ...
       'r', 'DisplayName', 'Upper Surface', 'LineStyle', '-', 'Marker', 'o');
  
   % Calculate lift, drag and L/D for the current angle
   lD = L_D(n) * ones(size(x)); % Assuming L/D is a constant value for the angle

   % Create a second y-axis for L/D
   yyaxis right;

   % Plot L/D with error bars
   errorbar(x, lD, 0.01 * ones(size(x)), ...
       'g', 'DisplayName', 'L/D', 'LineStyle', '-', 'Marker', 'o'); % Adjust error value as needed
end

% Setting the labels for y-axes
yyaxis left;
ylabel('C_L, C_D');
yyaxis right;
ylabel('L/D');

% Other plot settings
grid on;
axis ij;
xlabel('Dimensionless Chord');
title(['Pressure Distribution of NACA 0012 Airfoil']);
legend('Location', 'Best');

hold off;

% NACA 1422 Code
% Import Data
DUpp = importdata('Lab2_FriA730_4412U.txt');
DLow = importdata('Lab2_FriA730_4412L.txt');

% Extract Data
DPUpp = DUpp.data(:, 3); % Dynamic Pressure (Pa)
AAUpp = DUpp.data(:, 4); % Angle of Attack (degrees)
SRUpp = DUpp.data(:, 5) * 248.84; % Scannivalve Reading (in H20) converted to Pa
SPUpp = DUpp.data(:, 6); % Scannivalve Port number

DPLow = DLow.data(:, 3); % Dynamic Pressure (Pa)
AALow = DLow.data(:, 4); % Angle of Attack (degrees)
SRLow = DLow.data(:, 5) * 248.84; % Scannivalve Reading (in H20) converted to Pa
SPLow = DLow.data(:, 6); % Scannivalve Port number

% Define angles of attack and x-coordinates
alpha44r = [0, 4, 8, 10, 14];
x4412U = [0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.61, 0.81, 1];
x4412L = [0, 0.1, 0.2, 0.3, 0.4, 0.55, 0.7, 1];
z4412U = [0, 0.0449, 0.0643, 0.0874, 0.0975, 0.0980, 0.0797, 0.0458, 0];
z4412L = [0, -0.0293, -0.0274, -0.0225, -0.0180, -0.0119, -0.0063, 0];

% Error calculations for Scannivalve and Pressure Transducer
Scanni_err = 0.005 * 578.351151715723; % Error of scannivalve in Pa
PresTrans_err = 0.01 * 1365.79987240739; % Error of pressure transducer in Pa

ErrorU = sqrt((DPUpp.^(-1) * Scanni_err).^2 + (-(SRUpp ./ (DPUpp.^2)) * PresTrans_err).^2);
ErrorL = sqrt((DPLow.^(-1) * Scanni_err).^2 + (-(SRLow ./ (DPLow.^2)) * PresTrans_err).^2);

% Calculate Cp
CpU = SRUpp ./ DPUpp; % Calculate Cp for upper surface
CpL = SRLow ./ DPLow; % Calculate Cp for lower surface

% Initialize arrays
CpU44 = zeros(length(x4412U), 5);
CpL44 = zeros(length(x4412L), 5);
ErrorLowerMin = zeros(length(x4412L), 5); % Error bounds for lower Cp
ErrorLowerMax = zeros(length(x4412L), 5);
ErrorUpperMin = zeros(length(x4412U), 5); % Error bounds for upper Cp
ErrorUpperMax = zeros(length(x4412U), 5);

% Calculate Cp and error bounds for upper surface
for n = 1:5
   % Extracting corresponding Cp values
   CpU44(2:9, n) = CpU(((n-1)*8+1):(n*8));
end

% Calculate Cp and error bounds for lower surface
for n = 1:5
   % Extracting corresponding Cp values
   CpL44(1:7, n) = CpL(((n-1)*7+1):(n*7));
end

% Ensure the first row for both surfaces is set to 1
CpU44(1, :) = 1;
CpL44(1, :) = 1;

% Plotting results
for n = 1:5
   a = alpha44r(n);
   figure
   hold on
   plot(x4412L, CpL44(:, n), 'b', 'DisplayName', 'Lower Surface'); % Lower surface
   plot(x4412U, CpU44(:, n), 'm', 'DisplayName', 'Upper Surface'); % Upper surface
  
   % Adding error bounds to plots
   %plot(x4412L, ErrorLowerMin(:, n), 'b--', 'DisplayName', 'Lower Min Error');
   %plot(x4412L, ErrorLowerMax(:, n), 'b:', 'DisplayName', 'Lower Max Error');
   %plot(x4412U, ErrorUpperMin(:, n), 'm--', 'DisplayName', 'Upper Min Error');
   %plot(x4412U, ErrorUpperMax(:, n), 'm:', 'DisplayName', 'Upper Max Error');
  
   hold off
   grid on
   axis ij
   xlabel('Dimensionless Chord');
   ylabel('Pressure Coefficient');
   title(['Pressure Distribution of NACA 4412 Airfoil at \alpha = ' num2str(a) '\circ']);
   legend;
end
