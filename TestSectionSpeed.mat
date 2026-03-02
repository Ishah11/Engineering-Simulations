%This experiment aims to determine the pressure of laboratory air, calibrate a pressure transducer and scannivalve, and determine the test section speed as a function of the fan setting. 

%%Part 1

R = 287; % Universal Gas Constant J/kg.K
T = 32.8 + 273.15; % Temperature from Data Chart in K
P = 96.88 * 10^3; % Pressure from Data chart(Pa)
T_err = 1; % Temperature error in C
P_err = 0.1*10^3; % Pa
rho = P/(R*T); % kg/m^3
rho_err = sqrt((P_err/(R*T))^2 + (P*T_err/(R*T^2))^2);
disp(rho);
disp(rho_err);



%% Lab Part 2
clear;
format long;
% Import Data
D = importdata('Lab1_FriA730_2.txt');
WT = D.data(:,1); % Wind tunnel setting
T = D.data(:,2); % Air temperature C
M = D.data(:,3); % Maonometer(in. of H20)
PT = D.data(:,4); % Pressure Transducer(V)
SR = D.data(:,5); % Scannivalve Reading(V)
V = D.data(:,7); % Velocity (m/s)
AP = D.data(:,8); % Ambient Pressure (kPa)
% Errors from error table
AbsP_err = .1; % Error of absolute pressure reported in kPa
Temp_err = 1; % Error of thermometerreported in C
Mano_bias = -0.062; % Bias of manometer reported in of h2o
Scanni_err = 0.005; % Error of scannivalve reported in V
PresTrans_err = 0.01; % Error of the pressure transducer reported in V
conversion_factor = 249.0889;
% Convert Manometer readings to Pascals
M = M + Mano_bias;
Pa_readings = M * conversion_factor;
SR = SR*-1;
% Error for Manometer (Scannivalve)
p = polyfit(SR,Pa_readings,1);
slope = p(1);
Mano_err = slope * Scanni_err;
x_error = Scanni_err * ones(size(SR));
y_error = Mano_err * ones(size(Pa_readings));
% Plot Scannivalve with error bars and linear correlation
figure(1);
errorbar(SR, Pa_readings, y_error, y_error, x_error, x_error, 'o', 'DisplayName', 'Data Points with Error');
hold on;
% Plot linear fit (SR vs Pa_readings)
SR_fit = linspace(min(SR), max(SR), 100); % Generate points for a smoother line
Pa_fit = polyval(p, SR_fit); % Compute corresponding y-values from the fit
plot(SR_fit, Pa_fit, '-', 'DisplayName', 'Linear Fit', 'LineWidth', 1.5);
% Set axis limits to show a better linear correlation
xlim([min(SR) - 0.1, max(SR) + 0.1]); % Adjust as needed
ylim([min(Pa_readings) - 100, max(Pa_readings) + 100]); % Adjust as needed
% Axis formatting for clearer correlation
axis square; % Set equal scaling for x and y axes
grid on;
% Labels and title
title('Manometer Readings vs Scannivalve Readings');
xlabel('Scannivalve (Volts)');
ylabel('Manometer (Pa)');
% Add legend
legend('show', 'Location', 'Best');
hold off;
% Plot Transducer with error bars and linear correlation
% Error for Manometer (Transducer)
p = polyfit(PT,Pa_readings,1);
slope = p(1);
Mano_err = slope * PresTrans_err;
x_error = PresTrans_err * ones(size(PT));
y_error = Mano_err * ones(size(Pa_readings));
figure(2);
errorbar(PT, Pa_readings, y_error, y_error, x_error, x_error, 'o', 'DisplayName', 'Data Points with Error');
hold on;
% Plot linear fit (SR vs Pa_readings)
PT_fit = linspace(min(PT), max(PT), 100); % Generate points for a smoother line
Pa_fit = polyval(p, PT_fit); % Compute corresponding y-values from the fit
plot(PT_fit, Pa_fit, '-', 'DisplayName', 'Linear Fit', 'LineWidth', 1.5);
% Set axis limits to show a better linear correlation
xlim([min(PT) - 0.1, max(PT) + 0.1]); % Adjust as needed
ylim([min(Pa_readings) - 100, max(Pa_readings) + 100]); % Adjust as needed
% Axis formatting for clearer correlation
axis square; % Set equal scaling for x and y axes
grid on;
% Labels and title
title('Manometer Readings vs Pressure Transducer');
xlabel('Pressure Transducer (Volts)');
ylabel('Manometer (Pa)');
% Add legend
legend('show', 'Location', 'Best');
hold off;



%% Lab Part 3
clear;
format long;
% Import Data
D = importdata('Lab1_FriA730_3.txt');
WT1 = D.data(:,1); % Wind tunnel setting
T = D.data(:,2); % Air temperature C
M1 = D.data(:,3); % Maonometer(in. of H20)
PT1 = D.data(:,4); % Pressure Transducer(V)
SR1 = D.data(:,5); % Scannivalve Reading(V)
SP = D.data(:,6); % Scannivalve port (#)
V = D.data(:,7); % Velocity (m/s)
AP = D.data(:,8); % Ambient Pressure (kPa)
%Average the values from each port so there is only one value
%coresponding to wind speed-------------------------------------------------
% Unique Wind Tunnel Settings
WT = unique(WT1);
numSettings = length(WT);
% Initialize array to hold averaged manometer values
M = zeros(length(WT), 1);
% Loop through each unique wind tunnel setting (Manometer)
for i = 1:length(WT)
   startIndex = (i-1) * 3 + 1;
   endIndex = i * 3;
   M(i) = mean(M1(startIndex:endIndex));
end
% Initialize array to hold averaged Scannivalve values
SR = zeros(length(WT), 1);
% Loop through each unique wind tunnel setting (Scannivalve)
for i = 1:length(WT)
   startIndex = (i-1) * 3 + 1;
   endIndex = i * 3;
   SR(i) = mean(SR1(startIndex:endIndex));
end
% Initialize array to hold averaged Pressure Transducer values
PT = zeros(length(WT), 1);
% Loop through each unique wind tunnel setting (Pressure Transducer)
for i = 1:length(WT)
   startIndex = (i-1) * 3 + 1;
   endIndex = i * 3;
   PT(i) = mean(PT1(startIndex:endIndex));
end
% Initialize array to hold averaged Pressure Transducer values
V_avg = zeros(length(WT), 1);
% Loop through each unique wind tunnel setting (True Velocity)
for i = 1:length(WT)
   startIndex = (i-1) * 3 + 1;
   endIndex = i * 3;
   V_avg(i) = mean(V(startIndex:endIndex));
end
% Errors from error table--------------------------------------------------
T_err = 1; % Error of thermometer reported in C
Mano_bias = -0.062; % Bias of manometer reported in of h2o
SR_err = 0.005; % Error of scannivalve reported in V
PT_err = 0.01; % Error of the pressure transducer reported in V
% Air Density Calculation--------------------------------------------------
R = 287; % J/kg.K
T = 32.8 + 273.15;
P = 96.49 * 10^3; % Pressure (Pa)
rho = P/(R*T); % kg/m^3
% Convert Manometer readings to Pascals
conversion_factor = 249.0889;
M = M + Mano_bias;
M = M * conversion_factor;
% Linear fits of pressure transducer and scannivalve
PT_fit = polyfit(PT,M,1);
SR_fit = polyfit(SR,M,1);
% Convert SR and PT to pressure difference(PD) using the linear fit--------------------------------------------------------------------------
SR_PD = SR_fit(1) * SR + SR_fit(2); % in Pascals
PT_PD = PT_fit(1) * PT + PT_fit(2); % in Pascals
% Calculate velocity using Bernoulli's equation----------------------------------------------------------------------
SR_V = sqrt(2 * SR_PD / rho);
PT_V = sqrt(2 * PT_PD / rho);
M_V = sqrt(2 * M/(rho*(1 - 6.25^-2))); % 6.25 is the area ratio A1/A2
% Calculate the y-intercept and slope errors of Scannivalve-------------------------------------------------------------------
SR_res = M - (SR_fit(1) * SR + SR_fit(2));
% Calculate the standard deviation of the residuals
SR_sd_res = std(SR_res);
% Calculate the standard error of the slope and y-intercept
SR_n = length(SR); % Number of data points
SR_mean = mean(SR);
SR_var = sum((SR - SR_mean).^2);
SR_m_error = SR_sd_res/sqrt(SR_var * (SR_n - 1));
SR_c_error = SR_sd_res * sqrt(sum(SR.^2)/(SR_n * SR_var));
% Calculate the y-intercept and slope errors of Transducer--------------------------------------------------------------------
PT_res = M - (PT_fit(1) * PT + PT_fit(2));
% Calculate the standard deviation of the residuals
PT_sd_res = std(PT_res);
% Calculate the standard error of the slope and y-intercept
PT_n = length(PT); % Number of data points
PT_mean = mean(PT);
PT_var = sum((PT - PT_mean).^2);
PT_m_error = PT_sd_res/sqrt(PT_var * (PT_n - 1));
PT_c_error = PT_sd_res * sqrt(sum(PT.^2)/(PT_n * PT_var));
% Calculate Velocity Errors---------------------------------------------------
% Loop through each Scannivalve reading and calculate velocity and error
for i = 1:length(SR)
  
   % Error propagation for deltaP
   SR_deltaP_err(i) = sqrt((SR(i) * SR_m_error)^2 + SR_c_error^2 + (SR_fit(1) * SR_err)^2);
  
   % Error propagation for velocity
   SR_V_err(i) = (1 / (2 * SR_V(i))) * SR_deltaP_err(i);
end
% Loop through each Pressure Transducer reading and calculate velocity and error
for i = 1:length(PT)
  
   % Error propagation for deltaP
   PT_deltaP_err = sqrt((PT(i) * PT_m_error)^2 + PT_c_error^2 + (PT_fit(1) * PT_err)^2);
   % Error propagation for velocity
   PT_V_err(i) = (1 / (2 * SR_V(i))) * PT_deltaP_err;
end
% Loop through each Manometer reading and calculate velocity and error
for i = 1:length(M) 
   % Error propagation for velocity
   M_V_err(i) = (1 / sqrt(2 * M(i) * rho * (1 - (1/(6.25^2))))) * SR_fit(1) * SR_err;
end
disp(M_V_err)
% Plot SR_V vs WT with error bars------------------------------------------
figure(1);
errorbar(WT, SR_V, SR_V_err, 'o', 'DisplayName', 'Scannivalve Error Bars', 'MarkerSize', 8, 'LineWidth', 1.5);
hold on;
%disp(SR_V_err);
% Fit a linear model
p_SR = polyfit(WT, SR_V, 1);
SR_V_fit = polyval(p_SR, WT);
plot(WT, SR_V_fit, '-k', 'DisplayName', 'Line of Best Fit', 'LineWidth', 1.5);
% Add labels, title, and legend
xlabel('Wind Tunnel Setting (Hz)');
ylabel('Calculated Velocity from Scannivalve (m/s)');
title('Calculated Velocity from Scannivalve vs Wind Tunnel Setting');
legend('show');
grid on;
hold off;
% Plot PT_V vs WT with error bars------------------------------------------
figure(2);
errorbar(WT, PT_V, PT_V_err, 'o', 'DisplayName', 'Manometer Error Bars', 'MarkerSize', 8, 'LineWidth', 1.5);
hold on;
%disp(PT_V_err);
% Fit a linear model
p_PT = polyfit(WT, PT_V, 1);
PT_V_fit = polyval(p_PT, WT);
plot(WT, PT_V_fit, '-b', 'DisplayName', 'Line of Best Fit', 'LineWidth', 1.5);
% Add labels, title, and legend
xlabel('Wind Tunnel Setting (Hz)');
ylabel('Velocity from Pressure Transducer (m/s)');
title('Calculated Velocity from Pressure Transducer vs Wind Tunnel Setting');
legend('show');
grid on;
hold off;
% Plot M_V vs WT with error bars------------------------------------------
figure(3);
errorbar(WT, M_V, M_V_err, 'o', 'DisplayName', 'Manometer Error Bars', 'MarkerSize', 8, 'LineWidth', 1.5);
hold on;
% Fit a linear model
p_M = polyfit(WT, M_V, 1);
M_V_fit = polyval(p_M, WT);
plot(WT, M_V_fit, '-r', 'DisplayName', 'Line of Best Fit', 'LineWidth', 1.5);
% Add labels, title, and legend
xlabel('Wind Tunnel Setting (Hz)');
ylabel('Calculated Velocity from Manometer (m/s)');
title('Calculated Velocity from Manometer vs Wind Tunnel Setting');
legend('show');
grid on;
hold off;
% Plot all velocities on a graph------------------------------------------
figure(4);
plot(WT, M_V_fit, '-r', 'DisplayName', 'Calculated Manometer Velocity', 'LineWidth', 1.5);
hold on;
plot(WT, PT_V_fit, '-b', 'DisplayName', 'Calculated Pressure Transducer Velocity', 'LineWidth', 1.5);
plot(WT, SR_V_fit, ':b', 'DisplayName', 'Calculated Scannivalve Velocity', 'LineWidth', 1.5);
plot(WT,V_avg, '-k', 'DisplayName', 'True Velocity', 'LineWidth', 1.5);
xlabel('Wind Tunnel Setting (Hz)');
ylabel('Velocities (m/s)');
title('Velocities vs Wind Tunnel Setting');
legend('show');
grid on;
hold off;
