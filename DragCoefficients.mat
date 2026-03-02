%The objective of this experiment is to measure the drag coefficients for two airfoils at zero angle of attack. 
% Part 1
% Load the data and extract relevant columns
D = importdata('Lab3_FriA730_Cal.txt');
T = D.data(:,1); % Air temperature (C)
DP = D.data(:,2); % Dynamic Pressure (Pa)
HW = D.data(:,3); % Hotwire Output (V)
V = D.data(:,4); % Velocity (m/s)
loglog(HW,V);
xlabel("log(Hotwire Anemometer)");
ylabel("log(Velocity)");
title("log(Hotwire Anemometer) vs log(Velocity)");
logHW = log(HW); % log(Hotwire Output)
logV = log(V);   % log(Velocity)
% Fit a linear model to the log-log data
coefficients = polyfit(logHW, logV, 1);
% Extract the slope (m) and intercept (b)
m = coefficients(1); % Slope
b = coefficients(2); % Intercept
% Display the slope-intercept form in log-log scale
fprintf('The log-log equation is: log(V) = %.2f * log(HW) + %.2f\n', m, b);
% Find air density rho
Abs_P = 97.02 * 1000;
Temp = 29 + 273.15;
T_err = 1; % Temperature error in C
P_err = 0.1*10^3; % Pa
R = 287.05;
rho = Abs_P/(R*Temp);
rho_err = sqrt((P_err/(R*Temp))^2 + (Abs_P*T_err/(R*Temp^2))^2);

% Part 2
D = importdata('Lab3_FriA730_0012.txt');
HW = D.data(:,3); % Hotwire Output (V)
P = D.data(:,6); % Position (steps)
AV = D.data(:,5); % Air Velocity (m/s)
m = 11.80; %coefficients from the calibration curve
b = -11.03;
%convert hotwire data to m/s
HW = m .* HW + b;
%convert steps to norm
maxP = max(P);
P_norm = P ./ maxP;
% Convert Velocity to norm
maxV = max(HW);
V = HW ./ maxV;
% Plot
figure();
plot(V,P_norm);
xlabel("Hotwire Anemometer (m/s)");
ylabel("Normalized Position");
title("Calculated Velocity vs Position for 0012");
% Convert Steps to meters
P_meter = (P ./ 6400) * (1 ./ 39.37);
% Calculate h1
h1 = trapz(P_meter, HW)/mean(AV);
c = 100/1000;
% Calculate d and Cd
d = (rho * mean(AV.^2) * h1) - (rho * trapz(P_meter, HW.^2));
Cd = d / (mean(PT)*c);

% Part 3
D = importdata('Lab3_FriA730_4412.txt');
HW = D.data(:,3); % Hotwire Output (V)
P = D.data(:,6); % Position (steps)
AV = D.data(:,5); % Air Velocity (m/s)
m = 11.80; %coefficients from the calibration curve
b = -11.03;
%convert hotwire data to m/s
HW = m .* HW + b;
%convert steps to norm
maxP = max(P);
P_norm = P ./ maxP;
% Convert Velocity to norm
maxV = max(HW);
V = HW ./ maxV;
% Plot
figure();
plot(V,P_norm);
xlabel("Hotwire Anemometer (m/s)");
ylabel("Normalized Position");
title("Calculated Velocity vs Position for 4412");
% Convert Steps to meters
P_meter = (P ./ 6400) * (1 ./ 39.37);
% Calculate h1
h1 = trapz(P_meter, HW)/mean(AV);
c = 100/1000;
% Calculate d and Cd
d = (rho * mean(AV.^2) * h1) - (rho * trapz(P_meter, HW.^2));
Cd = d / (mean(PT)*c);
