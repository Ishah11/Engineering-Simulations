%The objective of this experiment was to measure the low-speed lift and drag characteristics of a cambered finite wing and to observe the effects of adding various winglet configurations. 
clear all;
gravityG1 = {'Lab4_G3_none_grav.csv','Lab4_G3_20deg_grav.csv','Lab4_G3_40deg_grav.csv','Lab4_G3_60deg_grav.csv'};
angleRe150 = {'Lab4_G3_none_150k.csv','Lab4_G3_20deg_150k.csv','Lab4_G3_40deg_150k.csv','Lab4_G3_60deg_150k.csv'};
angleRe300 = {'Lab4_G3_none_300k.csv','Lab4_G3_20deg_300k.csv','Lab4_G3_40deg_300k.csv','Lab4_G3_60deg_300k.csv'};
D = readmatrix('Lab4_G3_Conditions.txt');
P = D(:,3);
P = P * 1000;          % Convert pressure from kPa to Pa
T = D(:,2);       % Extract 2nd column for temperature in Celsius
T = T + 273.15;        % Convert temperature to Kelvin
R = 287.05;            % Specific gas constant for dry air in J/(kg*K)
rho = P ./ (R .* T);    % Calculate density using the ideal gas law
sizeT = length(T);      % Get the length of the temperature array
mu = zeros(1, sizeT);   % Preallocate mu as a numeric array (list) with zeros
for i = 1:sizeT
   mu(i) = 1.716e-5 * ((273.15 + 110.4) / (T(i) + 110.4)) * (T(i) / 273.15)^(3/2);
end
c = 0.1397;
b = 0.1524;
s = c*b;
Re150 = D(:,4);
Re300 = D(:,5);
V150 = Re150.*mu./(rho*c);
V300 = Re300.*mu./(rho*c);
for i = 1:sizeT
   qinf150(i) = 0.5*rho(i)*V150(i)^2;
end
for i = 1:sizeT
   qinf300(i) = 0.5*rho(i)*V300(i)^2;
end
%%Angles, Cl, and Cd with Re 150K
%------------------------------------------------------------------------
% Initialize cell arrays to store Cl and alpha values for each file
Cl_all = cell(1, length(angleRe150));
alpha_all = cell(1, length(angleRe150));
% Loop through each file and calculate Cl and Cd
for i = 1:length(angleRe150)
   D = readmatrix(angleRe150{i},'delimiter', ",");
   G = readmatrix(gravityG1{i},'delimiter', ",");
   Fa = D(:,7)-G(:,7);
   Fn = D(:,8)-G(:,8);
   alpha = D(:,10);
  
   % Calculate Lift (L) and Drag (D) forces
   L = Fn .* cos(deg2rad(alpha)) - Fa .* sin(deg2rad(alpha));
   D = Fn .* sin(deg2rad(alpha)) + Fa .* cos(deg2rad(alpha));
  
   % Calculate Cl and Cd
   Cl = L ./ (qinf150(i).*s);
   Cd = D ./ (qinf150(i).*s);
  
   del_Fn = 0.05;
   del_Fa = 0.05;
   del_alpha = deg2rad(0.05);
   del_qinf = 0.5;
   Cl_Fn = cos(deg2rad(alpha))./(qinf150(i).*s);
   Cl_Fa = -sin(deg2rad(alpha))./(qinf150(i).*s);
   Cl_alpha = -Fn .* sin(deg2rad(alpha)) - Fa .* cos(deg2rad(alpha));
   Cl_alpha = Cl_alpha./(qinf150(i).*s);
   Cl_qinf = Fn .* cos(deg2rad(alpha)) - Fa .* sin(deg2rad(alpha));
   Cl_qinf = Cl_qinf./(qinf150(i)^2.*s);
   Cl_err = Cl .* sqrt((Cl_Fn.*del_Fn).^2+(Cl_Fa.*del_Fa).^2+(Cl_alpha.*del_alpha).^2+(Cl_qinf.*del_qinf).^2);
   Cd_Fn = sin(deg2rad(alpha))./(qinf150(i).*s);
   Cd_Fa = cos(deg2rad(alpha))./(qinf150(i).*s);
   Cd_alpha = Fn .* cos(deg2rad(alpha)) - Fa .* sin(deg2rad(alpha));
   Cd_alpha = Cd_alpha./(qinf150(i).*s);
   Cd_qinf = -Fn .* cos(deg2rad(alpha)) - Fa .* sin(deg2rad(alpha));
   Cd_qinf = Cd_qinf./(qinf150(i)^2.*s);
   Cd_err = Cd .* sqrt((Cd_Fn.*del_Fn).^2+(Cd_Fa.*del_Fa).^2+(Cd_alpha.*del_alpha).^2+(Cd_qinf.*del_qinf).^2);
   % Store the Cl and alpha values for this file
   Cl_err_all{i} = Cl_err;
   Cd_err_all{i} = Cd_err;
   Cd_all{i} = Cd;
   Cl_all{i} = Cl;
   alpha_all{i} = alpha;
end
% Plot all three Cl vs. alpha lines on the same plot
figure(1);
hold on;  % Keep the plot open to overlay each line
for i = 1:length(angleRe150)
   plot(cell2mat(alpha_all(i)), cell2mat(Cl_all(i)), 'DisplayName', sprintf('File %d', i));  % Use legend to identify each line
   errorbar(cell2mat(alpha_all(i)), cell2mat(Cl_all(i)), cell2mat(Cl_err_all(i)));
end
hold off;
% Label the axes and add a legend
xlabel('\alpha (Angle of Attack)');
ylabel('C_L (Lift Coefficient)');
title('C_L vs. \alpha for Re 150K for Winglet Angles');
legend('None','20 degrees','40 degrees','60 degrees');  % Show legend with labels for each line
grid on;
% Plot all three Cl vs. Cd lines on the same plot
figure(2);
hold on;  % Keep the plot open to overlay each line
for i = 1:length(angleRe150)
   plot(cell2mat(Cd_all(i)), cell2mat(Cl_all(i)), 'DisplayName', sprintf('File %d', i));  % Use legend to identify each line
   errorbar(cell2mat(Cd_all(i)), cell2mat(Cl_all(i)), cell2mat(Cl_err_all(i)), 'both');
end
hold off;
% Label the axes and add a legend
xlabel('C_D (Drag Coefficient)');
ylabel('C_L (Lift Coefficient)');
title('C_L vs. C_D for Re 150K for Winglet Angles');
legend('None','20 degrees','40 degrees','60 degrees');  % Show legend with labels for each line
grid on;
%%Angles, Cl, and Cd with Re 300K
%-----------------------------------------------------------------------
% Initialize cell arrays to store Cl and alpha values for each file
Cl_all = cell(1, length(angleRe300));
alpha_all = cell(1, length(angleRe300));
% Loop through each file and calculate Cl and Cd
for i = 1:length(angleRe300)
   D = readmatrix(angleRe300{i},'delimiter', ",");
   G = readmatrix(gravityG1{i},'delimiter', ",");
   Fa = D(:,7)-G(:,7);
   Fn = D(:,8)-G(:,8);
   alpha = D(:,10);
  
   % Calculate Lift (L) and Drag (D) forces
   L = Fn .* cos(deg2rad(alpha)) - Fa .* sin(deg2rad(alpha));
   D = Fn .* sin(deg2rad(alpha)) + Fa .* cos(deg2rad(alpha));
  
   % Calculate Cl and Cd
   Cl = L ./ (qinf300(i).*s);
   Cd = D ./ (qinf300(i).*s);
  
   del_Fn = 0.05;
   del_Fa = 0.05;
   del_alpha = deg2rad(0.05);
   del_qinf = 0.5;
   Cl_Fn = cos(deg2rad(alpha))./(qinf300(i).*s);
   Cl_Fa = -sin(deg2rad(alpha))./(qinf300(i).*s);
   Cl_alpha = -Fn .* sin(deg2rad(alpha)) - Fa .* cos(deg2rad(alpha));
   Cl_alpha = Cl_alpha./(qinf300(i).*s);
   Cl_qinf = Fn .* cos(deg2rad(alpha)) - Fa .* sin(deg2rad(alpha));
   Cl_qinf = Cl_qinf./(qinf300(i)^2.*s);
   Cl_err = Cl .* sqrt((Cl_Fn.*del_Fn).^2+(Cl_Fa.*del_Fa).^2+(Cl_alpha.*del_alpha).^2+(Cl_qinf.*del_qinf).^2);
   Cd_Fn = sin(deg2rad(alpha))./(qinf300(i).*s);
   Cd_Fa = cos(deg2rad(alpha))./(qinf300(i).*s);
   Cd_alpha = Fn .* cos(deg2rad(alpha)) - Fa .* sin(deg2rad(alpha));
   Cd_alpha = Cd_alpha./(qinf300(i).*s);
   Cd_qinf = -Fn .* cos(deg2rad(alpha)) - Fa .* sin(deg2rad(alpha));
   Cd_qinf = Cd_qinf./(qinf300(i)^2.*s);
   Cd_err = Cd .* sqrt((Cd_Fn.*del_Fn).^2+(Cd_Fa.*del_Fa).^2+(Cd_alpha.*del_alpha).^2+(Cd_qinf.*del_qinf).^2);
   % Store the Cl and alpha values for this file
   Cl_err_all{i} = Cl_err;
   Cd_err_all{i} = Cd_err;
   Cd_all{i} = Cd;
   Cl_all{i} = Cl;
   alpha_all{i} = alpha;
end
% Plot all three Cl vs. alpha lines on the same plot
figure(3);
hold on;  % Keep the plot open to overlay each line
for i = 1:length(angleRe150)
   plot(cell2mat(alpha_all(i)), cell2mat(Cl_all(i)), 'DisplayName', sprintf('File %d', i));  % Use legend to identify each line
   errorbar(cell2mat(alpha_all(i)), cell2mat(Cl_all(i)), cell2mat(Cl_err_all(i)));
end
hold off;
% Label the axes and add a legend
xlabel('\alpha (Angle of Attack)');
ylabel('C_L (Lift Coefficient)');
title('C_L vs. \alpha for Re 300K for Winglet Angles');
legend('None','20 degrees','40 degrees','60 degrees');  % Show legend with labels for each line
grid on;
% Plot all three Cl vs. Cd lines on the same plot
figure(4);
hold on;  % Keep the plot open to overlay each line
for i = 1:length(angleRe300)
   plot(cell2mat(Cd_all(i)), cell2mat(Cl_all(i)), 'DisplayName', sprintf('File %d', i));  % Use legend to identify each line
   errorbar(cell2mat(Cd_all(i)), cell2mat(Cl_all(i)), cell2mat(Cl_err_all(i)), 'both');
end
hold off;
% Label the axes and add a legend
xlabel('C_D (Drag Coefficient)');
ylabel('C_L (Lift Coefficient)');
title('C_L vs. C_D for Re 300K for Winglet Angles');
legend('None','20 degrees','40 degrees','60 degrees');  % Show legend with labels for each line
grid on;
clear all;
gravityG2 ={'Lab4_G3_none_grav.csv','Lab4_G3_back_grav.csv','Lab4_G3_mid_grav.csv','Lab4_G3_front_grav.csv'};
posRe150 = {'Lab4_G3_none_150k.csv','Lab4_G3_back_150k.csv','Lab4_G3_mid_150k.csv','Lab4_G3_front_150k.csv'};
posRe300 = {'Lab4_G3_none_300k.csv','Lab4_G3_back_300k.csv','Lab4_G3_mid_300k.csv','Lab4_G3_front_300k.csv'};
D = readmatrix('Lab4_G3_Conditions.txt');
P = D(:,3);
P = P * 1000;          % Convert pressure from kPa to Pa
T = D(:,2);       % Extract 2nd column for temperature in Celsius
T = T + 273.15;        % Convert temperature to Kelvin
R = 287.05;            % Specific gas constant for dry air in J/(kg*K)
rho = P ./ (R .* T);    % Calculate density using the ideal gas law
sizeT = length(T);      % Get the length of the temperature array
mu = zeros(1, sizeT);   % Preallocate mu as a numeric array (list) with zeros
for i = 1:sizeT
   mu(i) = 1.716e-5 * ((273.15 + 110.4) / (T(i) + 110.4)) * (T(i) / 273.15)^(3/2);
end
c = 0.1397;
b = 0.1524;
s = c*b;
Re150 = D(:,4);
Re300 = D(:,5);
V150 = Re150.*mu./(rho*c);
V300 = Re300.*mu./(rho*c);
for i = 1:sizeT
   qinf150(i) = 0.5*rho(i)*V150(i)^2;
end
for i = 1:sizeT
   qinf300(i) = 0.5*rho(i)*V300(i)^2;
end
%%Angles, Cl, and Cd with Re 150K
%------------------------------------------------------------------------
% Initialize cell arrays to store Cl and alpha values for each file
Cl_all = cell(1, length(posRe150));
alpha_all = cell(1, length(posRe150));
% Loop through each file and calculate Cl and Cd
for i = 1:length(posRe150)
   D = readmatrix(posRe150{i},'delimiter', ",");
   G = readmatrix(gravityG2{i},'delimiter', ",");
   Fa = D(:,7)-G(:,7);
   Fn = D(:,8)-G(:,8);
   alpha = D(:,10);
  
   % Calculate Lift (L) and Drag (D) forces
   L = Fn .* cos(deg2rad(alpha)) - Fa .* sin(deg2rad(alpha));
   D = Fn .* sin(deg2rad(alpha)) + Fa .* cos(deg2rad(alpha));
  
   % Calculate Cl and Cd
   Cl = L ./ (qinf150(i).*s);
   Cd = D ./ (qinf150(i).*s);
  
   del_Fn = 0.05;
   del_Fa = 0.05;
   del_alpha = deg2rad(0.05);
   del_qinf = 0.5;
   Cl_Fn = cos(deg2rad(alpha))./(qinf150(i).*s);
   Cl_Fa = -sin(deg2rad(alpha))./(qinf150(i).*s);
   Cl_alpha = -Fn .* sin(deg2rad(alpha)) - Fa .* cos(deg2rad(alpha));
   Cl_alpha = Cl_alpha./(qinf150(i).*s);
   Cl_qinf = Fn .* cos(deg2rad(alpha)) - Fa .* sin(deg2rad(alpha));
   Cl_qinf = Cl_qinf./(qinf150(i)^2.*s);
   Cl_err = Cl .* sqrt((Cl_Fn.*del_Fn).^2+(Cl_Fa.*del_Fa).^2+(Cl_alpha.*del_alpha).^2+(Cl_qinf.*del_qinf).^2);
   Cd_Fn = sin(deg2rad(alpha))./(qinf150(i).*s);
   Cd_Fa = cos(deg2rad(alpha))./(qinf150(i).*s);
   Cd_alpha = Fn .* cos(deg2rad(alpha)) - Fa .* sin(deg2rad(alpha));
   Cd_alpha = Cd_alpha./(qinf150(i).*s);
   Cd_qinf = -Fn .* cos(deg2rad(alpha)) - Fa .* sin(deg2rad(alpha));
   Cd_qinf = Cd_qinf./(qinf150(i)^2.*s);
   Cd_err = Cd .* sqrt((Cd_Fn.*del_Fn).^2+(Cd_Fa.*del_Fa).^2+(Cd_alpha.*del_alpha).^2+(Cd_qinf.*del_qinf).^2);
   % Store the Cl and alpha values for this file
   Cl_err_all{i} = Cl_err;
   Cd_err_all{i} = Cd_err;
   Cd_all{i} = Cd;
   Cl_all{i} = Cl;
   alpha_all{i} = alpha;
end
% Plot all three Cl vs. alpha lines on the same plot
figure(1);
hold on;  % Keep the plot open to overlay each line
for i = 1:length(posRe150)
   plot(cell2mat(alpha_all(i)), cell2mat(Cl_all(i)), 'DisplayName', sprintf('File %d', i));  % Use legend to identify each line
   errorbar(cell2mat(alpha_all(i)), cell2mat(Cl_all(i)), cell2mat(Cl_err_all(i)));
end
hold off;
% Label the axes and add a legend
xlabel('\alpha (Angle of Attack)');
ylabel('C_L (Lift Coefficient)');
title('C_L vs. \alpha for Re 150K for Winglet Locations');
legend('None','Back','Mid','Front');  % Show legend with labels for each line
grid on;
% Plot all three Cl vs. Cd lines on the same plot
figure(2);
hold on;  % Keep the plot open to overlay each line
for i = 1:length(posRe150)
   plot(cell2mat(Cd_all(i)), cell2mat(Cl_all(i)), 'DisplayName', sprintf('File %d', i));  % Use legend to identify each line
   errorbar(cell2mat(Cd_all(i)), cell2mat(Cl_all(i)), cell2mat(Cl_err_all(i)), 'both');
end
hold off;
% Label the axes and add a legend
xlabel('C_D (Drag Coefficient)');
ylabel('C_L (Lift Coefficient)');
title('C_L vs. C_D for Re 150K for Winglet Locations');
legend('None','Back','Mid','Front');  % Show legend with labels for each line
grid on;
%%Angles, Cl, and Cd with Re 300K
%-----------------------------------------------------------------------
% Initialize cell arrays to store Cl and alpha values for each file
Cl_all = cell(1, length(posRe300));
alpha_all = cell(1, length(posRe300));
% Loop through each file and calculate Cl and Cd
for i = 1:length(posRe300)
   D = readmatrix(posRe300{i},'delimiter', ",");
   G = readmatrix(gravityG2{i},'delimiter', ",");
   Fa = D(:,7)-G(:,7);
   Fn = D(:,8)-G(:,8);
   alpha = D(:,10);
  
   % Calculate Lift (L) and Drag (D) forces
   L = Fn .* cos(deg2rad(alpha)) - Fa .* sin(deg2rad(alpha));
   D = Fn .* sin(deg2rad(alpha)) + Fa .* cos(deg2rad(alpha));
  
   % Calculate Cl and Cd
   Cl = L ./ (qinf300(i).*s);
   Cd = D ./ (qinf300(i).*s);
  
   del_Fn = 0.05;
   del_Fa = 0.05;
   del_alpha = deg2rad(0.05);
   del_qinf = 0.5;
   Cl_Fn = cos(deg2rad(alpha))./(qinf300(i).*s);
   Cl_Fa = -sin(deg2rad(alpha))./(qinf300(i).*s);
   Cl_alpha = -Fn .* sin(deg2rad(alpha)) - Fa .* cos(deg2rad(alpha));
   Cl_alpha = Cl_alpha./(qinf300(i).*s);
   Cl_qinf = Fn .* cos(deg2rad(alpha)) - Fa .* sin(deg2rad(alpha));
   Cl_qinf = Cl_qinf./(qinf300(i)^2.*s);
   Cl_err = Cl .* sqrt((Cl_Fn.*del_Fn).^2+(Cl_Fa.*del_Fa).^2+(Cl_alpha.*del_alpha).^2+(Cl_qinf.*del_qinf).^2);
   Cd_Fn = sin(deg2rad(alpha))./(qinf300(i).*s);
   Cd_Fa = cos(deg2rad(alpha))./(qinf300(i).*s);
   Cd_alpha = Fn .* cos(deg2rad(alpha)) - Fa .* sin(deg2rad(alpha));
   Cd_alpha = Cd_alpha./(qinf300(i).*s);
   Cd_qinf = -Fn .* cos(deg2rad(alpha)) - Fa .* sin(deg2rad(alpha));
   Cd_qinf = Cd_qinf./(qinf300(i)^2.*s);
   Cd_err = Cd .* sqrt((Cd_Fn.*del_Fn).^2+(Cd_Fa.*del_Fa).^2+(Cd_alpha.*del_alpha).^2+(Cd_qinf.*del_qinf).^2);
   % Store the Cl and alpha values for this file
   Cl_err_all{i} = Cl_err;
   Cd_err_all{i} = Cd_err;
   Cd_all{i} = Cd;
   Cl_all{i} = Cl;
   alpha_all{i} = alpha;
end
% Plot all three Cl vs. alpha lines on the same plot
figure(3);
hold on;  % Keep the plot open to overlay each line
for i = 1:length(posRe150)
   plot(cell2mat(alpha_all(i)), cell2mat(Cl_all(i)), 'DisplayName', sprintf('File %d', i));  % Use legend to identify each line
   errorbar(cell2mat(alpha_all(i)), cell2mat(Cl_all(i)), cell2mat(Cl_err_all(i)));
end
hold off;
% Label the axes and add a legend
xlabel('\alpha (Angle of Attack)');
ylabel('C_L (Lift Coefficient)');
title('C_L vs. \alpha for Re 300K for Winglet Locations');
legend('None','Back','Mid','Front');  % Show legend with labels for each line
grid on;
% Plot all three Cl vs. Cd lines on the same plot
figure(4);
hold on;  % Keep the plot open to overlay each line
for i = 1:length(posRe300)
   plot(cell2mat(Cd_all(i)), cell2mat(Cl_all(i)), 'DisplayName', sprintf('File %d', i));  % Use legend to identify each line
   errorbar(cell2mat(Cd_all(i)), cell2mat(Cl_all(i)), cell2mat(Cl_err_all(i)), 'both');
end
hold off;
% Label the axes and add a legend
xlabel('C_D (Drag Coefficient)');
ylabel('C_L (Lift Coefficient)');
title('C_L vs. C_D for Re 300K for Winglet Locations');
legend('None','Back','Mid','Front');  % Show legend with labels for each line
grid on;
clear all;
gravityG3 = {'Lab4_G3_none_grav.csv','Lab4_G3_long_grav.csv', 'Lab4_G3_med_grav.csv','Lab4_G3_short_grav.csv'};
sizeRe150 = {'Lab4_G3_none_150k.csv','Lab4_G3_long_150k.csv', 'Lab4_G3_med_150k.csv','Lab4_G3_short_150k.csv'};
sizeRe300 = {'Lab4_G3_none_300k.csv','Lab4_G3_long_300k.csv', 'Lab4_G3_med_300k.csv','Lab4_G3_small_300k.csv'};
D = readmatrix('Lab4_G3_Conditions.txt');
P = D(:,3);
P = P * 1000;          % Convert pressure from kPa to Pa
T = D(:,2);       % Extract 2nd column for temperature in Celsius
T = T + 273.15;        % Convert temperature to Kelvin
R = 287.05;            % Specific gas constant for dry air in J/(kg*K)
rho = P ./ (R .* T);    % Calculate density using the ideal gas law
sizeT = length(T);      % Get the length of the temperature array
mu = zeros(1, sizeT);   % Preallocate mu as a numeric array (list) with zeros
for i = 1:sizeT
   mu(i) = 1.716e-5 * ((273.15 + 110.4) / (T(i) + 110.4)) * (T(i) / 273.15)^(3/2);
end
c = 0.1397;
b = 0.1524;
s = c*b;
Re150 = D(:,4);
Re300 = D(:,5);
V150 = Re150.*mu./(rho*c);
V300 = Re300.*mu./(rho*c);
for i = 1:sizeT
   qinf150(i) = 0.5*rho(i)*V150(i)^2;
end
for i = 1:sizeT
   qinf300(i) = 0.5*rho(i)*V300(i)^2;
end
%%Angles, Cl, and Cd with Re 150K
%------------------------------------------------------------------------
% Initialize cell arrays to store Cl and alpha values for each file
Cl_all = cell(1, length(sizeRe150));
alpha_all = cell(1, length(sizeRe150));
% Loop through each file and calculate Cl and Cd
for i = 1:length(sizeRe150)
   D = readmatrix(sizeRe150{i},'delimiter', ",");
   G = readmatrix(gravityG3{i},'delimiter', ",");
   Fa = D(:,7)-G(:,7);
   Fn = D(:,8)-G(:,8);
   alpha = D(:,10);
  
   % Calculate Lift (L) and Drag (D) forces
   L = Fn .* cos(deg2rad(alpha)) - Fa .* sin(deg2rad(alpha));
   D = Fn .* sin(deg2rad(alpha)) + Fa .* cos(deg2rad(alpha));
  
   % Calculate Cl and Cd
   Cl = L ./ (qinf150(i).*s);
   Cd = D ./ (qinf150(i).*s);
  
   del_Fn = 0.05;
   del_Fa = 0.05;
   del_alpha = deg2rad(0.05);
   del_qinf = 0.5;
   Cl_Fn = cos(deg2rad(alpha))./(qinf150(i).*s);
   Cl_Fa = -sin(deg2rad(alpha))./(qinf150(i).*s);
   Cl_alpha = -Fn .* sin(deg2rad(alpha)) - Fa .* cos(deg2rad(alpha));
   Cl_alpha = Cl_alpha./(qinf150(i).*s);
   Cl_qinf = Fn .* cos(deg2rad(alpha)) - Fa .* sin(deg2rad(alpha));
   Cl_qinf = Cl_qinf./(qinf150(i)^2.*s);
   Cl_err = Cl .* sqrt((Cl_Fn.*del_Fn).^2+(Cl_Fa.*del_Fa).^2+(Cl_alpha.*del_alpha).^2+(Cl_qinf.*del_qinf).^2);
   Cd_Fn = sin(deg2rad(alpha))./(qinf150(i).*s);
   Cd_Fa = cos(deg2rad(alpha))./(qinf150(i).*s);
   Cd_alpha = Fn .* cos(deg2rad(alpha)) - Fa .* sin(deg2rad(alpha));
   Cd_alpha = Cd_alpha./(qinf150(i).*s);
   Cd_qinf = -Fn .* cos(deg2rad(alpha)) - Fa .* sin(deg2rad(alpha));
   Cd_qinf = Cd_qinf./(qinf150(i)^2.*s);
   Cd_err = Cd .* sqrt((Cd_Fn.*del_Fn).^2+(Cd_Fa.*del_Fa).^2+(Cd_alpha.*del_alpha).^2+(Cd_qinf.*del_qinf).^2);
   % Store the Cl and alpha values for this file
   Cl_err_all{i} = Cl_err;
   Cd_err_all{i} = Cd_err;
   Cd_all{i} = Cd;
   Cl_all{i} = Cl;
   alpha_all{i} = alpha;
end
% Plot all three Cl vs. alpha lines on the same plot
figure(1);
hold on;  % Keep the plot open to overlay each line
for i = 1:length(sizeRe150)
   plot(cell2mat(alpha_all(i)), cell2mat(Cl_all(i)), 'DisplayName', sprintf('File %d', i));  % Use legend to identify each line
   errorbar(cell2mat(alpha_all(i)), cell2mat(Cl_all(i)), cell2mat(Cl_err_all(i)));
end
hold off;
% Label the axes and add a legend
xlabel('\alpha (Angle of Attack)');
ylabel('C_L (Lift Coefficient)');
title('C_L vs. \alpha for Re 150K for Winglet Size');
legend('None','Long','Med','Short');  % Show legend with labels for each line
grid on;
% Plot all three Cl vs. Cd lines on the same plot
figure(2);
hold on;  % Keep the plot open to overlay each line
for i = 1:length(sizeRe150)
   plot(cell2mat(Cd_all(i)), cell2mat(Cl_all(i)), 'DisplayName', sprintf('File %d', i));  % Use legend to identify each line
   errorbar(cell2mat(Cd_all(i)), cell2mat(Cl_all(i)), cell2mat(Cl_err_all(i)), 'both');
end
hold off;
% Label the axes and add a legend
xlabel('C_D (Drag Coefficient)');
ylabel('C_L (Lift Coefficient)');
title('C_L vs. C_D for Re 150K for Winglet Size');
legend('None','Long','Med','Short');  % Show legend with labels for each line
grid on;
%%Angles, Cl, and Cd with Re 300K
%-----------------------------------------------------------------------
% Initialize cell arrays to store Cl and alpha values for each file
Cl_all = cell(1, length(sizeRe300));
alpha_all = cell(1, length(sizeRe300));
% Loop through each file and calculate Cl and Cd
for i = 1:length(sizeRe300)
   D = readmatrix(sizeRe300{i},'delimiter', ",");
   G = readmatrix(gravityG3{i},'delimiter', ",");
   Fa = D(:,7)-G(:,7);
   Fn = D(:,8)-G(:,8);
   alpha = D(:,10);
  
   % Calculate Lift (L) and Drag (D) forces
   L = Fn .* cos(deg2rad(alpha)) - Fa .* sin(deg2rad(alpha));
   D = Fn .* sin(deg2rad(alpha)) + Fa .* cos(deg2rad(alpha));
  
   % Calculate Cl and Cd
   Cl = L ./ (qinf300(i).*s);
   Cd = D ./ (qinf300(i).*s);
  
   del_Fn = 0.05;
   del_Fa = 0.05;
   del_alpha = deg2rad(0.05);
   del_qinf = 0.5;
   Cl_Fn = cos(deg2rad(alpha))./(qinf300(i).*s);
   Cl_Fa = -sin(deg2rad(alpha))./(qinf300(i).*s);
   Cl_alpha = -Fn .* sin(deg2rad(alpha)) - Fa .* cos(deg2rad(alpha));
   Cl_alpha = Cl_alpha./(qinf300(i).*s);
   Cl_qinf = Fn .* cos(deg2rad(alpha)) - Fa .* sin(deg2rad(alpha));
   Cl_qinf = Cl_qinf./(qinf300(i)^2.*s);
   Cl_err = Cl .* sqrt((Cl_Fn.*del_Fn).^2+(Cl_Fa.*del_Fa).^2+(Cl_alpha.*del_alpha).^2+(Cl_qinf.*del_qinf).^2);
   Cd_Fn = sin(deg2rad(alpha))./(qinf300(i).*s);
   Cd_Fa = cos(deg2rad(alpha))./(qinf300(i).*s);
   Cd_alpha = Fn .* cos(deg2rad(alpha)) - Fa .* sin(deg2rad(alpha));
   Cd_alpha = Cd_alpha./(qinf300(i).*s);
   Cd_qinf = -Fn .* cos(deg2rad(alpha)) - Fa .* sin(deg2rad(alpha));
   Cd_qinf = Cd_qinf./(qinf300(i)^2.*s);
   Cd_err = Cd .* sqrt((Cd_Fn.*del_Fn).^2+(Cd_Fa.*del_Fa).^2+(Cd_alpha.*del_alpha).^2+(Cd_qinf.*del_qinf).^2);
   % Store the Cl and alpha values for this file
   Cl_err_all{i} = Cl_err;
   Cd_err_all{i} = Cd_err;
   Cd_all{i} = Cd;
   Cl_all{i} = Cl;
   alpha_all{i} = alpha;
end
% Plot all three Cl vs. alpha lines on the same plot
figure(3);
hold on;  % Keep the plot open to overlay each line
for i = 1:length(sizeRe150)
   plot(cell2mat(alpha_all(i)), cell2mat(Cl_all(i)), 'DisplayName', sprintf('File %d', i));  % Use legend to identify each line
   errorbar(cell2mat(alpha_all(i)), cell2mat(Cl_all(i)), cell2mat(Cl_err_all(i)));
end
hold off;
% Label the axes and add a legend
xlabel('\alpha (Angle of Attack)');
ylabel('C_L (Lift Coefficient)');
title('C_L vs. \alpha for Re 300K for Winglet Size');
legend('None','Long','Med','Short');  % Show legend with labels for each line
grid on;
% Plot all three Cl vs. Cd lines on the same plot
figure(4);
hold on;  % Keep the plot open to overlay each line
for i = 1:length(sizeRe300)
   plot(cell2mat(Cd_all(i)), cell2mat(Cl_all(i)), 'DisplayName', sprintf('File %d', i));  % Use legend to identify each line
   errorbar(cell2mat(Cd_all(i)), cell2mat(Cl_all(i)), cell2mat(Cl_err_all(i)), 'both');
end
hold off;
% Label the axes and add a legend
xlabel('C_D (Drag Coefficient)');
ylabel('C_L (Lift Coefficient)');
title('C_L vs. C_D for Re 300K for Winglet Size');
legend('None','Long','Med','Short');  % Show legend with labels for each line
grid on;

