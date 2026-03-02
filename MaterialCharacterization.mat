%This experiment conducted a tensile test on 1018 carbon steel and 2024 aluminum to determine Young’s modulus, yield strength, ultimate strength, percent elongation, reduction in area, and power law plasticity coefficients. 
clear;
clc;
% Specify the file name
filename = 'W_0730_G1_Al.TXT';  % Replace with your actual file name
% Read the data
data = readmatrix(filename);
% Extract the columns
al_area = (pi*(2.53/2)^2); %mm
index = data(:,1);
x = data(:,2); %mm
x = x./(12.7); %mm
y = data(:,3); %kN
y = y./al_area; %Gpa
ym_al = polyfit(x(x<= 0.008),y(x<= 0.008),1);
disp(ym_al)
%Yield Stress
al_ysx = linspace(0,0.0255,10000);
al_ysy = al_ysx*ym_al(1)-1.547;
%Ultimate Strength
al_su = max(y);
%Fracture Stress
al_pf = 2.29000;
al_sf = al_pf/al_area;
%Percent Elongation
al_ef = 16.29000;
al_pe = 100*al_ef;
disp(al_pe);
%Reduction Area
al_farea = (pi*(2.13/2)^2);
al_ra = 100*((al_area-al_farea)/al_area);
%True Stress
y1 = data(:,3);
al_yts = abs(y1)./(al_farea./(1+x));
%True Strain
al_xts = log(1+x);
% Plotting the data
figure(1);
plot(x, y, '-b', 'LineWidth', 2);
hold on;
plot(al_ysx(al_ysy>=0),al_ysy(al_ysy>=0),'--g', 'LineWidth', 2)
xlabel('Strain, MPa');
ylabel('Stress, MPa');
title('Aluminum Stress-Strain Diagram');
legend('Original','2% offset');
grid on;
hold off;
figure(2);
plot(al_xts,al_yts,'-r','LineWidth', 2);
xlabel('Strain, MPa');
ylabel('Stress, MPa');
title('Aluminum True Stress-Strain Curve');
grid on;
figure(3);
loglog(al_xts,al_yts,'-r','LineWidth', 2);
xlabel('Strain, MPa');
ylabel('Stress, MPa');
title('Aluminum True Stress-Strain log-log Plot');
grid on;
%-------------------------------------------------------------------------
% Specify the file name
filename = 'AEE_325_W_0730_G1_St.TXT';  % Replace with your actual file name
% Read the data
data = readmatrix(filename);
% Extract the columns
st_area = (pi*(2.51/2)^2); %mm
index = data(:,1);
x = data(:,2); %mm
x = x./(12.7); %mm
y = data(:,3); %kN
y = y./st_area; %Gpa
ym_st = polyfit(x(x<= 0.0016),y(x<= 0.0016),1);
disp(ym_st);
%Yield Stress
st_ysx = linspace(0,0.0225,10000);
st_ysy = st_ysx*ym_st(1)-4.213;
%Ultimate Stress
st_su = max(y);
%Fracture Strength
st_pf = 1.77500;
st_sf = st_pf/st_area;
%Percent elongation
st_ef = 23.95000;
st_pe = 100*st_ef;
%Reduction Area
st_farea = (pi*(1.43/2)^2);
st_ra = 100*((st_area-st_farea)/st_area);
%True Stress
y1 = data(:,3);
st_yts = abs(y1)./(st_farea./(1+x));
%True Strain
st_xts = log(1+x);
% Plotting the data
figure(1);
plot(x, y, '-b', 'LineWidth', 2);
hold on;
plot(st_ysx(st_ysy>=0),st_ysy(st_ysy>=0),'--g','LineWidth', 2)
xlabel('Strain, MPa');
ylabel('Stress, MPa');
title('Steel Stress-Strain Diagram');
legend('Original','2% offset');
grid on;
hold off;
figure(2);
plot(st_xts,st_yts,'-r','LineWidth', 2);
xlabel('Strain, MPa');
ylabel('Stress, MPa');
title('Steel True Stress-Strain Curve');
grid on;
figure(3);
loglog(st_xts,st_yts,'-r','LineWidth', 2);
xlabel('Strain, MPa');
ylabel('Stress, MPa');
title('Steel True Stress-Strain log-log Plot');
grid on;
