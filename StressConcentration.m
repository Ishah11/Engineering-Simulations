%This experiment investigates stress concentrations in polycarbonate plates with varying hole diameters subjected to uniaxial tensile loading.
clear; clc; clear all;
files = {'00625.txt','0125.txt','0250.txt','0500.txt','0750.txt','1000.txt'};
P = [4e3,4e3,4e3,3.5e3,2.9e3,2.3e3];
t = [3.005e-3,2.99e-3,3.01e-3,3.01e-3,3.015e-3,3.02e-3];
w = [50.88e-3,50.87e-3,50.82e-3,50.98e-3,50.83e-3,50.61e-3];
d = [1.5875e-3,3.175e-3,6.35e-3,12.7e-3,19.05e-3,25.4e-3];
hole = {'1/16','1/8','1/4','1/2','3/4','1'};
E = zeros(1,length(files));
steve = zeros(1,length(files));
Ktexp = zeros(1,length(files));
Ktanalytical = [2.8943,2.7985,2.6422,2.4297,2.2877,2.1694];
for i = 1:length(files)
   % Read entire file into a cell array
   C = readcell(files{i}, 'Delimiter', '\t');
   farfieldRow = find(strcmp(C, 'Farfield')) + 1;
   stressRow = find(strcmp(C, 'Stress Concentration')) + 1;
   farfieldData = cell2mat(C(farfieldRow+1:farfieldRow+8, 1:2)); 
   indexFar = farfieldData(:,1);
   strainFar = farfieldData(:,2)./100;
   stressData = cell2mat(C(stressRow+1:stressRow+8, 1:3)); 
   indexStress = stressData(:,1);
   distanceStress = (w(i)/2)-stressData(:,2)./1000;
   strainStress = stressData(:,3)./100;
  
   A = w(i)*t(i);
   E1 = (P(i)/A)./strainFar;
  
   E(i) = mean(E1);
   steve(i) = std(E1);
   % a = radius of hole
   a = 0.5 * d(i);
  
   % y = distance from center of hole (excluding inside the hole, start at edge)
   y = linspace(a, 0.5 * w(i), 100);  % outward from hole edge to edge of plate
  
   % Equation 3: Theoretical stress at y
   sigma_theoretical = 0.5 * (P(i)/A) * (1 + (a^2 ./ y.^2)) + 0.5 * (P(i)/A) * (1 + (3 * a^4 ./ y.^4));
   stress = E(i)*strainStress;
   figure(i);
   plot(distanceStress, stress, 'o-', 'DisplayName', 'Experimental');
   hold on;
   plot(y, sigma_theoretical, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Theoretical');
   xlabel('Position from center (m)');
   ylabel('Stress (Pa)');
   title(['Stress Distribution for D = ', hole{i}]);
   grid on;
   legend('show');
   % Net section area and applied stress
   % Far-field nominal stress
   Area = (w(i) - d(i)) * t(i);
   sigma0 = P(i) / Area;
  
   % Max stress at hole (from stress region)
   sigma_max = max(strainStress) * E(i);  % strainStress is from near hole
  
   % Experimental Kt
   Ktexp(i) = sigma_max / sigma0;
end
percent_diff = abs(Ktanalytical-Ktexp)./Ktanalytical;
percent_diff = percent_diff .* 100;
d_w = d./w;
figure();
plot(d_w,Ktexp,'DisplayName','Experimental');
hold on;
plot(d_w,Ktanalytical,'DisplayName','Analytical')
xlabel('d/w');
ylabel('K_t');
title('Stress Concentration Factor');
legend('show');
grid on;

