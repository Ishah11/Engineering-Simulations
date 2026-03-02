%This experiment investigated the dynamic behavior of a simply supported mass-beam system subjected to both free and forced vibrations.
clear; clc; close all;
%% Data
Mt = 1.3129;
Mb = 0.6055;
l_star = 0.0375; %m
E = 7.31e10;
l = 0.540;
b = 0.05;
h = 0.006;
b_star = 0.084;
h_star = 0.006;
I = 1/12*b*h^3;
I_star = 1/12*b_star*h_star^3;
F = [1,2,3,4,5]*9.81;
LVDT = [0.39,0.77,1.16,1.54,1.92]*10^(-3);
m = 15.1e-3;
e = 31.75e-3;
M = Mt+17/35*Mb;
%% Part 1
keq_exp = mean(F./LVDT);
Rei = (E*I)/(E*I_star);
keq_theor = 48*E*I/(Rei*l^3+(1-Rei)*(l-l_star)^3);
k_err = 100*(keq_exp-keq_theor)/keq_theor;
%% Part 2
D = importdata('F_1200_G1_free.txt');
time = D.data(:,1);
pos = D.data(:,2);
figure(1);
plot(time,pos,'-b','LineWidth', 1.5);
xlabel('time (sec)');
ylabel('Amplitude (mm)');
title('Damped Free Response');
grid on;
[peaks, locs] = findpeaks(pos,time,'MinPeakDistance',0.045);
first10 = peaks(1:min(10,end));
deltaL = zeros(1,9);
for i = 1:9
   deltaL(i) = log(first10(i)/first10(i+1));
end
delta_avg = 1/9*sum(deltaL);
zeta = @(zeta)2*pi*zeta/sqrt(1-zeta^2)-delta_avg;
damping_ratio = fzero(zeta,0.5);
%% Part 3
first10_times = diff(locs(1:min(10, end)));
omega_d = 2*pi/mean(first10_times);
omega_n = omega_d/sqrt(1-damping_ratio^2);
omegan_theo = sqrt(keq_theor/(M));
omega_n_err = 100*(omega_n-omegan_theo)/omegan_theo;
%% Part 4
Data = importdata('F_1200_G1_forced.txt');
freq = Data.data(:,1).*(2.*pi)./60;
amp = Data.data(:,2).*10^(-3);
phi = Data.data(:,3).*(pi/180);
[amp_max, idx_max] = max(amp); 
omega_res = freq(idx_max);
omega_res_theor = omegan_theo/sqrt(1-2*damping_ratio^2);
omega_res_err = 100*(omega_res-omega_res_theor)/omega_res_theor;
%% Part 5
r = freq./omega_n;
c = 2*damping_ratio*m*omega_n;
X = amp;
MX_me = M*X/(m*e);
MX_me_theor = r.^2./sqrt((1-r.^2).^2+(2*damping_ratio*r).^2);
figure(2);
plot(r,MX_me,'DisplayName','Experimental','LineWidth',2);
hold on;
plot(r,MX_me_theor,'DisplayName','Theoretical','LineWidth',2);
xlabel('Dimensionless Excitation Frequency');
ylabel('Dimensionless Amplitude');
title('Damped Force Response Data');
legend('show')
grid on;
