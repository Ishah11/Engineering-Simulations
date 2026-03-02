%This lab investigates numerical solutions to the Small Disturbance Equations for both subsonic and transonic flows over a symmetric circular-arc airfoil at zero angle of attack. 
clear; clc; close all;
I = 102;
J = I;
Ks = [3,1.3];
tau = 0.1;
titles = {'Subsonic','Transonic'};
data = {'TSDEdata3.txt','TSDEdata1_3.txt'};
incomp = 'incompr.txt';
Di = importdata(incomp);
x_c = Di.data(:,1);
Cpi = Di.data(:,2)*tau^(-2/3);
Z = @(x) 2.*tau.*(0.25-(x-0.5).^2);
incompz = Z(x_c);
for n = 1:2
   K = @(M) (1-M^2)/((M^2*tau)^(2/3))-Ks(n);
   M(n) = fsolve(K,0.8);
   [Cp,x] = JacobiIteration(M(n),tau);
   D = importdata(data{n});
   TSDE_x = D.data(:,1);
   TSDE_cp = D.data(:,2);
   figure(n)
   hold on;
   axis ij;
   plot(x,Cp,'-k','DisplayName','Prandtle-Glauert Solution');
   plot(TSDE_x,TSDE_cp,'b','DisplayName','TSDE Data');
   plot(x_c,Cpi,'r','DisplayName','Incompressible Data');
   xlabel('Dimesionless Position x');
   ylabel('Coefficient of Pressure  C_p');
   title(sprintf('%s K = %.1f', titles{n}, Ks(n)));
   legend('show');
  
   z = Z(x);
   TSDE_z = Z(TSDE_x);
   Jacobi_cd(:,n) = 2*trapz(z,Cp);
   TSDE_cd(:,n) = 2*trapz(TSDE_z,TSDE_cp);
   incomp_cd(:,n) = 2*trapz(incompz,Cpi);
end
T = table([Jacobi_cd(1);TSDE_cd(1);incomp_cd(1)],[Jacobi_cd(2);TSDE_cd(2);incomp_cd(2)], RowNames={'Prandtl-Glauert', 'TSDE', 'Incompressible'}, ...
   VariableNames={'K = 3','K = 1.3'});
disp(T)
function[Cp,x]=  JacobiIteration(M,tau)
   B = sqrt(1-M^2);
   x = linspace(-0.5,1.5);
   dx = mean(diff(x));
   I = length(x);
   J = I;
   dz = dx;
   error = 1;
   Error = zeros(1,10^5);
   k = 0;
   phi = zeros(I,J);
   phi0 = zeros(I,J);
   MaxError = 10^(-6);
   while error > MaxError
       k = k+1;
       if k == 10000
           disp('did not converge')
           break
       end
       for j = 1:J
           phi(1,j) = phi(3,j);
           phi(I,j) = phi(I-2,j);
       end
       for i = 2:I-1
           phi(i,J) = phi(i,J-2);
           if  x(i) < 0 || x(i) > 1 %off body conditions
               phi(i,1) = phi(i,3);
           else
               z = 2*tau*(1-2*x(i));
               phi(i,1) = phi(i,3) - 2*dz*z; %on body boundary conditions
           end
       end
   phi_old = phi0;
%find new values for interiro poitns
   for i = 2:I-1
       for j = 2:J-1
           Jacobi = (B^2*(phi0(i-1,j)+phi0(i+1,j))/dx^2+(phi0(i,j+1)+phi0(i,j-1))/dz^2)/(2*B^2/dx^2 + 2/dz^2);
           phi(i,j) = Jacobi;
           %phi0(i,j) = phi(i,j); %save val
       end
   end
   %compute error
   error = norm(phi-phi_old);
   Error(k) = error;
   phi0 = phi;
   end
   u(2:I-1) = (phi(3:I,1)-phi(1:I-2,1))/(2*dx);
   Cp = -2*u/tau^(2/3);
   range = (x > 0 & x < 1);
   x = x(range);
   Cp = Cp(range);
end
