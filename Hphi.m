%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Author: M. Stumpf                %
%  Date: 26 May 2026                %
%  Revised: 22 October 2026         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This m-file evaluates the azimuthal component of the H-field outside
% a thin spherical shell as radiated by a vertical electric dipole 
% located inside the shell
%
% To get the field in the absence of the shell, set YL = 0.
%
clear; clc;
%
%%%%%%%%%%%%%%%%
% EM CONSTANTS %
%%%%%%%%%%%%%%%%
%
c0 = 299792458;
mu0 = 4*pi*1e-7;
%
ep0 = 1/c0^2/mu0;
%
Y0 = 1/c0/mu0;
%
%%%%%%%%%%%%%%%%%%%%%
% SPHERE PARAMETERS %
%%%%%%%%%%%%%%%%%%%%%
%
% Radius
R0 = 0.050; 
% Layer admittance
sigma = 5.65e+7; 
delta = 1.0e-6;
YL = sigma*delta;
%
% absence of the sphere
% YL = 0;   
%
% ----------------- DISTANCE from the sphere's center ---------------------
R = 1.10*R0;
%
% TIME AXIS
Tmax = 6*R0/c0; M = 2^12; dt = Tmax/M;
t = (0:M-1)*dt;
%
% REAL-FREQUENCY AXIS
df = 1/Tmax; f = (-M/2:M/2-1)*df;
omega = 2*pi*f;
%
% COMPLEX FREQUENCY
alpha = 12/Tmax; s = alpha + 1i*omega;
%
% -------------------------------------------------------------------------
% ----------------------- s-DOMAIN FUNCTION -------------------------------
% -------------------------------------------------------------------------
%
eta0 = s*ep0; zeta0 = s*mu0; gamma0 = s/c0;
%
Kp = gamma0.^2*R0^2 + gamma0*R0 + 1;
Km = gamma0.^2*R0^2 - gamma0*R0 + 1;
%
% TRANSMISSION coefficient (without the exp-delay term)
N = 2*eta0.*gamma0.^3.*R0^4;
%
D = N + YL*(Kp.*Km) - YL*Kp.^2.*exp(-2.*gamma0.*R0);
D(abs(D) < 1e-20) = 1e-20;
%
TT = N./D;
%
%
%
% REFLECTION coefficient (without the exp-delay term)
N = -YL*Kp.^2;
RR = N./D;
%
%
% ------------------------- ELECTRIC DIPOLE -------------------------------
% PE EXCITATION pulse shape
tw = 0.5*R0/c0; Im = 1; 
power = 4; ratio = power^(-power-1)*gamma(power+1)*exp(power);
tr = tw/ratio;     
I = Im*tr*gamma(power+1)*exp(power)./(s*tr + power).^(power+1);
%
% length of the dipole
ell = 1.0e-3;
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------------ s-domain POTENTIAL FUNCTIONS ----------------------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: With the total exp/delay factors
%
% ------------------ in R > R0 (outside the sphere) -----------------------
Phi1 = TT.*I*ell.*exp(-gamma0*R)/4/pi/R;
% 1st R-derivative of the potential function
dPhi1 = -Phi1.*(gamma0 + 1/R); 
% 1st R-derivative of the potential function
d2Phi1 = Phi1.*(gamma0.^2 + 2*gamma0/R + 2/R^2); 
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------- s-domain EM FIELDS ----------------------------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% polar angle of observation
theta = 1*pi/2; 
%
% ------------------ in R > R0 (outside the sphere) -----------------------
% azimuthal component of H-field
Hphi1 = - dPhi1*sin(theta);
%
% ---------------------------------------
% FFT-based Laplace-transform inversion %
% ---------------------------------------
%
Fifft = ifft(ifftshift(Hphi1));
Ft = real(Fifft) * (M * df);
Ft = Ft .* exp(alpha*t);
%
disp(['Time required: ' num2str(toc) ' (s)']);
%
figure(3); hold on; plot(c0*t/R0, Ft, 'LineWidth', 1.0); grid on
xlabel('c_0t/R0 (-)'); ylabel('Hphi(t)');
%
%
%
%
%
%
%
%
%
%
%