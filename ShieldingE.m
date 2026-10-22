%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Author: M. Stumpf                %
%  Date: 21 May 2026                %
%  Revised: 22 October 2026         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This m-file evaluates ANALYTICALLY E-type shielding efficiency 
% of a conductive spherical shell
%
close all
%
clc; clear all;
%
%
%%%%%%%%%%%%%%%%
% EM CONSTANTS %
%%%%%%%%%%%%%%%%
%
c0 = 299792458;
mu0 = 4*pi*1e-7;
%
ep0 = 1/c0^2/mu0;
Y0 = 1/c0/mu0;
%
%
% Radius of the spherical shield
R0 = 0.050;
%
% Layer admittance
sigma = 5.65e+7;                      % Copper
delta = 1*1.0e-6;
YL = sigma*delta;
%
% Frequency range
K = 5000; lambda = linspace(10, 1/2, K)*R0;
freq = c0./lambda; om = 2*pi*freq;
%
s = 1i*om + 1*1e-12;
%
eta0 = s*ep0; zeta0 = s*mu0;
gamma0 = s/c0;
%
% -------------------------------------------------------------------------
% ----------------- E-type TRANSMISSION COEFFICIENT -----------------------
% -------------------------------------------------------------------------
%
Kp = (gamma0*R0).^2 + gamma0*R0 + 1;
Km = (gamma0*R0).^2 - gamma0*R0 + 1;
%
KpKm = (gamma0*R0).^4 + (gamma0*R0).^2 + 1;
%
% TRANSMISSION COEFFICIENT (E-field)
D = 2*eta0*R0.*(gamma0*R0).^3 + YL*KpKm - YL*Kp.^2.*exp(-2*gamma0*R0);
TE = 2*eta0*R0.*(gamma0*R0).^3./D.*exp(-gamma0*R0);
%
%
figure(1); plot(R0./lambda, -20*log10(abs(TE))); grid on
xlabel('R0/lambda'); ylabel('SEdB')
%
%
%
%
%
%