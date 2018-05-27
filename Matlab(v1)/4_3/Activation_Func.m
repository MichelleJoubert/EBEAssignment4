
clc;
clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Temporal intervals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T       = 20; %Simulation time in milliseconds
N       = 1000; %Number of simulation points (colums in arrays)
dT      = T/N; %Magnitude of time segments 
t       = 0:dT:(T-dT); %Time axis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spacial intervals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X       = 2; %Length of the axon in cm
R       = 2000; %Amount of segments
dX      = X/R; %Magnitudes of spacial segments
x       = 0:dX:X; %Spacial axis



Ie1     = zeros(R,N+1); %Create the array for the stimulating current
Ie2     = zeros(R,N+1); %Create the array for the stimulating current

Ie1(:,1:201) = 3500; %Current stimulation amplitude 
Ie2(:,1:201) = -3500; %Current stimulation amplitude 

d       = 0.001;
z       = 0.05;
re      = 0.3;
ri      = 0.1;

f_e1    = zeros(R+1,N);
f_e2    = zeros(R+1,N);

x_e     = zeros(R+1,N);

for i = 1:1:N
    for j = 1:1:R
        
        x_e(j,i)     = dX*(j-ceil(R/2));

        f_e1(j,i)     = ((d*0.1)/(4*ri*0.1))*((re*Ie1(j,i))/(4*pi))*(((x_e(j,i))^2+z^2)^(-5/2))*(2*(x_e(j,i))^2-z^2);
        
        f_e2(j,i)     = ((d*0.1)/(4*ri*0.1))*((re*Ie2(j,i))/(4*pi))*(((x_e(j,i))^2+z^2)^(-5/2))*(2*(x_e(j,i))^2-z^2);
        
    end
    
end

figure(1);
plot(x,f_e1(:,101));
title('Activation function for an anodic input signals')
grid on;
xlabel('Distance along the axon [cm]');