
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
R       = 3000; %Amount of segments
dX      = X/R; %Magnitudes of spacial segments
x       = 0:dX:X; %Spacial axis

M       = ceil(R/2);

Ie1     = zeros(R,N+1); %Create the array for the stimulating current
Ie2     = zeros(R,N+1); %Create the array for the stimulating current

Ie1(:,1:201) = -3500; %Current stimulation amplitude 
Ie2(:,1:201) = +3500; %Current stimulation amplitude 

d       = 0.001;
z       = 0.05;
re      = 0.3;
ri      = 0.1;
ElPos1  = M*dX - 0.01; 
ElPos2  = M*dX + 0.01;

f_e     = zeros(R+1,N);
f_e1    = zeros(R+1,N);
f_e2    = zeros(R+1,N);

x_e1    = zeros(R+1,N);
x_e2    = zeros(R+1,N);

for i = 1:1:N
    for j = 1:1:R
        
        x_e1(j,i)    = dX*j-ElPos1;
        x_e2(j,i)    = dX*j-ElPos2;

        f_e1(j,i)     = ((d*dX)/(4*ri*dX))*((re*Ie1(j,i))/(4*pi))*(((x_e1(j,i))^2+z^2)^(-5/2))*(2*(x_e1(j,i))^2-z^2);
        
        f_e2(j,i)     = ((d*dX)/(4*ri*dX))*((re*Ie2(j,i))/(4*pi))*(((x_e2(j,i))^2+z^2)^(-5/2))*(2*(x_e2(j,i))^2-z^2);
        
        f_e(j,i)      = f_e1(j,i) + f_e2(j,i);
        
    end
    
end

figure(1);
plot(x,f_e(:,101));
title('Activation function for a bipolar electrode with z = 0.05 and A = 3.5 mA/cm^2 ')
grid on;
xlabel('Distance along the axon [cm].');