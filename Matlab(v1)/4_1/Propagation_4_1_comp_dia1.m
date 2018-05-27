% Title:     Hodgkin-Huxley propagation model (effect of diameter) 
% Author:    Bernard Bauermeister (u10367561)
% Goal: The goal of this script is to create an action potential propagation
% model to test the differenct in propagaion patterns for two axon diameters

clc; 
clear all; 
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Temporal intervals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T       = 20; %Simulation time in milliseconds
N       = 20000; %Number of simulation points (colums in arrays)
dT      = T/N; %Magnitude of time segments 
t       = 0:dT:T; %Time axis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spacial intervals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X       = 10; %Simulation time in milliseconds
R       = 25; %Amount of segments
dX      = X/R; %Magnitudes of spacial segments
x       = 0:dX:X; %Spacial axis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constant Variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%For the potassium channel
Ek      = -12; %Nernst potential of potassium
Gk_max  = 36; %Maximuim conductance of potassium (when all channels are open)

%for the sodium channel
Ena     = 115; %Nernst potential of sodium
Gna_max = 120; %maximuim conductance of sodium (when all channels are open)

%For the leakage channel
El      = 10.6; %Nernst potential of the leakage components
Gl      = 0.3; %Corresponding conductance of the leakage components

%For the capacitance current
Cm      = 1; %The capacitance per unit area

%State variable initial values
Vrest   = -70; %Cell membrane resting potential
h0      = 0.6; %Sodium activiating subunit initial condition
m0      = 0.05; %Sodium inactivating subunit initial condition
n0      = 0.32;%Potassium subuinit initial condition
Vm0     = Vrest; %Membrane potential initial condition
vm0     = Vm0 - Vrest; %vm initial condition

%State variables arrays
h       = zeros(R,N+1);
m       = zeros(R,N+1);
n       = zeros(R,N+1);
Vm      = zeros(R+1,N);

vm      = zeros(R+1,N);

%General constants
Temp    = 6.3; %Temp in Celsius
k       = 3^((0.1*Temp)-0.63); 

%intracellular resistance per cm
ri      = 0.1; 

%extracellular resistance per cm
re      = 0.3; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%diameter of axon
d       = 0.005;

%length of axon simulation
L       = dX; %Unmylinated case

%Position of stimulus
M       = ceil(R/2);

%Set initial conditions for next amplitude
for j = 1:1:R
    
%Set initial conditions for next amplitude
%State variables arrays
h(j,1)    = h0;
m(j,1)    = m0;
n(j,1)    = n0;
Vm(j,1)   = Vm0;

%Other arrays
vm(j,1)   = vm0;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stimulating current
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Creates a stimulation current in the form of a block wave with a specific
%amplitude and specific duration
Is      = zeros(R,N+1); %Create the array for the stimulating current

Is(M,2000:2400) = 150;  %Current Stimulation starting in the middle of the axon

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hodgkin-Huxley Calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Array creations
Gk      = zeros(R,N+1);
Gna     = zeros(R,N+1);

Ik      = zeros(R,N+1);
Ina     = zeros(R,N+1);
IL      = zeros(R,N+1);

ah      = zeros(R,N+1);
bh      = zeros(R,N+1);

am      = zeros(R,N+1);
bm      = zeros(R,N+1);

an      = zeros(R,N+1);
bn      = zeros(R,N+1);

dn      = zeros(R,N+1);
dm      = zeros(R,N+1);
dh      = zeros(R,N+1);

dVm     = zeros(R,N+1);

Im      = zeros(R,N+1);

%calculates for the distances after the stimulus
for i = 1:1:N
    for j = 2:1:R  
    
    %Update vm
    vm(j,i)   = Vm(j,i) - Vrest;  
        
    %Membrane Current: Internal Stimulation
    Im(j,i) = (d/(4*ri*dX*L))*(vm(j-1,i) - 2*vm(j,i) + vm(j+1,i))+ Is(j,i);%(Is(j,i)/(pi*d*L));
    
    %Potassium current density
    Gk(j,i)   = Gk_max*(n(j,i)^4);
    Ik(j,i)   = Gk(j,i)*(vm(j,i)-Ek);
    
    %Sodium current density
    Gna(j,i)  = Gna_max*(m(j,i)^3)*h(j,i);
    Ina(j,i)  = Gna(j,i)*(vm(j,i)-Ena);
    
    %Leakage current density
    IL(j,i)   = Gl*(vm(j,i)-El);
    
    %Alpha and beta 
    an(j,i)   = (0.1-0.01*vm(j,i))/(exp(1-0.1*vm(j,i))-1);
    bn(j,i)   = 0.125*exp(-vm(j,i)/80);
    
    am(j,i)   = (2.5-0.1*vm(j,i))/(exp(2.5-0.1*vm(j,i))-1);
    bm(j,i)   = 4*exp(-vm(j,i)/18);
    
    ah(j,i)   = 0.07*exp(-vm(j,i)/20);
    bh(j,i)   = 1/(exp(3-0.1*vm(j,i))+1);
    
    %changes in variables m,n and h.
    dn(j,i)   = (an(j,i)*(1-n(j,i))-bn(j,i)*n(j,i))*k*dT;
    dm(j,i)   = (am(j,i)*(1-m(j,i))-bm(j,i)*m(j,i))*k*dT;
    dh(j,i)   = (ah(j,i)*(1-h(j,i))-bh(j,i)*h(j,i))*k*dT;
    
    %Change in membrane voltage calculation.
    dVm(j,i)  = dT*(Im(j,i)-Ik(j,i)-Ina(j,i)-IL(j,i))/Cm;
    
    %Calculation of new state variable values.
    Vm(j,i+1) = Vm(j,i) + dVm(j,i);
    vm(j,i+1) = vm(j,i) + dVm(j,i);
    n(j,i+1)  = n(j,i) + dn(j,i);
    m(j,i+1)  = m(j,i) + dm(j,i);
    h(j,i+1)  = h(j,i) + dh(j,i);
    
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
for p = 2:1:R
    
    plot(t,vm(p,:)+(p*120),'b');
    hold on;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code is repeated to calculate and compare the propagation for a different
%diameter. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constant Variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%For the potassium channel
Ek      = -12; %Nernst potential of potassium
Gk_max  = 36; %Maximuim conductance of potassium (when all channels are open)

%for the sodium channel
Ena     = 115; %Nernst potential of sodium
Gna_max = 120; %maximuim conductance of sodium (when all channels are open)

%For the leakage channel
El      = 10.6; %Nernst potential of the leakage components
Gl      = 0.3; %Corresponding conductance of the leakage components

%For the capacitance current
Cm      = 1; %The capacitance per unit area

%State variable initial values
Vrest   = -70; %Cell membrane resting potential
h0      = 0.6; %Sodium activiating subunit initial condition
m0      = 0.05; %Sodium inactivating subunit initial condition
n0      = 0.32;%Potassium subuinit initial condition
Vm0     = Vrest; %Membrane potential initial condition
vm0     = Vm0 - Vrest; %vm initial condition

%State variables arrays
h       = zeros(R,N+1);
m       = zeros(R,N+1);
n       = zeros(R,N+1);
Vm      = zeros(R+1,N);

vm      = zeros(R+1,N);

%General constants
Temp    = 6.3; %Temp in Celsius
k       = 3^((0.1*Temp)-0.63); 

%intracellular resistance per cm
ri      = 0.1; 

%extracellular resistance per cm
re      = 0.3; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%diameter of axon
d       = 0.01;

%length of axon simulation
L       = dX; %Unmylinated case

%Position of stimulus
M       = ceil(R/2);

%Set initial conditions for next amplitude
for j = 1:1:R
    
%Set initial conditions for next amplitude
%State variables arrays
h(j,1)    = h0;
m(j,1)    = m0;
n(j,1)    = n0;
Vm(j,1)   = Vm0;

%Other arrays
vm(j,1)   = vm0;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stimulating current
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Creates a stimulation current in the form of a block wave with a specific
%amplitude and specific duration
Is      = zeros(R,N+1); %Create the array for the stimulating current

Is(M,2000:4000) = 0.5;  %Current Stimulation starting in the middle of the axon

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hodgkin-Huxley Calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Array creations
Gk      = zeros(R,N+1);
Gna     = zeros(R,N+1);

Ik      = zeros(R,N+1);
Ina     = zeros(R,N+1);
IL      = zeros(R,N+1);

ah      = zeros(R,N+1);
bh      = zeros(R,N+1);

am      = zeros(R,N+1);
bm      = zeros(R,N+1);

an      = zeros(R,N+1);
bn      = zeros(R,N+1);

dn      = zeros(R,N+1);
dm      = zeros(R,N+1);
dh      = zeros(R,N+1);

dVm     = zeros(R,N+1);

Im      = zeros(R,N+1);

%calculates for the distances after the stimulus
for i = 1:1:N
    for j = 2:1:R
    
    %Update vm
    vm(j,i)   = Vm(j,i) - Vrest;  
        
    %Membrane Current: Internal Stimulation
    Im(j,i) = (d/(4*ri*dX*L))*(vm(j-1,i) - 2*vm(j,i) + vm(j+1,i))+ (Is(j,i)/(pi*d*L));
    
    %Potassium current density
    Gk(j,i)   = Gk_max*(n(j,i)^4);
    Ik(j,i)   = Gk(j,i)*(vm(j,i)-Ek);
    
    %Sodium current density
    Gna(j,i)  = Gna_max*(m(j,i)^3)*h(j,i);
    Ina(j,i)  = Gna(j,i)*(vm(j,i)-Ena);
    
    %Leakage current density
    IL(j,i)   = Gl*(vm(j,i)-El);
    
    %Alpha and beta 
    an(j,i)   = (0.1-0.01*vm(j,i))/(exp(1-0.1*vm(j,i))-1);
    bn(j,i)   = 0.125*exp(-vm(j,i)/80);
    
    am(j,i)   = (2.5-0.1*vm(j,i))/(exp(2.5-0.1*vm(j,i))-1);
    bm(j,i)   = 4*exp(-vm(j,i)/18);
    
    ah(j,i)   = 0.07*exp(-vm(j,i)/20);
    bh(j,i)   = 1/(exp(3-0.1*vm(j,i))+1);
    
    %changes in variables m,n and h.
    dn(j,i)   = (an(j,i)*(1-n(j,i))-bn(j,i)*n(j,i))*k*dT;
    dm(j,i)   = (am(j,i)*(1-m(j,i))-bm(j,i)*m(j,i))*k*dT;
    dh(j,i)   = (ah(j,i)*(1-h(j,i))-bh(j,i)*h(j,i))*k*dT;
    
    %Change in membrane voltage calculation.
    dVm(j,i)  = dT*(Im(j,i)-Ik(j,i)-Ina(j,i)-IL(j,i))/Cm;
    
    %Calculation of new state variable values.
    Vm(j,i+1) = Vm(j,i) + dVm(j,i);
    vm(j,i+1) = vm(j,i) + dVm(j,i);
    n(j,i+1)  = n(j,i) + dn(j,i);
    m(j,i+1)  = m(j,i) + dm(j,i);
    h(j,i+1)  = h(j,i) + dh(j,i);
    
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
for p = 2:1:R
    
    plot(t,vm(p,:)+(p*120),'r');
    hold on;
    
end
title('Difference in propagation patterns.');
xlabel('Simulation time [ms]');



