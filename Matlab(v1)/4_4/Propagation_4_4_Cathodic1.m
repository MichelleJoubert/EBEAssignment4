% Title:     Hodgkin-Huxley propagation model monopolar Cathodic stimulation
%            Script 1           
% Author:    Bernard Bauermeister (u10367561)
% Goal: The goal of this script is to Simulate the propagation of action 
%potentials when a monopolar electrode (i.e. return electrode is far from 
%stimulation electrode) at a distance z from the nerve fibre provides a 
%single cathodic (negative) monophasic stimulation pulse. The simulation is
%made to repeat for z = 0.05 cm and z = 0.1 cm to show the difference in
%activation area and attenuation of the original stimulation pulse. The
%results are graphed on the same plot to make it easy to compare the two.

clc; 
clear all; 
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Temporal intervals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T       = 13; %Simulation time in milliseconds
N       = 13000; %Number of simulation points (colums in arrays)
dT      = T/N; %Magnitude of time segments 
t       = 0:dT:T; %Time axis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spacial intervals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X       = 2; %Length of the axon in cm
R       = 101; %Amount of segments
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
m0      = 0.08; %Sodium inactivating subunit initial condition
n0      = 0.32;%Potassium subuinit initial condition
Vm0     = Vrest; %Membrane potential initial condition
vm0     = Vm0 - Vrest; %vm initial condition

%State variables arrays
h       = zeros(R,N+1);
m       = zeros(R,N+1);
n       = zeros(R,N+1);
Vm      = zeros(R+1,N);

vm      = zeros(R+1,N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%General constants
Temp    = 6.3; %Temp in Celsius
k       = 3^((0.1*Temp)-0.63); 
ri      = 0.1; %intracellular resistance per cm
re      = 0.3; %extracellular resistance per cm
d       = 0.001; %diameter of axon
L       = dX; %length of axon
M       = ceil(R/2); %Position of stimulus

%Constants that describe electrode position  
ElPos   = M; %Represents half the axon length 
z       = 0.05; %Distance from the axon membrane surface
A       = -3500; %Stimulation amplitude

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
% Hodgkin-Huxley Calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Array creations
f_e     = zeros(R+1,N);
x_e     = zeros(R,N+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stimulating current
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Creates a stimulation current in the form of a block wave with a specific
%amplitude and specific duration
Ie      = zeros(R,N+1); %Create the array for the stimulating current

Ie(:,1:201) = A; %Current stimulation amplitude 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:1:N
    for j = 1:1:R
    
    x_e(j,i)    = dX*(j-ElPos);
    
    %Update vm
    vm(j,i)   = Vm(j,i) - Vrest;  
        
    %Potassium current density
    Gk        = Gk_max*(n(j,i)^4);
    Ik        = Gk*(vm(j,i)-Ek);
    
    %Sodium current density
    Gna       = Gna_max*(m(j,i)^3)*h(j,i);
    Ina       = Gna*(vm(j,i)-Ena);
    
    %Leakage current density
    IL        = Gl*(vm(j,i)-El);
    
    %Alpha and beta 
    an        = (0.1-0.01*vm(j,i))/(exp(1-0.1*vm(j,i))-1);
    bn        = 0.125*exp(-vm(j,i)/80);
    
    am        = (2.5-0.1*vm(j,i))/(exp(2.5-0.1*vm(j,i))-1);
    bm        = 4*exp(-vm(j,i)/18);
    
    ah        = 0.07*exp(-vm(j,i)/20);
    bh        = 1/(exp(3-0.1*vm(j,i))+1);
    
    %changes in variables m,n and h.
    dn        = (an*(1-n(j,i))-bn*n(j,i))*k*dT;
    dm        = (am*(1-m(j,i))-bm*m(j,i))*k*dT;
    dh        = (ah*(1-h(j,i))-bh*h(j,i))*k*dT;
    
    %Membrane Current: Internal Stimulation
    if j == 1 %Provision for the first position
        Im = (d/(4*ri*dX*L))*(-2*vm(j,i) + 2*vm(j+1,i));
    else if j == N %Provision for the last position
        Im = (d/(4*ri*dX*L))*(2*vm(j-1,i) - 2*vm(j,i));
        else 
            Im = (d/(4*ri*dX*L))*(vm(j-1,i) - 2*vm(j,i) + vm(j+1,i));
        end
    end
    
    f_e(j,i) = ((d*dX)/(4*ri*L))*((re*Ie(j,i))/(4*pi))*(((x_e(j,i)^2)+(z^2))^(-5/2))*(2*(x_e(j,i)^2)-(z^2)); 
    
    %Change in membrane voltage calculation.
    dVm      = dT*(Im+f_e(j,i)-Ik-Ina-IL)/Cm;
    
    %Calculation of new state variable values.
    Vm(j,i+1) = Vm(j,i) + dVm;
    vm(j,i+1) = vm(j,i) + dVm; 
    n(j,i+1)  = n(j,i) + dn;
    m(j,i+1)  = m(j,i) + dm;
    h(j,i+1)  = h(j,i) + dh;
    
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);

subplot(2,2,1);
for p = 201:200:N
    
   hold on;
   plot(x,vm(:,p)+(((p-1)/200)*50),'b');
    
end
ylim([-100,3200]);
title('AP propagation with z = 0.05 cm and A = -3.5 mA/cm^2');
xlabel('Distacne along the axon [cm].');

subplot(2,2,3)
plot(x,f_e(:,1));
title('Activation function with z = 0.05 cm and with A = -3.5 mA/cm^2');
xlabel('Distacne along the axon [cm].');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Temporal intervals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T       = 13; %Simulation time in milliseconds
N       = 13000; %Number of simulation points (colums in arrays)
dT      = T/N; %Magnitude of time segments 
t       = 0:dT:T; %Time axis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spacial intervals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X       = 2; %Length of the axon in cm
R       = 101; %Amount of segments
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
m0      = 0.08; %Sodium inactivating subunit initial condition
n0      = 0.32;%Potassium subuinit initial condition
Vm0     = Vrest; %Membrane potential initial condition
vm0     = Vm0 - Vrest; %vm initial condition

%State variables arrays
h       = zeros(R,N+1);
m       = zeros(R,N+1);
n       = zeros(R,N+1);
Vm      = zeros(R+1,N);

vm      = zeros(R+1,N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%General constants
Temp    = 6.3; %Temp in Celsius
k       = 3^((0.1*Temp)-0.63); 
ri      = 0.1; %intracellular resistance per cm
re      = 0.3; %extracellular resistance per cm
d       = 0.001; %diameter of axon
L       = dX; %length of axon
M       = ceil(R/2); %Position of stimulus

%Constants that describe electrode position  
ElPos   = M; %Represents half the axon length 
z       = 0.1; %Distance from the axon membrane surface
A       = -28000; %Stimulation amplitude

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
% Hodgkin-Huxley Calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Array creations
f_e     = zeros(R+1,N);
x_e     = zeros(R,N+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stimulating current
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Creates a stimulation current in the form of a block wave with a specific
%amplitude and specific duration
Ie      = zeros(R,N+1); %Create the array for the stimulating current

Ie(:,1:201) = A; %Current stimulation amplitude 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:1:N
    for j = 1:1:R
    
    x_e(j,i)    = dX*(j-ElPos);
    
    %Update vm
    vm(j,i)   = Vm(j,i) - Vrest;  
        
    %Potassium current density
    Gk        = Gk_max*(n(j,i)^4);
    Ik        = Gk*(vm(j,i)-Ek);
    
    %Sodium current density
    Gna       = Gna_max*(m(j,i)^3)*h(j,i);
    Ina       = Gna*(vm(j,i)-Ena);
    
    %Leakage current density
    IL        = Gl*(vm(j,i)-El);
    
    %Alpha and beta 
    an        = (0.1-0.01*vm(j,i))/(exp(1-0.1*vm(j,i))-1);
    bn        = 0.125*exp(-vm(j,i)/80);
    
    am        = (2.5-0.1*vm(j,i))/(exp(2.5-0.1*vm(j,i))-1);
    bm        = 4*exp(-vm(j,i)/18);
    
    ah        = 0.07*exp(-vm(j,i)/20);
    bh        = 1/(exp(3-0.1*vm(j,i))+1);
    
    %changes in variables m,n and h.
    dn        = (an*(1-n(j,i))-bn*n(j,i))*k*dT;
    dm        = (am*(1-m(j,i))-bm*m(j,i))*k*dT;
    dh        = (ah*(1-h(j,i))-bh*h(j,i))*k*dT;
    
    %Membrane Current: Internal Stimulation
    if j == 1 %Provision for the first position
        Im = (d/(4*ri*dX*L))*(-2*vm(j,i) + 2*vm(j+1,i));
    else if j == N %Provision for the last position
        Im = (d/(4*ri*dX*L))*(2*vm(j-1,i) - 2*vm(j,i));
        else 
            Im = (d/(4*ri*dX*L))*(vm(j-1,i) - 2*vm(j,i) + vm(j+1,i));
        end
    end
    
    f_e(j,i) = ((d*dX)/(4*ri*L))*((re*Ie(j,i))/(4*pi))*(((x_e(j,i)^2)+(z^2))^(-5/2))*(2*(x_e(j,i)^2)-(z^2)); 
    
    %Change in membrane voltage calculation.
    dVm      = dT*(Im+f_e(j,i)-Ik-Ina-IL)/Cm;
    
    %Calculation of new state variable values.
    Vm(j,i+1) = Vm(j,i) + dVm;
    vm(j,i+1) = vm(j,i) + dVm; 
    n(j,i+1)  = n(j,i) + dn;
    m(j,i+1)  = m(j,i) + dm;
    h(j,i+1)  = h(j,i) + dh;
    
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);

subplot(2,2,2);
for p = 201:200:N
    
   hold on;
   plot(x,vm(:,p)+(((p-1)/200)*50),'b');
    
end
ylim([-100,3200]);
title('AP propagation with z = 0.1 cm and A = -28 mA/cm^2');
xlabel('Distacne along the axon [cm].');

subplot(2,2,4)
plot(x,f_e(:,1));
title('Activation function with z = 0.1 cm and with A = -28 mA/cm^2');
xlabel('Distacne along the axon [cm].');

