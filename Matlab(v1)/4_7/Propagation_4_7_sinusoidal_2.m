% Title:     Hodgkin-Huxley propagation model sinusoidal stimulation
% Author:    Bernard Bauermeister (u10367561)
%            Script 2
% Goal: The goal of this script is to Simulate the propagation of action 
%potentials when a sinusoidal stimulation signal is used to achieve
%entertainment on the nerve axon.

clc; 
clear all; 
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Temporal intervals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T       = 100; %Simulation time in milliseconds
N       = 50000; %Number of simulation points (colums in arrays)
dT      = T/N; %Magnitude of time segments 
t       = 0:dT:T; %Time axis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spacial intervals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X       = 10; %Simulation time in milliseconds
R       = 100; %Amount of segments 
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
Temp    = 10; %Temp in Celsius
k       = 3^((0.1*Temp)-0.63); 

%intracellular resistance per cm
ri      = 0.1; 

%extracellular resistance per cm
re      = 0.3; 

%diameter of axon
d       = 0.01;

L       = dX; %Unmylinated

M       = ceil(R/2);


%Sinusoidal variables
A       = 100; %Amplitude of stimulation
f1      = 140; %In Hz

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

Is      = zeros(R,N+1); %Create the array for the stimulating current

Is(M,:) = A*sin(2*pi*f1*(10^(-3))*t);  %Current Stimulation starting in the middle of the axon 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hodgkin-Huxley Calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:1:N
    for j = 2:1:R
    
    %Update vm
    vm(j,i)   = Vm(j,i) - Vrest;  
        
    %Membrane Current: Internal Stimulation
    Im = (d/(4*ri*dX*L))*(vm(j-1,i) - 2*vm(j,i) + vm(j+1,i))+ (Is(j,i));
    
    %Potassium current density
    Gk   = Gk_max*(n(j,i)^4);
    Ik   = Gk*(vm(j,i)-Ek);
    
    %Sodium current density
    Gna  = Gna_max*m(j,i)^3*h(j,i);
    Ina  = Gna*(vm(j,i)-Ena);
    
    %Leakage current density
    IL   = Gl*(vm(j,i)-El);
    
    %Alpha and beta 
    an   = (0.1-0.01*vm(j,i))/(exp(1-0.1*vm(j,i))-1);
    bn   = 0.125*exp(-vm(j,i)/80);
    
    am   = (2.5-0.1*vm(j,i))/(exp(2.5-0.1*vm(j,i))-1);
    bm   = 4*exp(-vm(j,i)/18);
    
    ah   = 0.07*exp(-vm(j,i)/20);
    bh   = 1/(exp(3-0.1*vm(j,i))+1);
    
    %changes in variables m,n and h.
    dn   = (an*(1-n(j,i))-bn*n(j,i))*k*dT;
    dm   = (am*(1-m(j,i))-bm*m(j,i))*k*dT;
    dh   = (ah*(1-h(j,i))-bh*h(j,i))*k*dT;
    
    %Change in membrane voltage calculation.
    dVm  = dT*(Im-Ik-Ina-IL)/Cm;
    
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
subplot(2,1,1);
plot(t,vm(M,:),t,Is(M,:),'--');
title('Sinusoidal input with A = 100 uA/cm^2  and T = 10^{\circ}C');

subplot(2,1,2);
for p = 2:1:R
    
    plot(t,vm(p,:)+(p*20),'b');
    hold on;
    
end
ylim([0, 2100]);
title('AP propagation with a sinusoidal input where A = 100 uA/cm^2  and T = 10^{\circ}C.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Temporal intervals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T       = 100; %Simulation time in milliseconds
N       = 50000; %Number of simulation points (colums in arrays)
dT      = T/N; %Magnitude of time segments 
t       = 0:dT:T; %Time axis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spacial intervals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X       = 10; %Simulation time in milliseconds
R       = 100; %Amount of segments 
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
Temp    = 11; %Temp in Celsius
k       = 3^((0.1*Temp)-0.63); 

%intracellular resistance per cm
ri      = 0.1; 

%extracellular resistance per cm
re      = 0.3; 

%diameter of axon
d       = 0.01;

L       = dX; %Unmylinated

M       = ceil(R/2);


%Sinusoidal variables
A       = 65; %Amplitude of stimulation
f1      = 140; %In Hz

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

Is      = zeros(R,N+1); %Create the array for the stimulating current

Is(M,:) = A*sin(2*pi*f1*(10^(-3))*t);  %Current Stimulation starting in the middle of the axon (33)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hodgkin-Huxley Calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:1:N
    for j = 2:1:R
    
    %Update vm
    vm(j,i)   = Vm(j,i) - Vrest;  
        
    %Membrane Current: Internal Stimulation
    Im = (d/(4*ri*dX*L))*(vm(j-1,i) - 2*vm(j,i) + vm(j+1,i))+ (Is(j,i));
    
    %Potassium current density
    Gk   = Gk_max*(n(j,i)^4);
    Ik   = Gk*(vm(j,i)-Ek);
    
    %Sodium current density
    Gna  = Gna_max*m(j,i)^3*h(j,i);
    Ina  = Gna*(vm(j,i)-Ena);
    
    %Leakage current density
    IL   = Gl*(vm(j,i)-El);
    
    %Alpha and beta 
    an   = (0.1-0.01*vm(j,i))/(exp(1-0.1*vm(j,i))-1);
    bn   = 0.125*exp(-vm(j,i)/80);
    
    am   = (2.5-0.1*vm(j,i))/(exp(2.5-0.1*vm(j,i))-1);
    bm   = 4*exp(-vm(j,i)/18);
    
    ah   = 0.07*exp(-vm(j,i)/20);
    bh   = 1/(exp(3-0.1*vm(j,i))+1);
    
    %changes in variables m,n and h.
    dn   = (an*(1-n(j,i))-bn*n(j,i))*k*dT;
    dm   = (am*(1-m(j,i))-bm*m(j,i))*k*dT;
    dh   = (ah*(1-h(j,i))-bh*h(j,i))*k*dT;
    
    %Change in membrane voltage calculation.
    dVm  = dT*(Im-Ik-Ina-IL)/Cm;
    
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
figure(2);
subplot(2,1,1);
plot(t,vm(M,:),t,Is(M,:),'--');
title('Sinusoidal input with A = 65 uA/cm^2 and T = 11^{\circ}C');

subplot(2,1,2);
for p = 2:1:R
    
    plot(t,vm(p,:)+(p*20),'b');
    hold on;
    
end
ylim([0,2100]);
title('AP propagation with a sinusoidal input where A = 65 uA/cm^2  and T = 11^{\circ}C.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Temporal intervals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T       = 100; %Simulation time in milliseconds
N       = 50000; %Number of simulation points (colums in arrays)
dT      = T/N; %Magnitude of time segments 
t       = 0:dT:T; %Time axis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spacial intervals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X       = 10; %Simulation time in milliseconds
R       = 100; %Amount of segments 
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
Temp    = 10; %Temp in Celsius
k       = 3^((0.1*Temp)-0.63); 

%intracellular resistance per cm
ri      = 0.1; 

%extracellular resistance per cm
re      = 0.3; 

%diameter of axon
d       = 0.01;

L       = dX; %Unmylinated

M       = ceil(R/2);


%Sinusoidal variables
A       = 65; %Amplitude of stimulation
f1      = 140; %In Hz

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

Is      = zeros(R,N+1); %Create the array for the stimulating current

Is(M,:) = A*sin(2*pi*f1*(10^(-3))*t);  %Current Stimulation starting in the middle of the axon (33)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hodgkin-Huxley Calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:1:N
    for j = 2:1:R
    
    %Update vm
    vm(j,i)   = Vm(j,i) - Vrest;  
        
    %Membrane Current: Internal Stimulation
    Im = (d/(4*ri*dX*L))*(vm(j-1,i) - 2*vm(j,i) + vm(j+1,i))+ (Is(j,i));
    
    %Potassium current density
    Gk   = Gk_max*(n(j,i)^4);
    Ik   = Gk*(vm(j,i)-Ek);
    
    %Sodium current density
    Gna  = Gna_max*m(j,i)^3*h(j,i);
    Ina  = Gna*(vm(j,i)-Ena);
    
    %Leakage current density
    IL   = Gl*(vm(j,i)-El);
    
    %Alpha and beta 
    an   = (0.1-0.01*vm(j,i))/(exp(1-0.1*vm(j,i))-1);
    bn   = 0.125*exp(-vm(j,i)/80);
    
    am   = (2.5-0.1*vm(j,i))/(exp(2.5-0.1*vm(j,i))-1);
    bm   = 4*exp(-vm(j,i)/18);
    
    ah   = 0.07*exp(-vm(j,i)/20);
    bh   = 1/(exp(3-0.1*vm(j,i))+1);
    
    %changes in variables m,n and h.
    dn   = (an*(1-n(j,i))-bn*n(j,i))*k*dT;
    dm   = (am*(1-m(j,i))-bm*m(j,i))*k*dT;
    dh   = (ah*(1-h(j,i))-bh*h(j,i))*k*dT;
    
    %Change in membrane voltage calculation.
    dVm  = dT*(Im-Ik-Ina-IL)/Cm;
    
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
figure(3);
subplot(2,1,1);
plot(t,vm(M,:),t,Is(M,:),'--');
title('Sinusoidal input with A = 65 uA/cm^3 and T = 10^{\circ}C');

subplot(2,1,2);
for p = 2:1:R
    
    plot(t,vm(p,:)+(p*20),'b');
    hold on;
    
end
ylim([0,2100]);
title('AP propagation with a sinusoidal input where A = 65 uA/cm^2  and T = 10^{\circ}C.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%