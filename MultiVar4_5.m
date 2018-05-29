function [Done] = MultiVar4_5(Distance1,Distance2,I)
%   This function simulates the propagation of an action potential along   
%   nerve fibre of various distances away from the nerve fibre stimulated by    
%	a bipolar stimulation pulse.
%
%   MultiVar4_5(Distance1,Distance2,I) function simulates the Hodgkin-Huxley model 
%   for the squid giant axon constants and user specified values of 
%   stimulation distance from the fibre in cm for specified input stimulus. 
%   As output it plots voltage(membrane potential)time series for each of  
%   the various fibre diameters. 
%   Distance1 is the first electrode's distance away from the axon (in cm).
%   Distance2 is the second electrode's distance away from the axon (in cm).
%   I is the applied current or external stimulation (in µA)
%   
%   Example:
%   MultiVar4_5(0.05,0.05,3) for symmetric bipolar stimulation (for Question 4.5)
%   MultiVar4_5(0.05,0.055,3) for asymetric bipolar stimulation (for Question 4.5)
%
%%	Simulation timing variables
    t = 0;
    loop = 0;
    dt = 0.01;
    tspan = 200;
    [t,loop] = FindT(tspan,dt);
    
%%	Simulation spatial variables
    x = 0;
    xloop = 0;
    dx = 0.02;
    xspan = 2;
	[x,xloop] = FindX(xspan,dx);

	Done = 'Done';
    
%%  Running through the various input diameters 
    [data,t,f,x] = HHsim(Distance1,Distance2,I,t,x,loop,xloop,dx);
    
%%  Plots of membrane potential for various input diameters   
	figure
    for n = 2:1:xloop-1
        plot(t,data(n,:)+(n*15),'b');
        hold on
    end
    xlabel('Time (msec)');
    ylabel('Membrane Potential (mV)');
    label=strcat('AP propagation with electrodes a distance ',{' '},num2str(Distance1),' and',{' '},num2str(Distance2),'cm from axon');
    title(label);
        
    figure
    plot(x,f(:,190),'b');
    xlabel('Fibre length (cm)');
    ylabel('Activating function (mV)');
    label=strcat('Activation function with electrodes a distance ',{' '},num2str(Distance1),' and',{' '},num2str(Distance2),'cm from axon');
    title(label);
end

%%	Creating simulation timeframe
function [t,loop] = FindT(tspan,dt)
    loop  = ceil(tspan/dt);
    t = (0:loop)*dt;
end

%%	Creating simulation fibre length
function [x,xloop] = FindX(xspan,dx)
    xloop  = ceil(xspan/dx);
    x = (0:xloop)*dx;
end

%%  Hodgkin-Huxley model
function [v,t,f,x] = HHsim(z1,z2,defI,t,x,loop,xloop,dx)    
    defTemp = 6.3;      % environmental temperature (in deg Celsius)
    dt = 0.001;         % time steps for simulation
    Diameter = 0.01;    % Diameter of fibre (in cm)0.012

%%	Constants and intial values for squid giant axon  
    gNa = 120;                  % Conductance of sodium channels (in m.mho/cm^2)
	vNa= 115;                   % Ionic potential for sodium channels (in mV)
	gK = 36;                    % Conductance of potassium channels (in m.mho/cm^2)
	vK= -12;                    % Ionic potential for potassium channels (in mV)
	gL= 0.3;                    % Conductance of leakage channels (in m.mho/cm^2)
	vL= 10.6;                   % Ionic potential for leakage channels (in mV)
	Cm = 1;                     % Capacitance of membrane (in µF/cm^2)
    L = dx;                     % unmyelinated fibre
    ri = 0.1;                   % specific resistance of axoplasm (in kOhm.cm)
    re = 0.3;                   % specific extracellular resistance (in kOhm.cm)
    m0 = 0.05;                  % Initial value of gating variable m
    n0 = 0.32;                  % Initial value of gating variable n
    h0 = 0.6;                   % Initial value of gating variable h
    Vrest = -70;                % Membrane potential at rest (in mV)
    V0 = Vrest;                 % Initial value of actual membrane potential (in mV)
    v0 = V0 - Vrest;            % Initial value of normalised membrane potential (in mV)    
    
%%	Initializing variable vectors  
    V = zeros(xloop+1,loop+1);      % Membrane potential
    v = zeros(xloop+1,loop+1);      % Normalised membrane potential
    Im = zeros(xloop,loop+1);       % Transmembrane current
    Iion = zeros(xloop,loop+1);     % Ionic currents
    xe1 = zeros(xloop,loop+1);      % Distance along axon away from electrode 1
    xe2 = zeros(xloop,loop+1);      % Distance along axon away from electrode 2
    f1 = zeros(xloop+1,loop+1);     % Activating function
    f2 = zeros(xloop+1,loop+1);     % Activating function
    f = zeros(xloop+1,loop+1);      % Activating function
	m = zeros(xloop,loop+1);        % Gating variable m
	h = zeros(xloop,loop+1);        % Gating variable h
	n = zeros(xloop,loop+1);        % Gating variable n
    dm = zeros(xloop,loop+1);       % Change in Gating variable m
	dh = zeros(xloop,loop+1);       % Change in Gating variable h
	dn = zeros(xloop,loop+1);       % Change in Gating variable h
	Ispan1 = zeros(xloop,loop+1);   % Stimulating current
    Ispan2 = zeros(xloop,loop+1);   % Stimulating current
    alphaN = zeros(xloop,loop+1);   % Gating variable rate constant
    betaN = zeros(xloop,loop+1);    % Gating variable rate constant
    alphaM = zeros(xloop,loop+1);   % Gating variable rate constant
    betaM = zeros(xloop,loop+1);    % Gating variable rate constant
    alphaH = zeros(xloop,loop+1);   % Gating variable rate constant
    betaH = zeros(xloop,loop+1);    % Gating variable rate constant

%%	Ispan is the applied current vector to hold all instances of the external 
    p = xloop/2;
    Ispan2(p+1:p+11,100:150) = +defI;
    Ispan1(p-11:p-1,100:150) = -defI;
    
%%	Phi is the temperature adjusting factor to be applied to the gating variables
    phi = 3^((0.1*defTemp)-0.63);
    
%%	Set initial values for the membrane potentials and gating variables   
    for s = 1:xloop
        V(s,1) = V0;
        v(s,1) = v0;
        n(s,1) = n0;
        m(s,1) = m0;
        h(s,1) = h0;
    end
    
%%	Solver: Euler method
    for i=1:1:loop-1
        for s=2:1:xloop-1
            electrodePos1 = p*dx-0.01;
            xe1(s,i) = dx*s-electrodePos1;
            electrodePos2 = p*dx+0.01;
            xe2(s,i) = dx*s-electrodePos2;
        end
    end
 
    for i=1:1:loop-1
        for s=2:1:xloop-1  
            v(s,i) = V(s,i) - Vrest;
            alphaN(s,i) = (0.1-0.01*v(s,i))/(exp(1-0.1*v(s,i))-1);
            betaN(s,i) = 0.125*exp(-v(s,i)/80);
            alphaM(s,i) = (2.5-0.1*v(s,i))/(exp(2.5-0.1*v(s,i))-1);
            betaM(s,i) = 4*exp(-v(s,i)/18);
            alphaH(s,i) = 0.07*exp(-v(s,i)/20);
            betaH(s,i) = 1/(exp(3-0.1*v(s,i))+1);

            Iion(s,i) = gNa*(m(s,i)^3)*h(s,i)*(v(s,i)-vNa) + gK*(n(s,i)^4)*(v(s,i)-vK) + gL*(v(s,i)-vL);
            Im(s,i) = (Diameter/(4*ri*dx*L))*(v(s-1,i)-2*v(s,i)+v(s+1,i));
            f1(s,i) = (Diameter*dx/(4*ri*L))*((re*Ispan1(s,i))/(4*pi))*(((xe1(s,i)^2)+(z1^2))^(-5/2))*(2*(xe1(s,i)^2)-(z1^2));
            f2(s,i) = (Diameter*dx/(4*ri*L))*((re*Ispan2(s,i))/(4*pi))*(((xe2(s,i)^2)+(z2^2))^(-5/2))*(2*(xe2(s,i)^2)-(z2^2));
            dn(s,i) = phi*dt*(alphaN(s,i)*(1-n(s,i)) - betaN(s,i)*n(s,i));
            dm(s,i) = phi*dt*(alphaM(s,i)*(1-m(s,i)) - betaM(s,i)*m(s,i));
            dh(s,i) = phi*dt*(alphaH(s,i)*(1-h(s,i)) - betaH(s,i)*h(s,i));
            
            f(s,i+1) =  f(s,i) + f1(s,i)+ f2(s,i);
            V(s,i+1) = V(s,i) + dt*(Im(s,i) + f(s,i) - Iion(s,i))/Cm; 
            v(s,i+1) = v(s,i) + dt*(Im(s,i) + f(s,i) - Iion(s,i))/Cm;
            n(s,i+1) = n(s,i) + dn(s,i);
            m(s,i+1) = m(s,i) + dm(s,i);
            h(s,i+1) = h(s,i) + dh(s,i);
        end
    end
end


