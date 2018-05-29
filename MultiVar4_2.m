function [Done] = MultiVar4_2(Diameter,I,sType)
%   This function simulates the propagation of an action potential across a synapse. 
% 
%   MultiVar4_2(Diameter,I,sType) function is a simple synapse model that uses  
%   the Hodgkin-Huxley modelfor the squid giant axon constants and user specified  
%   values of thefibre diameter in cm. As output it plots voltage (membrane   
%   potential) time series for each of the various fibre diameters. 
%   Diameter is the fibre diameter (in cm). 
%   I is the applied current or external stimulation (in µA)
%   Selecting "sType" = 1 applies inhibitory synapse characteristics
%   Selecting "sType" = 2 applies excitatory synapse characteristics
%
%   Example:
%   MultiVar4_2(0.01,1,1) for inhibitory synapse (for Question 4.2)
%   MultiVar4_2(0.01,1,2) for excitatory synapse (for Question 4.2)
%
%%	Simulation timing variables
    t = 0;
    loop = 0;
    dt = 0.01;
    tspan = 400;
    [t,loop] = FindT(tspan,dt);
    
%%	Simulation spatial variables
    x = 0;
    xloop = 0;
    dx = 0.05;
    xspan = 10;
	[x,xloop] = FindX(xspan,dx);

	Done = 'Done';
    
%%  Running through the various input diameters 
        [data,t,synapsePos] = HHsim(Diameter,I,sType,t,x,loop,xloop,dx);
    
%%  Plots of membrane potential for various input diameters   
        figure
        for n = 2:7:xloop-1
            hold on
            if n == synapsePos 
                plot(t,data(n,:)+(n*10),'r');
            else
                plot(t,data(n,:)+(n*10),'b');
            end   
        end
        xlabel('Time (msec)');
        ylabel('Membrane Potential (mV)');
        if sType == 1
            label=strcat('Action potential propagation across an inhibitory synapse');
        elseif sType == 2
            label=strcat('Action potential propagation across an excitatory synapse');
        end
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
function [v,t,synapsePos] = HHsim(defDiameter,defI,defsType,t,x,loop,xloop,dx)    
    defTemp = 6.3;      % environmental temperature (in deg Celsius)
    dt = 0.001;         % time steps
    
%%	Constants and intial values for squid giant axon  
    gNa = 120;                  % Conductance of sodium channels (in m.mho/cm^2)
	vNa= 115;                   % Ionic potential for sodium channels (in mV)
	gK = 36;                    % Conductance of potassium channels (in m.mho/cm^2)
	vK= -12;                    % Ionic potential for potassium channels (in mV)
	gL= 0.3;                    % Conductance of leakage channels (in m.mho/cm^2)
	vL= 10.6;                   % Ionic potential for leakage channels (in mV)
	Cm = 1;                     % Capacitance of membrane (in micro.F/cm^2)
    L = dx;                     % unmyelinated fibre
    ri = 0.1;                   % specific resistance of axoplasm (in kOhm.cm)
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
	m = zeros(xloop,loop+1);        % Gating variable m
	h = zeros(xloop,loop+1);        % Gating variable h
	n = zeros(xloop,loop+1);        % Gating variable n
    dm = zeros(xloop,loop+1);       % Change in Gating variable m
	dh = zeros(xloop,loop+1);       % Change in Gating variable h
	dn = zeros(xloop,loop+1);       % Change in Gating variable n
	Ispan = zeros(xloop,loop+1);    % Stimulating current
    Isynapse = zeros(xloop,loop);   % Synaptic current
    gSyn = zeros(xloop,loop);       % Synaptic Conductance
    alphaN = zeros(xloop,loop+1);   % Gating variable rate constant
    betaN = zeros(xloop,loop+1);    % Gating variable rate constant
    alphaM = zeros(xloop,loop+1);   % Gating variable rate constant
    betaM = zeros(xloop,loop+1);	% Gating variable rate constant
    alphaH = zeros(xloop,loop+1);   % Gating variable rate constant
    betaH = zeros(xloop,loop+1);	% Gating variable rate constant
    
%%	Ispan is the applied current vector to hold all instances of the external 
    Ispan(2,3000:4000) = defI;
    
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
    for i=1:loop-1
        for s=2:xloop-1
            v(s,i) = V(s,i) - Vrest;
            alphaN(s,i) = (0.1-0.01*v(s,i))/(exp(1-0.1*v(s,i))-1);
            betaN(s,i) = 0.125*exp(-v(s,i)/80);
            alphaM(s,i) = (2.5-0.1*v(s,i))/(exp(2.5-0.1*v(s,i))-1);
            betaM(s,i) = 4*exp(-v(s,i)/18);
            alphaH(s,i) = 0.07*exp(-v(s,i)/20);
            betaH(s,i) = 1/(exp(3-0.1*v(s,i))+1);
            
            synapsePos = (xloop/2);     % Position (spatial) of synapse
            t_0 = 17700;
            if s == synapsePos && i >= t_0      
                t0 = 10*t_0*dt;         % presynaptic action potential arrival time (from running MultiVar4_1 for fibre of half the length of total xspan)
                if defsType == 1        % Inhibitory synapse
                    Esyn = -75;         % Synatic reversal potential for inhibitory synapses (in mV) (-75mV)
                elseif defsType == 2    % Excitatory synapse
                    Esyn = 0;           % Synatic reversal potential for Excitatory synapses (in mV) (0mV)
                end
                gSynMax = 40;           % maximum Synatic conductance for general synapses (in pS) (40pS)
                tau = 5;                % Time constant for general synapses (in msec) (5msec)
                gSyn(s,i) = gSynMax*((t(1,i)-t0)/tau)*exp((t0-t(1,i))/tau);
                Isynapse(s,i) = gSyn(s,i)*(v(s,i) - Esyn);    
            end

            Iion(s,i) = gNa*(m(s,i)^3)*h(s,i)*(v(s,i)-vNa) + gK*(n(s,i)^4)*(v(s,i)-vK) + gL*(v(s,i)-vL);
            Im(s,i) = (defDiameter/(4*ri*dx*L))*(v(s-1,i)-2*v(s,i)+v(s+1,i)) + (Ispan(s,i)/(pi*defDiameter*L));
            dn(s,i) = phi*dt*(alphaN(s,i)*(1-n(s,i)) - betaN(s,i)*n(s,i));
            dm(s,i) = phi*dt*(alphaM(s,i)*(1-m(s,i)) - betaM(s,i)*m(s,i));
            dh(s,i) = phi*dt*(alphaH(s,i)*(1-h(s,i)) - betaH(s,i)*h(s,i));
           
            V(s,i+1) = V(s,i) + dt*(Im(s,i) - Iion(s,i) - Isynapse(s,i))/Cm;
            v(s,i+1) = v(s,i) + dt*(Im(s,i) - Iion(s,i))/Cm; 
            n(s,i+1) = n(s,i) + dn(s,i);
            m(s,i+1) = m(s,i) + dm(s,i);
            h(s,i+1) = h(s,i) + dh(s,i);
        end
    end
end


