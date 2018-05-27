function [Done] = MultiVar4_7(Diameter,I,PType)
%   This function simulates the propagation of an action potential along   
%   nerve fibre for various waveform stimulus.
%
%   MultiVar4_7(Diameter,I,PType) function simulates the Hodgkin-Huxley model 
%   for the squid giant axon constants and user specified waveforms. As output    
%   it plots voltage (membrane potential) time series for each of the various  
%   fibre diameters.
%   Diameter is the fibre diameter (in cm). This variable is a list of
%   values to allow for investigation of multiple diameters at once.
%   I is the applied current or external stimulation (in µA)
%   Selecting "PType" = 1 applies a sine wave stimulus
%   Selecting "PType" = 2 applies a square wave stimulus
%   Selecting "PType" = 3 applies a ramp function stimulus 
%
%   Example:
%   MultiVar4_7([0.01],0.4,1) for sine wave with freq 12.5Hz, max of 13.7hz for current
%   length of fibre for timespan 800msec
%   MultiVar4_7([0.01],0.4,2) for square wave
%   MultiVar4_7([0.01],1,3) for ramp function
%
%%	Simulation timing variables
    t = 0;
    loop = 0;
    dt = 0.01;
    tspan = 800;
    [t,loop] = FindT(tspan,dt);
    
%%	Simulation spatial variables
    x = 0;
    xloop = 0;
    dx = 0.1;
    xspan = 20;
	[x,xloop] = FindX(xspan,dx);

	Done = 'Done';
    
%%  Running through the various input diameters 
    for i=1:length(Diameter)
        Diameter1 = Diameter(i);
        [data,t,m,h,n,Ispan] = HHsim(Diameter1,I,PType,t,x,loop,xloop,dx);
        dataholder{1,i} = data;
    end
    
%%  Plots of membrane potential for various input diameters   
    for i=1:length(Diameter)
        figure
        x = dataholder{1,i};
        for n = 2:5:xloop
            plot(t,x(n,:)+(n*15),'b');
            hold on
        end
        xlabel('Time (msec)');
        ylabel('Membrane Potential (mV)');
        if PType == 1
            label=strcat('Action potential propagation for fibre under sine wave stimulation');
        elseif PType == 2
            label=strcat('Action potential propagation for fibre under square wave stimulation');
        elseif PType == 3
            label=strcat('Action potential propagation for fibre under ramp function stimulation');
        end
        title(label);   
        
        if PType == 3
            figure
            yyaxis left
            x = dataholder{1,i};
            plot(t,x(2,:));
            xlabel('Time (msec)');
            ylabel('Membrane Potential (mV)');
            
            yyaxis right
            plot(t,Ispan(2,:));
            ylabel('Stimulus (µAmps)');
            maxIspan = max(Ispan(:));
            ylim([0 (maxIspan+0.01)]);
        end
    end
end

%%	Creating simulation timeframe
function [t,loop] = FindT(tspan,dt)
    loop  = ceil(tspan/dt);
    t = (1:loop)*dt;
end

%%	Creating simulation fibre length
function [x,xloop] = FindX(xspan,dx)
    xloop  = ceil(xspan/dx);
    x = (1:xloop)*dx;
end

%%  Hodgkin-Huxley model
function [V,t,m,h,n,Ispan] = HHsim(defDiameter,defI,defPType,t,x,loop,xloop,dx)    
    Idur = 1;           % duration of stimulation application (in msec)
    defTemp = 6.3;      % environmental temperature (in deg Celsius)
    PulsePos = [1];     % PulsePos is a vector dictating the start time in msec of each successive pulse, eg. [1,17,26,40]
    dt = 0.001;         % time steps for simulation
    freq = 10;          % 0.05
    
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
    V = zeros(xloop,loop);          % Membrane potential
    v = zeros(xloop,loop);          % Normalised membrane potential
    Im = zeros(xloop,loop+1);       % Transmembrane current
    Iion = zeros(xloop,loop+1);     % Ionic currents
	m = zeros(xloop,loop+1);        % Gating variable m
	h = zeros(xloop,loop+1);        % Gating variable h
	n = zeros(xloop,loop+1);        % Gating variable n
    dm = zeros(xloop,loop+1);       % Change in Gating variable m
	dh = zeros(xloop,loop+1);       % Change in Gating variable h
	dn = zeros(xloop,loop+1);       % Change in Gating variable n
	Ispan = zeros(xloop,loop);      % Stimulating current
    alphaN = zeros(xloop,loop+1);   % Gating variable rate constant
    betaN = zeros(xloop,loop+1);    % Gating variable rate constant
    alphaM = zeros(xloop,loop+1);   % Gating variable rate constant
    betaM = zeros(xloop,loop+1);    % Gating variable rate constant
    alphaH = zeros(xloop,loop+1);   % Gating variable rate constant
    betaH = zeros(xloop,loop+1);    % Gating variable rate constant
   
%%	Ispan is the applied current vector to hold all instances of the external 
%   current applied over the simulation period as specified by PulsePos vector    
    if defPType == 1        % sinusoidal stimulus
        Ispan(2,:) = defI*sin(2*pi*freq*(10^(-3))*t);
    elseif defPType == 2    % square periodic pulsatile stimulus  
        Ispan(2,:) = 0.5*defI*((square(2*pi*freq*(10^(-3))*t))+1);
    elseif defPType == 3	% ramp function
        unitstep = t>=0;
        ramp = t.*unitstep*(10^(-3));
        Ispan(2,:) = defI*ramp;
    end
    
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
        for s=1:xloop-1
            v(s,i) = V(s,i) - Vrest;
            alphaN(s,i) = (0.1-0.01*v(s,i))/(exp(1-0.1*v(s,i))-1);
            betaN(s,i) = 0.125*exp(-v(s,i)/80);
            alphaM(s,i) = (2.5-0.1*v(s,i))/(exp(2.5-0.1*v(s,i))-1);
            betaM(s,i) = 4*exp(-v(s,i)/18);
            alphaH(s,i) = 0.07*exp(-v(s,i)/20);
            betaH(s,i) = 1/(exp(3-0.1*v(s,i))+1);

            Iion(s,i) = gNa*(m(s,i)^3)*h(s,i)*(v(s,i)-vNa) + gK*(n(s,i)^4)*(v(s,i)-vK) + gL*(v(s,i)-vL);
            if s-1<=0
                Im(s,i) = (defDiameter/(4*ri*dx*L))*(2*v(s+1,i)-2*v(s,i)) + (Ispan(s,i)/(pi*defDiameter*L));
                dn(s,i) = phi*dt*(alphaN(s,i)*(1-n(s,i)) - betaN(s,i)*n(s,i));
                dm(s,i) = phi*dt*(alphaM(s,i)*(1-m(s,i)) - betaM(s,i)*m(s,i));
                dh(s,i) = phi*dt*(alphaH(s,i)*(1-h(s,i)) - betaH(s,i)*h(s,i));
                V(s,i+1) = V(s,i) + dt*(Im(s,i) - Iion(s,i))/Cm;  
            else
                Im(s,i) = (defDiameter/(4*ri*dx*L))*(v(s-1,i)-2*v(s,i)+v(s+1,i)) + (Ispan(s,i)/(pi*defDiameter*L));
                dn(s,i) = phi*dt*(alphaN(s,i)*(1-n(s,i)) - betaN(s,i)*n(s,i));
                dm(s,i) = phi*dt*(alphaM(s,i)*(1-m(s,i)) - betaM(s,i)*m(s,i));
                dh(s,i) = phi*dt*(alphaH(s,i)*(1-h(s,i)) - betaH(s,i)*h(s,i));
                V(s,i+1) = V(s,i) + dt*(Im(s,i) - Iion(s,i))/Cm;
            end
            n(s,i+1) = n(s,i) + dn(s,i);
            m(s,i+1) = m(s,i) + dm(s,i);
            h(s,i+1) = h(s,i) + dh(s,i);
        end
    end
end


