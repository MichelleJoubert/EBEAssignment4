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
%   MultiVar4_7([0.01],0.5,1) for sine wave (for Question 4.7 i)
%   MultiVar4_7([0.01],0.5,2) for square wave (for Question 4.7 ii)
%   MultiVar4_7([0.01],0.3,3) for ramp function (for Question 4.8)
%
%%	Simulation timing variables
    t = 0;
    loop = 0;
    dt = 0.001;     % time steps for simulation
    tspan = 100;    % total simulation time
    [t,loop] = FindT(tspan,dt);
    
%%	Simulation spatial variables
    x = 0;
    xloop = 0;
    dx = 0.1;       % spatial steps
    xspan = 10;     % total fibre length
	[x,xloop] = FindX(xspan,dx);

	Done = 'Done';
    
%%  Running through the various input diameters 
    for i=1:length(Diameter)
        Diameter1 = Diameter(i);
        [data,t,m,h,n,Ispan,p] = HHsim(Diameter1,I,PType,t,x,loop,dt,xloop,dx);
        dataholder{1,i} = data;
    end
    
%%  Plots of membrane potential for various input diameters   
    for i=1:length(Diameter)
        figure
        subplot(2,1,1)
        x = dataholder{1,i};
        for k = 2:1:xloop-1
            plot(t,x(k,:)+(k*15),'b');
            hold on
        end
        xlabel('Time (msec)');
        if PType == 1
            label=strcat('Action potential propagation for fibre under sine wave stimulation');
        elseif PType == 2
            label=strcat('Action potential propagation for fibre under square wave stimulation');
        elseif PType == 3
            label=strcat('Action potential propagation for fibre under ramp function stimulation');
        end
        title(label);   
        
        subplot(2,1,2)
        yyaxis left
        x = dataholder{1,i};
        plot(t,x(p,:),'LineWidth',2);
        xlabel('Time (msec)');
        ylabel('Membrane Potential (mV)');
            
        yyaxis right
        plot(t,Ispan(p,:),'LineWidth',2);
        ylabel('Stimulus (µA)');
        maxIspan = max(Ispan(:));
        minIspan = min(Ispan(:));
        ylim([(minIspan-0.01) (maxIspan+0.01)]);
        if PType == 1
            label=strcat('Applied current for fibre under sine wave stimulation');
        elseif PType == 2
            label=strcat('Applied current for fibre under square wave stimulation');
        elseif PType == 3
            label=strcat('Applied current for fibre under ramp function stimulation');
        end
        title(label); 
        
        if PType == 3
            figure
            plot(t,m(p,:),t,n(p,:),t,h(p,:),'LineWidth',2);
            xlabel('Time (msec)');
            ylabel('Magnitude of gating variables');
            title('Gating variables at stimulus site for ramp function stimulation');  
            legend('m','n', 'h'); 
        end
    end
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
function [v,t,m,h,n,Ispan,p] = HHsim(defDiameter,defI,defPType,t,x,loop,dt,xloop,dx)    
    defTemp = 6.3;      % environmental temperature (in deg Celsius)
    freq = 120;         % Frequency for sine and square waves
    
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
    alphaN = zeros(xloop,loop+1);   % Gating variable rate constant
    betaN = zeros(xloop,loop+1);    % Gating variable rate constant
    alphaM = zeros(xloop,loop+1);   % Gating variable rate constant
    betaM = zeros(xloop,loop+1);    % Gating variable rate constant
    alphaH = zeros(xloop,loop+1);   % Gating variable rate constant
    betaH = zeros(xloop,loop+1);    % Gating variable rate constant
   
%%	Ispan is the applied current vector to hold all instances of the external   
    p = xloop/2;
    fs = 1000;              % adjusting sampling frequency by 1000 samples per second (taking into account dt, that gives fs of 1,000,000 per second)
    if defPType == 1        % sinusoidal stimulus
        Ispan(p,:) = defI*sin(2*pi*freq*(1/fs)*t);
    elseif defPType == 2    % square periodic pulsatile stimulus  
        Ispan(p,:) = defI*((square(2*pi*freq*(1/fs)*t))+1);
    elseif defPType == 3	% ramp function
        unitstep = t>=0;
        ramp = t.*unitstep;
        Ispan(p,:) = defI*0.01*ramp;
%         Ispan(p,3000:4000) = defI;    % use to demonstate the 0.3µA is
%         sufficient to elicit AP under normal conditions
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
        for s=2:xloop-1
            v(s,i) = V(s,i) - Vrest;
            alphaN(s,i) = (0.1-0.01*v(s,i))/(exp(1-0.1*v(s,i))-1);
            betaN(s,i) = 0.125*exp(-v(s,i)/80);
            alphaM(s,i) = (2.5-0.1*v(s,i))/(exp(2.5-0.1*v(s,i))-1);
            betaM(s,i) = 4*exp(-v(s,i)/18);
            alphaH(s,i) = 0.07*exp(-v(s,i)/20);
            betaH(s,i) = 1/(exp(3-0.1*v(s,i))+1);

            Iion(s,i) = gNa*(m(s,i)^3)*h(s,i)*(v(s,i)-vNa) + gK*(n(s,i)^4)*(v(s,i)-vK) + gL*(v(s,i)-vL);
            Im(s,i) = (defDiameter/(4*ri*dx*L))*(v(s-1,i)-2*v(s,i)+v(s+1,i)) + (Ispan(s,i)/(pi*defDiameter*L));
            dn(s,i) = phi*dt*(alphaN(s,i)*(1-n(s,i)) - betaN(s,i)*n(s,i));
            dm(s,i) = phi*dt*(alphaM(s,i)*(1-m(s,i)) - betaM(s,i)*m(s,i));
            dh(s,i) = phi*dt*(alphaH(s,i)*(1-h(s,i)) - betaH(s,i)*h(s,i));
                
            V(s,i+1) = V(s,i) + dt*(Im(s,i) - Iion(s,i))/Cm;
            v(s,i+1) = v(s,i) + dt*(Im(s,i) - Iion(s,i))/Cm;
            n(s,i+1) = n(s,i) + dn(s,i);
            m(s,i+1) = m(s,i) + dm(s,i);
            h(s,i+1) = h(s,i) + dh(s,i);
        end
    end
end


