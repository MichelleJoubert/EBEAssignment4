function [Done] = MultiVar4_1(Diameter,I,Idur)
%   This function simulates the propagation of an action potential along   
%   nerve fibre of two different diameters.
%
%   MultiVar4_1(Diameter,I,Idur) function simulates the Hodgkin-Huxley model 
%   for the squid giant axon constants and user specified values of the 
%   fibre diameter in cm. As output it plots voltage (membrane potential)   
%   time series for each of the various fibre diameters.
%   Diameter is the fibre diameter (in cm). This variable is a list of
%   values to allow for investigation of multiple diameters at once.
%   I is the applied current or external stimulation (in µA)
%   Idur is the duration of stimulation application (in msec) or the
%   stimulus duration type
%
%   Example:
%   MultiVar4_1([0.01],0.5,1) for single pulse
%   MultiVar4_1([0.02,0.01,0.005,0.002],0.5,1) for single pulse for multiple
%   fibre diameters (for Question 4.1)
%   MultiVar4_1([0.01],0.5,2) for constant pulse (for Question 4.6)
%   MultiVar4_1([0.01],0.5,3) for low frequency pulse (for Question 4.6)
%
%%	Simulation timing variables
    t = 0;
    loop = 0;
    dt = 0.001;     % time steps for simulation
    tspan = 50;     % total simulation time (change tspan = 100 for Question 4.6 simulation and tspan = 50 for Question 4.1) 
    [t,loop] = FindT(tspan,dt);
    
%%	Simulation spatial variables
    x = 0;
    xloop = 0;
    dx = 0.05;      % spatial steps
    xspan = 5;      % total fibre length
	[x,xloop] = FindX(xspan,dx);

	Done = 'Done';
    
%%  Running through the various input diameters 
    for i=1:length(Diameter)
        Diameter1 = Diameter(i);
        [data,t,p,Ispan] = HHsim(Diameter1,I,Idur,t,x,loop,dt,xloop,dx);
        dataholder{1,i} = data;
    end
    
%%  Plots of membrane potential for various input diameters   
    for i=1:length(Diameter)
        figure
        k = dataholder{1,i};
        for n = 2:2:xloop-1
            plot(t,k(n,:)+(n*15),'b');
            hold on
        end
        xlabel('Time (msec)');
        label=strcat('Action potential propagation for fibre of diameter ',{' '},num2str(Diameter(i)),'cm');
        title(label);
        
        if Idur == 2 || 3
            figure
            yyaxis left
            x = dataholder{1,i};
            plot(t,x(p,:));
            xlabel('Time (msec)');
            ylabel('Membrane Potential (mV)');

            yyaxis right
            plot(t,Ispan(p,:));
            ylabel('Stimulus (µA)');
            maxIspan = max(Ispan(:));
            minIspan = min(Ispan(:));
            ylim([(minIspan-0.01) (maxIspan+0.01)]);
            label=strcat('Applied current and membrane response for maximum firing rate of fibre');
            title(label); 
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
function [v,t,p,Ispan] = HHsim(defDiameter,defI,Idur,t,x,loop,dt,xloop,dx)    
    defTemp = 6.3;              % environmental temperature (in deg Celsius)
    
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
    if Idur == 1
        Ispan(p,3000:4000) = defI;
    elseif Idur == 2
        Ispan(p,:) = defI;
	elseif Idur == 3
        fs = 1000;  
        freq = 10;  
        Ispan(p,:) = defI*sin(2*pi*freq*(1/fs)*t);  
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


