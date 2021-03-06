function [Done] = MultiVar4_3(Distance,I)
%   This function simulates the propagation of an action potential along   
%   nerve fibre for various electrode distances away from the nerve fibre,    
%	stimulated by a single anodic or cathodic monophasic pulse.
%   
%   MultiVar4_3(Distance,I) function simulates the Hodgkin-Huxley model 
%   for the squid giant axon constants and user specified values of 
%   electrode distance away from the fibre for specified input stimulus.      
%   As output it plots voltage(membrane potential)time series   
%   for the specified electrode distances from the fibre and current. 
%   Distance is the electrode distance away from the axon (in cm). This 
%   variable is a list of values to allow for investigation of multiple 
%   distances at once.
%   I is the applied current or external stimulation (in �A)
%   The default simulation values for anodic stimulation are: z = 0.08cm, I = 0.5�A, Idur = 1msec
%   The default simulation values for cathodic stimulation are: z = 0.075cm, I = -0.5�A, Idur = 1msec
%   
%   Example:
%   ANODIC
%   MultiVar4_3([0.06,0.08,0.1],0.5) for anodic stimulation with differing electrode distances z (for Question 4.3)
%   MultiVar4_3([0.08],0.25) and MultiVar4_3([0.08],0.75) for anodic stimulation with differing current amplitudes I (for Question 4.3)
%   MultiVar4_3([0.08],0.5) for anodic stimulation with differing current application durations Idur with changes to duration made in code (for Question 4.3)
%   
%   CATHODIC
%   MultiVar4_3([0.05,0.075,0.1],-0.5) for cathodic stimulation with differing electrode distances z (for Question 4.4)
%   MultiVar4_3([0.075],-0.25) and MultiVar4_3([0.075],-0.75) for cathodic stimulation with differing current amplitudes I (for Question 4.4)
%   MultiVar4_3([0.075],-0.5) for cathodic stimulation with differing current application durations Idur with changes to duration made in code (for Question 4.4)
%
%%	Simulation timing variables
    t = 0;
    loop = 0;
    dt = 0.001;     % time steps for simulation
    tspan = 25;     % total simulation time         
    [t,loop] = FindT(tspan,dt);
    
%%	Simulation spatial variables
    x = 0;
    xloop = 0;
    dx = 0.02;      % spatial steps
    xspan = 2;      % total fibre length
	[x,xloop] = FindX(xspan,dx);

	Done = 'Done';
    
%%  Running through the various input diameters 
    for i=1:length(Distance)
        Distance1 = Distance(i);
        [data,t,f,x] = HHsim(Distance1,I,t,x,loop,dt,xloop,dx);
        dataholder{1,i} = data;
        fholder{1,i} = f;
    end
    
%%  Plots of membrane potential for various input diameters   
    for i=1:length(Distance)
        figure
        k = dataholder{1,i};
        for n = 2:2:xloop-1
            plot(t,k(n,:)+(n*15),'b');
            hold on
        end
        xlabel('Time (msec)');
        label=strcat('AP propagation with electrode a distance of ',{' '},num2str(Distance(i)),'cm from axon');
        title(label);
      
        figure
        k = fholder{1,i};
        plot(x,k(:,3500));
        xlabel('Fibre length (cm)');
        ylabel('Activating function (mV)');
        label=strcat('Activation function with electrode a distance of ',{' '},num2str(Distance(i)),'cm from axon');
        title(label);
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
function [v,t,f,x] = HHsim(z,defI,t,x,loop,dt,xloop,dx)    
    defTemp = 6.3;      % environmental temperature (in deg Celsius)
    Diameter = 0.0012;  % Diameter of fibre (in cm)
    
%%	Constants and intial values for squid giant axon  
    gNa = 120;                  % Conductance of sodium channels (in m.mho/cm^2)
	vNa= 115;                   % Ionic potential for sodium channels (in mV)
	gK = 36;                    % Conductance of potassium channels (in m.mho/cm^2)
	vK= -12;                    % Ionic potential for potassium channels (in mV)
	gL= 0.3;                    % Conductance of leakage channels (in m.mho/cm^2)
	vL= 10.6;                   % Ionic potential for leakage channels (in mV)
	Cm = 1;                     % Capacitance of membrane (in �F/cm^2)
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
    xe = zeros(xloop,loop+1);       % Distance along axon away from electrode
    f = zeros(xloop+1,loop+1);      % Activating function
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
    q = xloop/6;
    p = xloop/2;
    Ispan(p-q:p+q,3000:4000) = defI;    % current application duration of 1msec (DEFAULT)
%     Ispan(p-q:p+q,3000:8000) = defI;    % current application duration of 5msec
%     Ispan(p-q:p+q,3000:13000) = defI;   % current application duration of 10msec
    
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
            electrodePos = p*dx;            % position of electrode
            xe(s,i) = s*dx-electrodePos;	% gives fibre distance away from electrode for different dx segments
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
            if s-1<=0
                Im(s,i) = (Diameter/(4*ri*dx*L))*(2*v(s+1,i)-2*v(s,i));
                f(s,i) = (Diameter/(4*ri*dx*L))*((re*Ispan(s,i))/(4*pi))*(((xe(s,i)^2)+(z^2))^(-5/2))*(2*(xe(s,i)^2)-(z^2));
                dn(s,i) = phi*dt*(alphaN(s,i)*(1-n(s,i)) - betaN(s,i)*n(s,i));
                dm(s,i) = phi*dt*(alphaM(s,i)*(1-m(s,i)) - betaM(s,i)*m(s,i));
                dh(s,i) = phi*dt*(alphaH(s,i)*(1-h(s,i)) - betaH(s,i)*h(s,i));
                V(s,i+1) = V(s,i) + dt*(Im(s,i) + f(s,i) - Iion(s,i))/Cm;
                v(s,i+1) = v(s,i) + dt*(Im(s,i) + f(s,i) - Iion(s,i))/Cm;
            else
                Im(s,i) = (Diameter/(4*ri*dx*L))*(v(s-1,i)-2*v(s,i)+v(s+1,i));
                f(s,i) = (Diameter/(4*ri*dx*L))*((re*Ispan(s,i))/(4*pi))*(((xe(s,i)^2)+(z^2))^(-5/2))*(2*(xe(s,i)^2)-(z^2));
                dn(s,i) = phi*dt*(alphaN(s,i)*(1-n(s,i)) - betaN(s,i)*n(s,i));
                dm(s,i) = phi*dt*(alphaM(s,i)*(1-m(s,i)) - betaM(s,i)*m(s,i));
                dh(s,i) = phi*dt*(alphaH(s,i)*(1-h(s,i)) - betaH(s,i)*h(s,i));
                V(s,i+1) = V(s,i) + dt*(Im(s,i) + f(s,i) - Iion(s,i))/Cm;
                v(s,i+1) = v(s,i) + dt*(Im(s,i) + f(s,i) - Iion(s,i))/Cm;
            end
            n(s,i+1) = n(s,i) + dn(s,i);
            m(s,i+1) = m(s,i) + dm(s,i);
            h(s,i+1) = h(s,i) + dh(s,i);
        end
    end
end


