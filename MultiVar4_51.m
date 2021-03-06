function [Done] = MultiVar4_5(Distance,I)
%   This function simulates the propagation of an action potential along   
%   nerve fibre of various distances away from the nerve fibre stimulated by    
%	a bipolar stimulation pulse.
%
%   MultiVar4_5(Distance,I) function simulates the Hodgkin-Huxley model 
%   for the squid giant axon constants and user specified values of 
%   stimulation distance from the fibre in cm for specified input stimulus. 
%   As output it plots voltage(membrane potential)time series for each of  
%   the various fibre diameters. 
%   Distance is the electrode distance away from the axon (in mm). This 
%   variable is a list of values to allow for investigation of multiple 
%   distances at once.
%   I is the applied current or external stimulation (in �A)
%   
%   Example:
%   MultiVar4_5([0.025],10000) for bipolar stimulation with anode then
%   cathode placed along fibre 
%   MultiVar4_5([0.025],2000)for bipolar stimulation with cathode then
%   anode placed along fibre for cathode propagation (swop signs of defI in Ispan and electrodePos 1 and 2)
%
%%	Simulation timing variables
    t = 0;
    loop = 0;
    dt = 0.01;
    tspan = 200; %400 for full prop to end of fibre
    [t,loop] = FindT(tspan,dt);
    
%%	Simulation spatial variables
    x = 0;
    xloop = 0;
    dx = 0.05;
    xspan = 10;
	[x,xloop] = FindX(xspan,dx);

	Done = 'Done';
    
%%  Running through the various input diameters 
    for i=1:length(Distance)
        Distance1 = Distance(i);
        [data,m,h,n,Ispan,f] = HHsim(Distance1,I,t,x,loop,xloop,dx);
        dataholder{1,i} = data;
    end
    
%%  Plots of membrane potential for various input diameters   
    for i=1:length(Distance)
        figure
        k = dataholder{1,i};
        for n = 1:1:xspan
            plot(t,k(n,:)+(n*15),'b');
            hold on
        end
        xlabel('Time (msec)');
        ylabel('Membrane Potential (mV)');
        label=strcat('Action potential propagation for fibre with electrode a distance ',{' '},num2str(Distance(i)),'mm away from axon');
        title(label);
        
        figure
        plot(x,f(:,1500));
        xlabel('Fibre length (mm)');
        ylabel('Activating function (mV)');
        xlim([0 5])
        label=strcat('Action potential propagation for fibre with electrode a distance of ',{' '},num2str(Distance(i)),'mm away from axon');
        title(label);
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
function [V,t,m,h,n,Ispan,f] = HHsim(z,defI,t,x,loop,xloop,dx)    
    defIdur = 1;         % duration of stimulation application (in msec)
    defTemp = 6.3;      % environmental temperature (in deg Celsius)
    PulsePos = [1];     % PulsePos is a vector dictating the start time in msec of each successive pulse, eg. [1,17,26,40]
    dt = 0.001;         % time steps for simulation
    Diameter = 0.012;   % Diameter of fibre (in cm)
    
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
    re = 0.3;                   % specific extracellular resistance (in kOhm.cm)
    m0 = 0.05;                  % Initial value of gating variable m
    n0 = 0.32;                  % Initial value of gating variable n
    h0 = 0.6;                   % Initial value of gating variable h
    Vrest = -70;                % Membrane potential at rest (in mV)
    V0 = Vrest;                 % Initial value of actual membrane potential (in mV)
    v0 = V0 - Vrest;            % Initial value of normalised membrane potential (in mV)    
    
%%	Initializing variable vectors  
    V = zeros(xloop+1,loop);        % Membrane potential
    v = zeros(xloop+1,loop);        % Normalised membrane potential
    Im = zeros(xloop,loop+1);       % Transmembrane current
    Iion = zeros(xloop,loop+1);     % Ionic currents
    xe1 = zeros(xloop,loop+1);      % Distance along axon away from electrode 1
    xe2 = zeros(xloop,loop+1);      % Distance along axon away from electrode 2
    f = zeros(xloop,loop);          % Activating function
	m = zeros(xloop,loop+1);        % Gating variable m
	h = zeros(xloop,loop+1);        % Gating variable h
	n = zeros(xloop,loop+1);        % Gating variable n
    dm = zeros(xloop,loop+1);       % Change in Gating variable m
	dh = zeros(xloop,loop+1);       % Change in Gating variable h
	dn = zeros(xloop,loop+1);       % Change in Gating variable h
	Ispan = zeros(xloop,loop+1);    % Stimulating current
    alphaN = zeros(xloop,loop+1);   % Gating variable rate constant
    betaN = zeros(xloop,loop+1);    % Gating variable rate constant
    alphaM = zeros(xloop,loop+1);   % Gating variable rate constant
    betaM = zeros(xloop,loop+1);    % Gating variable rate constant
    alphaH = zeros(xloop,loop+1);   % Gating variable rate constant
    betaH = zeros(xloop,loop+1);    % Gating variable rate constant

%%	Ispan is the applied current vector to hold all instances of the external 
%   current applied over the simulation period as specified by PulsePos vector    
	IspanEnd = ceil(defIdur/dt); 
    for i=1:length(PulsePos)
        IspanPulsePos = (ceil(PulsePos(i)/dt));
        if IspanPulsePos + IspanEnd < length(Ispan)
            PlusSpan = IspanPulsePos + IspanEnd;
        else
            PlusSpan = length(Ispan);
        end
        Ispan(14,IspanPulsePos:PlusSpan) = Ispan(14,IspanPulsePos:PlusSpan)+defI;
        Ispan(16,IspanPulsePos:PlusSpan) = Ispan(16,IspanPulsePos:PlusSpan)-defI;
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
    for s=1:1:xloop-1
        for i=1:1:loop-1
            electrodePos1 = ceil(loop/10000)+0.01;
            xe1(s,i) = dx*(s-electrodePos1);
            electrodePos2 = ceil(loop/10000)-0.01;
            xe2(s,i) = dx*(s-electrodePos2);
        end
    end
 
%   Solution for fibre section after second electrode
    for i=1:1:loop-1
        for s=15:1:xloop-1
            v(s,i) = V(s,i) - Vrest;
            alphaN(s,i) = (0.1-0.01*v(s,i))/(exp(1-0.1*v(s,i))-1);
            betaN(s,i) = 0.125*exp(-v(s,i)/80);
            alphaM(s,i) = (2.5-0.1*v(s,i))/(exp(2.5-0.1*v(s,i))-1);
            betaM(s,i) = 4*exp(-v(s,i)/18);
            alphaH(s,i) = 0.07*exp(-v(s,i)/20);
            betaH(s,i) = 1/(exp(3-0.1*v(s,i))+1);

            xe1(s,i) = dx*(s-electrodePos1);
            xe2(s,i) = dx*(s-electrodePos2);
            Iion(s,i) = gNa*(m(s,i)^3)*h(s,i)*(v(s,i)-vNa) + gK*(n(s,i)^4)*(v(s,i)-vK) + gL*(v(s,i)-vL);
            if s-1<=0
                Im(s,i) = (Diameter/(4*ri*dx*L))*(2*v(s+1,i)-2*v(s,i));
                f(s,i) = ((re*Ispan(s,i))/(4*pi))*((((xe1(s,i)^2)+(z^2))^(-5/2))*(2*(xe1(s,i)^2)-(z^2))+(((xe2(s,i)^2)+(z^2))^(-5/2))*(2*(xe2(s,i)^2)-(z^2)));
                dn(s,i) = phi*dt*(alphaN(s,i)*(1-n(s,i)) - betaN(s,i)*n(s,i));
                dm(s,i) = phi*dt*(alphaM(s,i)*(1-m(s,i)) - betaM(s,i)*m(s,i));
                dh(s,i) = phi*dt*(alphaH(s,i)*(1-h(s,i)) - betaH(s,i)*h(s,i));
                V(s,i+1) = V(s,i) + dt*(Im(s,i) + f(s,i) - Iion(s,i))/Cm;  
            else
                Im(s,i) = (Diameter/(4*ri*dx*L))*(v(s-1,i)-2*v(s,i)+v(s+1,i));
                f(s,i) = ((re*Ispan(s,i))/(4*pi))*((((xe1(s,i)^2)+(z^2))^(-5/2))*(2*(xe1(s,i)^2)-(z^2))+(((xe2(s,i)^2)+(z^2))^(-5/2))*(2*(xe2(s,i)^2)-(z^2)));
                dn(s,i) = phi*dt*(alphaN(s,i)*(1-n(s,i)) - betaN(s,i)*n(s,i));
                dm(s,i) = phi*dt*(alphaM(s,i)*(1-m(s,i)) - betaM(s,i)*m(s,i));
                dh(s,i) = phi*dt*(alphaH(s,i)*(1-h(s,i)) - betaH(s,i)*h(s,i));
                V(s,i+1) = V(s,i) + dt*(Im(s,i) + f(s,i) - Iion(s,i))/Cm;  
            end
            n(s,i+1) = n(s,i) + dn(s,i);
            m(s,i+1) = m(s,i) + dm(s,i);
            h(s,i+1) = h(s,i) + dh(s,i);
        end
    
%   Solution for fibre section before first electrode        
        for s=14:-1:1
            v(s,i) = V(s,i) - Vrest;
            alphaN(s,i) = (0.1-0.01*v(s,i))/(exp(1-0.1*v(s,i))-1);
            betaN(s,i) = 0.125*exp(-v(s,i)/80);
            alphaM(s,i) = (2.5-0.1*v(s,i))/(exp(2.5-0.1*v(s,i))-1);
            betaM(s,i) = 4*exp(-v(s,i)/18);
            alphaH(s,i) = 0.07*exp(-v(s,i)/20);
            betaH(s,i) = 1/(exp(3-0.1*v(s,i))+1);

            xe1(s,i) = dx*(s-electrodePos1);
            xe2(s,i) = dx*(s-electrodePos2);
            Iion(s,i) = gNa*(m(s,i)^3)*h(s,i)*(v(s,i)-vNa) + gK*(n(s,i)^4)*(v(s,i)-vK) + gL*(v(s,i)-vL);
            if s-1<=0
                Im(s,i) = (Diameter/(4*ri*dx*L))*(2*v(s+1,i)-2*v(s,i));
                f(s,i) = ((re*Ispan(s,i))/(4*pi))*((((xe1(s,i)^2)+(z^2))^(-5/2))*(2*(xe1(s,i)^2)-(z^2))+(((xe2(s,i)^2)+(z^2))^(-5/2))*(2*(xe2(s,i)^2)-(z^2)));
                dn(s,i) = phi*dt*(alphaN(s,i)*(1-n(s,i)) - betaN(s,i)*n(s,i));
                dm(s,i) = phi*dt*(alphaM(s,i)*(1-m(s,i)) - betaM(s,i)*m(s,i));
                dh(s,i) = phi*dt*(alphaH(s,i)*(1-h(s,i)) - betaH(s,i)*h(s,i));
                V(s,i+1) = V(s,i) + dt*(Im(s,i) + f(s,i) - Iion(s,i))/Cm;  
            else
                Im(s,i) = (Diameter/(4*ri*dx*L))*(v(s-1,i)-2*v(s,i)+v(s+1,i));
                f(s,i) = ((re*Ispan(s,i))/(4*pi))*((((xe1(s,i)^2)+(z^2))^(-5/2))*(2*(xe1(s,i)^2)-(z^2))+(((xe2(s,i)^2)+(z^2))^(-5/2))*(2*(xe2(s,i)^2)-(z^2)));
                dn(s,i) = phi*dt*(alphaN(s,i)*(1-n(s,i)) - betaN(s,i)*n(s,i));
                dm(s,i) = phi*dt*(alphaM(s,i)*(1-m(s,i)) - betaM(s,i)*m(s,i));
                dh(s,i) = phi*dt*(alphaH(s,i)*(1-h(s,i)) - betaH(s,i)*h(s,i));
                V(s,i+1) = V(s,i) + dt*(Im(s,i) + f(s,i) - Iion(s,i))/Cm;  
            end
            n(s,i+1) = n(s,i) + dn(s,i);
            m(s,i+1) = m(s,i) + dm(s,i);
            h(s,i+1) = h(s,i) + dh(s,i);
        end
    end
end


