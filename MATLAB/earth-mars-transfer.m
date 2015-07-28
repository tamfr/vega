clc

clear all
tic
% Important Remarks: 
% 1. J2000 Used to Establish All Constants
% 2. For Analysis, Earth's Perihelion is the Epoch (JDN 2455565.547 [i.e. CE 2011 January 04 01:07:40.8 UT  Tuesday])
% 3. Inclination and Orbital Precession are Neglected

% Heliocentric Orbital System Constants (Source: NASA JPL Horizons Web Interface JDN: 2455562.5-2455927.5)

muSun = 2.9591230378107664*10^(-4)*149597870^3/86400^2; % Sun Standard Gravitational Parameter [km^3/s^2]
LonPerMars = 336.04084*pi/180; % Mars Longitude of Perihelion [radians] (http://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html)
LonPerEarth = 102.94719*pi/180; % Earth Longitude of Perihelion [radians] (http://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html)
Theta = LonPerMars-LonPerEarth; % Angle Between Earth's Perihelion and Mars's Perihelion [radians]
JDN0 = 2455565.2736; % Julian Day Number of Last Epoch (i.e. CE 2011 January 03 18:34:00 UT  Tuesday) 

% Constant Orbital Elements for Mars (Source: NASA JPL Horizons Web Interface JDN: 2455562.5-2455927.5)

eMars = 0.093439031; % Mars Eccentricity
ApMars = 249234450.9; % Mars Aphelion [km]
PerMars = ApMars*(1-eMars)/(1+eMars); % Mars Perihelion [km]
aMars = ApMars/(1+eMars); % Mars Semi-Major Axis [km]
TMars = 2*pi*(aMars^3/muSun)^(1/2); % Mars Orbital Period [s]
rMars = 3389.92; % Mean Radius of Mars [km]
muMars = 42828.3; % Mars Standard Graviational Parameter [km^3/s^2]
gMars = 3.71*10^(-3); % Gravitational Acceleration of Mars

% Constant Orbital Elements for Earth (Source: NASA JPL Horizons Web Interface JDN: 2455562.5-2455927.5)

eEarth = 0.016741967; % Earth Eccentricity
ApEarth = 152105805.7; % Earth Aphelion [km]
PerEarth = ApEarth*(1-eEarth)/(1+eEarth); % Earth Perihelion [km]
aEarth = ApEarth/(1+eEarth); % Earth Semi-Major Axis [km]
TEarth = 2*pi*(aEarth^3/muSun)^(1/2); % Earth Orbital Period [s]
rEarth = 6378.136; % Radius of Earth at Equator [km]
muEarth = 398600.440; % Earth Standard Graviational Parameter [km^3/s^2]
gEarth = 9.81*10^(-3); % Gravitational Acceleration of Earth

% Initilization of Variable Parameters

RM = 0; % Orbital Radius of Mars
RE = 0; % Orbital Radius of Earth
RT = 0; % Orbital Radius of Transfer Orbit
fM = 0; % True Anomaly of Mars
fE = 0; % True Anomaly of Earth
mfM = 0; % True Anomaly of Mars Relative to Earth's Perihelion
fT = 0; % True Anomaly of Transfer Orbit
aT = 0; % Transfer Orbit Semi-Major Axis
eT = 0; % Transfer Orbit Eccentricity 
TT = 0; % 1/2 Transfer Orbit Period
tT = 0; % Transfer Time
MM = 0; % Mars Mean Anomaly
MT = 0; % Transfer Orbit Mean Anomaly
EM = 0; % Mars Eccentric Anomaly
ET = 0; % Transfer Orbit Eccentric Anomaly
VTarrv = 0; % Transfer Orbit Arrival Velocity
VMarsArrv = 0; % Mars Orbit Velocity Upon Arrival
phi = 0; % Angle With Respect to Orbit Radius of Velocity Vector of Transfer Orbit

% Approximate Minimum Delta V Required for Transit (Minimum Hohmann)

RE = aEarth*(1-eEarth^2)/(1+eEarth*cos(Theta-pi));
aT = (PerMars+RE)/2;
eT = (PerMars-RE)/(PerMars+RE);
DVmin0 = (((1+eT)/(1-eT))*muSun/aT)^(1/2)-(muSun*(2/RE-1/aEarth))^(1/2);
TT = pi*(aT^3/muSun)^(1/2);
MM = (2*pi/TMars)*(TMars-TT);

    % Iteration Loop to Find Eccentric Anomaly
    
    E0 = pi/3; % Initilize Eccentric Anomaly for Loop
    E1 = -pi/2; % Initilize Eccentric Anomaly for Loop

    while abs(E1/E0-1) > 0.0005 % Condition On Which to Run Loop
        if E1 > 2*pi
            E1 = E1-2*pi;
        elseif E1 < 0
            E1 = E1+pi;
        else
            E0 = E1;
            E1 = E0-(E0-eMars*sin(E0)-MM)/(1-eMars*cos(E0)); % Newton's Root Finding Method
        end
    end
    
    fM1 = 2*atan(((1+eMars)/(1-eMars))^(1/2)*tan(E1/2)); % Mars True Anomaly for Min Hohmann Transfer
    
    % Test to Make Sure True Anomaly is Between 0 and 360 Degrees
    
    if fM1 < 0 
        fM1 = fM1+2*pi;
    elseif fM1 > 2*pi;
        fM1 = fM1-2*pi;
    end
    
    % The Following Uses the Method of Patched Conics to Determine the Real Delta V
    
    aPark = 330; % Altitude of Parking Orbit Above Earth [km]
    PerHyper = aPark + rEarth; % Periapsis Distance From Center of Earth
    ah = 0; % Initializes Hyperbolic Semi-Major Axis
    eh = 0; % Initializes Hyperbolic Eccentricity
    DVactual0 = 0; % Initializes Actual Delta V
    
    ah = muEarth/DVmin0^2;
    eh = PerHyper/ah+1; % Determines Hyperbolic Eccentricity
    DVactual0 = ((eh+1)/(eh-1))^(1/2)*DVmin0-(muEarth/PerHyper)^(1/2);
    
    % Delta V to Low Mars Orbit With Patched Conics
    
    RM = aMars*(1-eMars^2)/(1+eMars*cos(0));
    DVminArrv0 = (muSun*(2/RM-1/aMars))^(1/2)-(((1-eT)/(1+eT))*muSun/aT)^(1/2);
    aPark = 330; % Altitude of Parking Orbit Above Mars [km]
    PerHyper = aPark + rMars; % Periapsis Distance From Center of Mars
    ah = 0; % Initializes Hyperbolic Semi-Major Axis
    eh = 0; % Initializes Hyperbolic Eccentricity
    DVactualArrv0 = 0; % Initializes Actual Delta V
    ah = muMars/DVminArrv0^2;
    eh = PerHyper/ah+1; % Determines Hyperbolic Eccentricity
    DVactualArrv0 = ((eh+1)/(eh-1))^(1/2)*DVminArrv0-(muMars/PerHyper)^(1/2);
    
    % Total DeltaV for Transit to Mars
    
    DVtotal0 = DVactual0+DVactualArrv0;
    
% Parabolic Delta V at Point of Min Hohmann Transfer

DVpara0 = (2*muSun/RE)^(1/2)-(muSun*(2/RE-1/aEarth))^(1/2);

    % The Following Uses the Method of Patched Conics to Determine the Real Delta V
    
    aPark = 330; % Altitude of Parking Orbit Above Earth [km]
    PerHyper = aPark + rEarth; % Periapsis Distance From Center of Earth
    ah = 0; % Initializes Hyperbolic Semi-Major Axis
    eh = 0; % Initializes Hyperbolic Eccentricity
    DVactual2 = 0; % Initializes Actual Delta V
    
    ah = muEarth/DVpara0^2;
    eh = PerHyper/ah+1; % Determines Hyperbolic Eccentricity
    DVactual2 = ((eh+1)/(eh-1))^(1/2)*DVpara0-(muEarth/PerHyper)^(1/2);

% Approximate Minimum Delta V Required for Transit (Maximum Hohmann)

RE = aEarth*(1-eEarth^2)/(1+eEarth*cos(Theta));
aT = (ApMars+RE)/2;
eT = (ApMars-RE)/(ApMars+RE);
DVmin1 = (((1+eT)/(1-eT))*muSun/aT)^(1/2)-(muSun*(2/RE-1/aEarth))^(1/2);
TT = pi*(aT^3/muSun)^(1/2);
MM = (2*pi/TMars)*(TMars/2-TT); % Mars Mean Anomaly for Transfer

    % Iteration Loop to Find Eccentric Anomaly
    
    E2 = pi/3; % Initilize Eccentric Anomaly for Loop
    E3 = -pi/2; % Initilize Eccentric Anomaly for Loop

    while abs(E3/E2-1) > 0.0005 % Condition On Which to Run Loop
        if E3 > 2*pi
            E3 = E3-2*pi;
        elseif E3 < 0
            E3 = E3+pi;
        else    
            E2 = E3;
            E3 = E2-(E2-eMars*sin(E2)-MM)/(1-eMars*cos(E2)); % Newton's Root Finding Method
        end
    end
    
    fM2 = 2*atan(((1+eMars)/(1-eMars))^(1/2)*tan(E3/2)); % Mars True Anomaly for Maximum Hohmann Transfer
    
    % Test to Make Sure True Anomaly is Between 0 and 360 Degrees
    
    if fM2 < 0
        fM2 = fM2+2*pi;
    elseif fM2 > 2*pi;
        fM2 = fM2-2*pi;
    end
    
     % The Following Uses the Method of Patched Conics to Determine the Real Delta V
    
    aPark = 330; % Altitude of Parking Orbit Above Earth [km]
    PerHyper = aPark + rEarth; % Periapsis Distance From Center of Earth
    ah = 0; % Initializes Hyperbolic Semi-Major Axis
    eh = 0; % Initializes Hyperbolic Eccentricity
    DVactual1 = 0; % Initializes Actual Delta V
    
    ah = muEarth/DVmin1^2;
    eh = PerHyper/ah+1; % Determines Hyperbolic Eccentricity
    DVactual1 = ((eh+1)/(eh-1))^(1/2)*DVmin1-(muEarth/PerHyper)^(1/2);
    
    % Delta V to Low Mars Orbit With Patched Conics
    
    RM = aMars*(1-eMars^2)/(1+eMars*cos(pi));
    DVminArrv1 = (muSun*(2/RM-1/aMars))^(1/2)-(((1-eT)/(1+eT))*muSun/aT)^(1/2);
    aPark = 330; % Altitude of Parking Orbit Above Mars [km]
    PerHyper = aPark + rMars; % Periapsis Distance From Center of Mars
    ah = 0; % Initializes Hyperbolic Semi-Major Axis
    eh = 0; % Initializes Hyperbolic Eccentricity
    DVactualArrv1 = 0; % Initializes Actual Delta V
    ah = muMars/DVminArrv1^2;
    eh = PerHyper/ah+1; % Determines Hyperbolic Eccentricity
    DVactualArrv1 = ((eh+1)/(eh-1))^(1/2)*DVminArrv1-(muMars/PerHyper)^(1/2);
    
    % Total DeltaV for Transit to Mars
    
    DVtotal1 = DVactual1+DVactualArrv1;
    
% Parabolic Delta V at Point of Max Hohmann Transfer
    
DVpara1 = (2*muSun/RE)^(1/2)-(muSun*(2/RE-1/aEarth))^(1/2);

    % The Following Uses the Method of Patched Conics to Determine the Real Delta V
    
    aPark = 330; % Altitude of Parking Orbit Above Earth [km]
    PerHyper = aPark + rEarth; % Periapsis Distance From Center of Earth
    ah = 0; % Initializes Hyperbolic Semi-Major Axis
    eh = 0; % Initializes Hyperbolic Eccentricity
    DVactual3 = 0; % Initializes Actual Delta V
    
    ah = muEarth/DVpara1^2;
    eh = PerHyper/ah+1; % Determines Hyperbolic Eccentricity
    DVactual3 = ((eh+1)/(eh-1))^(1/2)*DVpara1-(muEarth/PerHyper)^(1/2);
    
% Optimize Conjunction Orbit Transfer for Set Departure Delta V Based On Transfer Time

fE = Theta-pi; % Reinitializes True Anomaly of Earth
DVmax(1) = 5; % Set Maximum Delta V for Earth Departure (Not Considering Planetary Escape)
N = 20; % Amount of Delta V Solutions
M=40;
StepSizeN = (DVmax(1)-DVmin1-.0007)/N;
StepSize = (2*pi)/M; % Initializes Step Size
tT = zeros(N,M);
tTRP = zeros(N,M);
DVactualArrv2 = zeros(N,M);
DVactualArrv2RP = zeros(N,M);
fM = zeros(N,M);
fMRP = zeros(N,M);
for j = 1:N
    fE = zeros(1,M);
    for i = 1:M
       
        RE = aEarth*(1-eEarth^2)/(1+eEarth*cos(fE(i))); % Determines Orbital Radius For Given True Anomaly of Earth
        eT(j,i) = (DVmax(j)+(muSun*(2/RE-1/aEarth))^(1/2))^2*(RE/muSun)-1; % Determines Eccentricity of Transfer Orbit
        aT(j,i) = RE/(1-eT(j,i));
        Y = aT(j,i)*(1-eT(j,i)^2)/(aMars*(1-eMars^2));

        fT0 = pi/3; % Initilize True Anomaly for Loop
        fT1 = -pi/2; % Initilize True Anomaly for Loop
        
        % Loop That Solves the Transfer Orbit True Anomaly By Setting It's and Mars's Radii Equal
        
        while abs(fT1/fT0-1) > 0.0005 % Condition On Which to Run Loop
        %     if fE <= pi/2    
                if fT1 > pi
                    fT1 = fT1-pi;
                elseif fT1 < 0
                    fT1 = fT1+pi;
                else
                    fT0 = fT1;
                    fM(j,i) = fT0+fE(i)+2*pi-Theta;
                    fT1 = fT0-(1+eT(j,i)*cos(fT0)-Y-Y*eMars*cos(fM(j,i)))/(Y*eMars*sin(fM(j,i))-eT(j,i)*sin(fT0)); % Newton's Root Finding Method
                end
        %     elseif fE > pi/2
        %         if fT1 > pi
        %             fT1 = fT1-pi;
        %         elseif fT1 < 0
        %             fT1 = fT1+pi;
        %         else 
        %             fT0 = fT1;
        %             fM = fT0+fE+2*pi-Theta;
        %             fT1 = fT0-(1+eT*cos(fT0)-Y-Y*eMars*cos(fM))/(Y*eMars*sin(fM)-eT*sin(fT0)); % Newton's Root Finding Method
        %         end
        %     end
        end
        
        RT = aT(j,i)*(1-eT(j,i)^2)/(1+eT(j,i)*cos(fT1)); % Determines Orbit Radius of Transfer Orbit
        VTarrv = (muSun*(2/RT-1/aT(j,i)))^(1/2); % Determines Velocity Upon Arrival for Transfer Orbit
        RM = aMars*(1-eMars^2)/(1+eMars*cos(fM(j,i))); % Determines Orbit Radius of Mars
        VMarsArrv = (muSun*(2/RM-1/aMars))^(1/2); % Determines Velocity Upon Arrival for Mars
        phi = asin(aT(j,i)^2*(1-eT(j,i)^2)/(RT*(2*aT(j,i)-RT)));
        DVmaxArrv = (VTarrv^2+VMarsArrv^2-2*VTarrv*VMarsArrv*cos(phi-pi/2))^(1/2); % Determines Delta V to Orbit With Mars
        ET = 2*atan(((1-eT(j,i))/(1+eT(j,i)))^(1/2)*tan(fT1/2)); % Determines Transfer Orbit Eccentric Anomaly
        MT = ET-eT(j,i)*sin(ET);
        tT(j,i) = MT*(aT(j,i)^3/muSun)^(1/2);

        % Delta V to Low Mars Orbit With Patched Conics

        aPark = 330; % Altitude of Parking Orbit Above Mars [km]
        PerHyper = aPark + rMars; % Periapsis Distance From Center of Mars
        ah = 0; % Initializes Hyperbolic Semi-Major Axis
        eh = 0; % Initializes Hyperbolic Eccentricity

        ah = muMars/DVmaxArrv^2;
        eh = PerHyper/ah+1; % Determines Hyperbolic Eccentricity
        DVactualArrv2(j,i) = ((eh+1)/(eh-1))^(1/2)*DVmaxArrv-(muMars/PerHyper)^(1/2);
        
        % Loop for 2nd Mars Intersection (i.e. On Return Path of Transfer Orbit)
        
        fTRP0 = 3*pi/2; % Initilize True Anomaly Return Path for Loop
        fTRP1 = -pi; % Initilize True Anomaly Return Path for Loop
        
        while abs(fTRP1/fTRP0-1) > 0.0005 % Condition On Which to Run Loop
        %     if fE <= pi/2    
                if fTRP1 > 2*pi
                    fTRP1 = fTRP1-pi;
                elseif fTRP1 < pi
                    fTRP1 = fTRP1+pi;
                else
                    fTRP0 = fTRP1;
                    fMRP(j,i) = fTRP0+fE(i)+2*pi-Theta;
                    fTRP1 = fTRP0-(1+eT(j,i)*cos(fTRP0)-Y-Y*eMars*cos(fMRP(j,i)))/(Y*eMars*sin(fMRP(j,i))-eT(j,i)*sin(fTRP0)); % Newton's Root Finding Method
                end
        %     elseif fE > pi/2
        %         if fT1 > pi
        %             fT1 = fT1-pi;
        %         elseif fT1 < 0
        %             fT1 = fT1+pi;
        %         else 
        %             fT0 = fT1;
        %             fM = fT0+fE+2*pi-Theta;
        %             fT1 = fT0-(1+eT*cos(fT0)-Y-Y*eMars*cos(fM))/(Y*eMars*sin(fM)-eT*sin(fT0)); % Newton's Root Finding Method
        %         end
        %     end
        end
        
        RTRP = aT(j,i)*(1-eT(j,i)^2)/(1+eT(j,i)*cos(fTRP1)); % Determines Orbit Radius of Transfer Orbit
        VTarrvRP = (muSun*(2/RTRP-1/aT(j,i)))^(1/2); % Determines Velocity Upon Arrival for Transfer Orbit
        RMRP = aMars*(1-eMars^2)/(1+eMars*cos(fMRP(j,i))); % Determines Orbit Radius of Mars
        VMarsArrvRP = (muSun*(2/RMRP-1/aMars))^(1/2); % Determines Velocity Upon Arrival for Mars
        phi = asin(aT(j,i)^2*(1-eT(j,i)^2)/(RTRP*(2*aT(j,i)-RTRP)));
        DVmaxArrvRP = (VTarrvRP^2+VMarsArrvRP^2-2*VTarrvRP*VMarsArrvRP*cos(phi-pi/2))^(1/2); % Determines Delta V to Orbit With Mars
        ETRP = 2*atan(((1-eT(j,i))/(1+eT(j,i)))^(1/2)*tan(fTRP1/2)); % Determines Transfer Orbit Eccentric Anomaly
        
        if ETRP < 0 % Tests to Make Sure Eccentric Anomaly is Positive
            ETRP = ETRP + 2*pi;
        end
        
        MTRP = ETRP-eT(j,i)*sin(ETRP);
        tTRP(j,i) = MTRP*(aT(j,i)^3/muSun)^(1/2);

        % Delta V to Low Mars Orbit With Patched Conics

        aPark = 330; % Altitude of Parking Orbit Above Mars [km]
        PerHyper = aPark + rMars; % Periapsis Distance From Center of Mars
        ahRP = 0; % Initializes Hyperbolic Semi-Major Axis
        ehRP = 0; % Initializes Hyperbolic Eccentricity

        ahRP = muMars/DVmaxArrvRP^2;
        ehRP = PerHyper/ahRP+1; % Determines Hyperbolic Eccentricity
        DVactualArrv2RP(j,i) = ((ehRP+1)/(ehRP-1))^(1/2)*DVmaxArrvRP-(muMars/PerHyper)^(1/2);
        
        % Determines When To Increment Earth True Anomaly
        
        if i < M
            fE(i+1)=fE(i)+StepSize;
        end 
        
    end
    
    % Calculates Optimization Parameter
    
    for k = 1:M
        NuT = (tT(j,k)-min(tT(j,:)))/max(tT(j,:));
        NuDV = (DVactualArrv2(j,k)-min(DVactualArrv2(j,:)))/max(DVactualArrv2(j,:));
        Nut(j,k) = (NuT+NuDV)/2;
    end
   
    if j < N
        DVmax(j+1) = DVmax(j)-StepSizeN;
    end
end

% Finds Optimum Transfer

[val,ind]=min(Nut);
[MIN,ind2]=min(val);   
l = ind(ind2);
m = ind2;

% x = -2*pi:.1:2*pi;
% y = 1+eT*cos(x)-Y-Y*eMars*cos(x+fM-fT1);
% plot(x,y)


% subplot(2,3,1),[AX,H1,H2] = plotyy(fE,DVactualArrv2(1,:),fE,tT(1,:)/(TEarth),'plot'); % Plot of Arrival Delta V and Time of Transfer vs. Earth True Anomaly
% set(get(AX(1),'Ylabel'),'String','Arrival Delta V [km/s]')
% set(get(AX(2),'Ylabel'),'String','Time of Transfer [Earth Normalized]') 
% axis normal
% 
% subplot(2,3,2),[AX,H1,H2] = plotyy(fE,DVactualArrv2(5,:),fE,tT(5,:)/(TEarth),'plot'); % Plot of Arrival Delta V and Time of Transfer vs. Earth True Anomaly
% set(get(AX(1),'Ylabel'),'String','Arrival Delta V [km/s]')
% set(get(AX(2),'Ylabel'),'String','Time of Transfer [Earth Normalized]') 
% axis normal
% 
% subplot(2,3,3),[AX,H1,H2] = plotyy(fE,DVactualArrv2(10,:),fE,tT(10,:)/(TEarth),'plot'); % Plot of Arrival Delta V and Time of Transfer vs. Earth True Anomaly
% set(get(AX(1),'Ylabel'),'String','Arrival Delta V [km/s]')
% set(get(AX(2),'Ylabel'),'String','Time of Transfer [Earth Normalized]') 
% axis normal
% 
% subplot(2,3,5),[AX,H1,H2] = plotyy(fE,DVactualArrv2RP(ind(ind2),:),fE,tTRP(ind(ind2),:)/(TEarth),'plot'); % Plot of Arrival Delta V and Time of Transfer vs. Earth True Anomaly
% set(get(AX(1),'Ylabel'),'String','Arrival Delta V [km/s]')
% set(get(AX(2),'Ylabel'),'String','Time of Transfer [Earth Normalized]') 
% axis normal
% 
% subplot(2,3,6),[AX,H1,H2] = plotyy(fE,DVactualArrv2(ind(ind2),:),fE,tT(ind(ind2),:)/(TEarth),'plot'); % Plot of Arrival Delta V and Time of Transfer vs. Earth True Anomaly
% set(get(AX(1),'Ylabel'),'String','Arrival Delta V [km/s]')
% set(get(AX(2),'Ylabel'),'String','Time of Transfer [Earth Normalized]') 
% axis normal
% 
% % Creates Plot of Optimum Transfer
% 
% f = 0:(2*pi/998):2*pi;
% for k=1:999
%     r(k) = aT(l,m)*(1-eT(l,m)^2)/(1+eT(l,m)*cos(f(k)));
% end
% 
% for k=1:999
%     R(k) = aEarth*(1-eEarth^2)/(1+eEarth*cos(f(k)));
% end
% 
% for k=1:999
%     RR(k) = aMars*(1-eMars^2)/(1+eMars*cos(f(k)));
% end
% 
% subplot(2,3,4)
% hold on
% polar(f+fE(m),r,'r')
% polar(f,R)
% polar(f+Theta,RR)
% polar(fM(l,m)+Theta,aMars*(1-eMars^2)/(1+eMars*cos(fM(l,m))),'ro')
% polar(fMRP(l,m)+Theta,aMars*(1-eMars^2)/(1+eMars*cos(fMRP(l,m))),'ro')
% polar(fE(m),aEarth*(1-eEarth^2)/(1+eEarth*cos(fE(m))),'bo')
% rlim = aT(l,m)*(1+eT(l,m));
% axis([-1 1 -1 1]*rlim);
% axis square
% 
% disp('Optimum Transfer Total Delta V to Mars [km/s]')
% disp(DVactualArrv2(l,m)+DVmax(l))
% disp('Optimum Transfer Earth True Anomaly [Radians]')
% disp(fE(m))
% disp('Optimum Transfer Time [yrs]')
% disp(tT(l,m)/TEarth)

% % Creates Plot of Set Transfer
% 
% l = 4;
% m = 5;
% f = 0:(2*pi/998):2*pi;
% for k=1:999
%     r(k) = aT(l,m)*(1-eT(l,m)^2)/(1+eT(l,m)*cos(f(k)));
% end
% 
% for k=1:999
%     R(k) = aEarth*(1-eEarth^2)/(1+eEarth*cos(f(k)));
% end
% 
% for k=1:999
%     RR(k) = aMars*(1-eMars^2)/(1+eMars*cos(f(k)));
% end
% 
% subplot(2,3,5)
% hold on
% polar(f+fE(m),r,'r')
% polar(f,R)
% polar(f+Theta,RR)
% rlim = aT(l,m)*(1+eT(l,m));
% axis([-1 1 -1 1]*rlim);
% axis square

% Earth and Mars Actual Positions Beginning at Chosen Epoch

% fE0 = (4.54323446*10^(-4))*pi/180; fM0 = 319.4234417*pi/180; JDN0 = 2455565.2736; % Date of Epoch: 03 Jan 2011 18:34:00 UT
% fE0 = (7.8331757*10^(-1))*pi/180; fM0 = 199.18562*pi/180; JDN0 = 2457391.2736; % Date of Epoch: 03 Jan 2016 18:34:00 UT
% fE0 = (2.498736*10^(-1))*pi/180; fM0 = 90.5240655*pi/180; JDN0 = 2459217.2736; % Date of Epoch: 02 Jan 2021 18:34:00 UT
% fE0 = (7.558275*10^(-2))*pi/180; fM0 = 309.35329*pi/180; JDN0 = 2461044.2736; % Date of Epoch: 03 Jan 2026 18:34:00 UT
 fE0 = (3.5988321*10^(2))*pi/180; fM0 = 192.14974*pi/180; JDN0 = 2462871.2736; % Date of Epoch: 04 Jan 2031 18:34:00 UT

tmax = TEarth*4; % Set Max Time for Loop to Run For
t = 0; % Relative Time at Epoch

EE0 = 2*atan(tan(fE0/2)*((1-eEarth)/(1+eEarth))^(1/2)); % Earth Eccentric Anomaly at Epoch
EM0 = 2*atan(tan(fM0/2)*((1-eMars)/(1+eMars))^(1/2)); % Mars Eccentric Anomaly at Epoch
ME = zeros(1,int32(tmax/86400));
fEA = zeros(1,int32(tmax/86400));
MMA = zeros(1,int32(tmax/86400));
fMA = zeros(1,int32(tmax/86400));
fMM = zeros(1,int32(tmax/86400));
TrajNum = 0;
TrajNumOld = 0;
TrajNumER = 0;
TrajNumOldER = 0;
C = 0;
CER = 0;
v1 = 0;

if EM0 < 0 % Tests to Make Sure Eccentric Anomaly is Positive
    EM0 = EM0 + 2*pi;
end

MM0 = EM0-eMars*sin(EM0); % Mars Mean Anamoly at Epoch
TM0 = MM0/(2*pi/TMars); % Mars Period Advance at Epoch

if EE0 < 0 % Tests to Make Sure Eccentric Anomaly is Positive
    EE0 = EE0 + 2*pi;
end

ME0 = EE0-eEarth*sin(EE0); % Earth Mean Anamoly at Epoch
TE0 = ME0/(2*pi/TEarth); % Earth Period Advance at Epoch

i = 0;

while t < tmax

    % Increment Time and Counter
    
    t = t+86400;
    i = i+1;
    z = 0;    
% Earth Position Over Time
    
    ME(i) = (2*pi/TEarth)*((TE0+t)-fix((TE0+t)/TEarth)*TEarth); % Earth Mean Anomaly   
    
    % Iteration Loop to Find Earth Eccentric Anomaly
    
    EE0 = pi/3; % Initilize Eccentric Anomaly for Loop
    EE1 = -pi/2; % Initilize Eccentric Anomaly for Loop

    while abs(EE1-EE0) > 0.0005 % Condition On Which to Run Loop
        if EE1 > 2*pi
            EE1 = EE1-2*pi;
        elseif E1 < 0
            EE1 = EE1+2*pi;
        else
            EE0 = EE1;
            EE1 = EE0-(EE0-eEarth*sin(EE0)-ME(i))/(1-eEarth*cos(EE0)); % Newton's Root Finding Method
        end
    end
    
    fEA(i) = 2*atan(((1+eEarth)/(1-eEarth))^(1/2)*tan(EE1/2)); % Earth True Anomaly 
    
    % Test to Make Sure True Anomaly is Between 0 and 360 Degrees
    
    if fEA(i) < 0 
        fEA(i) = fEA(i)+2*pi;
    elseif fEA(i) > 2*pi;
        fEA(i) = fEA(i)-2*pi;
    end
    
% Mars Position Over Time

    MMA(i) = (2*pi/TMars)*((TM0+t)-fix((TM0+t)/TMars)*TMars); % Mars Mean Anomaly    
    
    % Iteration Loop to Find Mars Eccentric Anomaly
    
    EM0 = pi/3; % Initilize Eccentric Anomaly for Loop
    EM1 = -pi/2; % Initilize Eccentric Anomaly for Loop

    while abs(EM1-EM0) > 0.0005 % Condition On Which to Run Loop
        if EM1 > 2*pi
            EM1 = EM1-2*pi;
        elseif E1 < 0
            EM1 = EM1+2*pi;
        else
            EM0 = EM1;
            EM1 = EM0-(EM0-eMars*sin(EM0)-MMA(i))/(1-eMars*cos(EM0)); % Newton's Root Finding Method
        end
    end
    
    fMA(i) = 2*atan(((1+eMars)/(1-eMars))^(1/2)*tan(EM1/2)); % Mars True Anomaly 
    
    % Test to Make Sure True Anomaly is Between 0 and 360 Degrees
    
    if fMA(i) < 0 
        fMA(i) = fMA(i)+2*pi;
    elseif fMA(i) > 2*pi;
        fMA(i) = fMA(i)-2*pi;
    end
    
    % Mars Modified True Anomaly (i.e Mars Angular Position With Respect to Epoch)
    
    fMM(i) = fMA(i)-2*pi+Theta;
    
    if fMM(i) < 0 
        fMM(i) = fMM(i)+2*pi;
    elseif fMM(i) > 2*pi;
        fMM(i) = fMM(i)-2*pi;
    end
    
% Hohmann Transfer Options
    
    fMH(i) = fEA(i)-Theta+pi;
    RMH = aMars*(1-eMars^2)/(1+eMars*cos(fMH(i)));
    REH = aEarth*(1-eEarth^2)/(1+eEarth*cos(fEA(i)));
    aH = (RMH+REH)/2;
    eH = (RMH-REH)/(RMH+REH);
    TTH = pi*(aH^3/muSun)^(1/2);
    EMH = 2*atan(tan(fMH(i)/2)*((1-eMars)/(1+eMars))^(1/2)); % Mars Eccentric Anomaly for Hohmann Transfer
  
    if EMH < 0 % Tests to Make Sure Eccentric Anomaly is Positive
        EMH = EMH + 2*pi;
    end
    
    MMH = EMH-eMars*sin(EMH); % Mars Mean Anomaly On Arrival for Hohmann Transfer
    MTOA = MMH/(2*pi/TMars); % Mars Time On Arrival
    MTOED = MTOA-TTH; % Mars Time On Earth Departure
    
    if MTOED < 0
        MTOED = MTOED + TMars;
    end
    
    MMOED = (2*pi/TMars)*(MTOED); % Mars Mean Anomaly On Earth Departure
    
    if MMOED > 2*pi % Resets the Mean Anamoly if it Exceeds 2*pi
        MMOED = MMOED-2*pi;
    end
    
    % Iteration Loop to Find Mars Eccentric Anomaly
    
    EM0 = pi/3; % Initilize Eccentric Anomaly for Loop
    EM1 = -pi/2; % Initilize Eccentric Anomaly for Loop

    while abs(EM1-EM0) > 0.0005 % Condition On Which to Run Loop
        if EM1 > 2*pi
            EM1 = EM1-2*pi;
        elseif E1 < 0
            EM1 = EM1+pi;
        else
            EM0 = EM1;
            EM1 = EM0-(EM0-eMars*sin(EM0)-MMOED)/(1-eMars*cos(EM0)); % Newton's Root Finding Method
        end
    end
    
    fMOED(i) = 2*atan(((1+eMars)/(1-eMars))^(1/2)*tan(EM1/2)); % Mars True Anomaly On Earth Departure for Hohmann Transfer
    
    % Test to Make Sure True Anomaly is Between 0 and 360 Degrees
    
    if fMOED(i) < 0 
        fMOED(i) = fMOED(i)+2*pi;
    elseif fMOED(i) > 2*pi;
        fMOED(i) = fMOED(i)-2*pi;
    end
    
    k = 0;
    v2 = v1;
    
    if abs(fMOED(i)-fMA(i)) < pi/180
        k = k+1;
        v1 = i;
        ETOMA = t+TTH; % Earth Time On Mars Arrival With Respect To Perihelion
        MEOMA = (2*pi/TEarth)*ETOMA; % Earth Mean Anomaly of Earth On Mars Arrival
        fEOED(k) = fEA(i);
        
        % Earth Position On Mars Arrival
        
        if MEOMA > 2*pi % Resets the Mean Anamoly if it Exceeds 2*pi
            MEOMA = MEOMA-2*pi*(abs(round((t-TEarth)/TEarth))+1);
        end

        % Iteration Loop to Find Earth Eccentric Anomaly

        EE0 = pi/3; % Initilize Eccentric Anomaly for Loop
        EE1 = -pi/2; % Initilize Eccentric Anomaly for Loop

        while abs(EE1-EE0) > 0.0005 % Condition On Which to Run Loop
            if EE1 > 2*pi
                EE1 = EE1-2*pi;
            elseif E1 < 0
                EE1 = EE1+pi;
            else
                EE0 = EE1;
                EE1 = EE0-(EE0-eEarth*sin(EE0)-MEOMA)/(1-eEarth*cos(EE0)); % Newton's Root Finding Method
            end
        end

        fEOMA = 2*atan(((1+eEarth)/(1-eEarth))^(1/2)*tan(EE1/2)); % Earth True Anomaly On Mars Arrival

        % Test to Make Sure True Anomaly is Between 0 and 360 Degrees

        if fEOMA < 0 
            fEOMA = fEOMA+2*pi;
        elseif fEOMA > 2*pi;
            fEOMA = fEOMA-2*pi;
        end
        
        x1 = REH*cos(fEA(i));
        y1 = REH*sin(fEA(i));
        x2 = RMH*cos(fMM(i));
        y2 = RMH*sin(fMM(i));
        
        slope = (y2-y1)/(x2-x1);
        
        xPER = (slope*x1-y1)/(slope+1/slope);
        perDistance = xPER*(1+(-1/slope)^2)^(1/2);
        
        JDN = (t+JDN0*86400)/86400; % Computes Julian Day Number
        format long
        disp('Date of Hohmann Departure [Julian Day Number]')
        disp(JDN)
        format short
        disp('Iteration Number and Transfer Time [Days]')
        disp(i)
        disp(TTH/(86400))
        disp(fEA(i))
        disp(fMM(i))
        disp(fEOMA)
        disp(REH); disp(RMH)
        disp(perDistance)
      
        % Creates Plot of Transfer
        if (v1-v2) > 1
            f = 0:(2*pi/998):2*pi;
            eH =(RMH-REH)/(RMH+REH);
            %subplot(1,3,1)
            for k=1:999
                r(k) = aH*(1-eH^2)/(1+eH*cos(f(k)));
            end

            for k=1:999
                R(k) = aEarth*(1-eEarth^2)/(1+eEarth*cos(f(k)));
            end

            for k=1:999
                RR(k) = aMars*(1-eMars^2)/(1+eMars*cos(f(k)));
            end

            hold on
            polar(fMM(i),aMars*(1-eMars^2)/(1+eMars*cos(fMA(i))),'ro')
            polar(fEA(i),aEarth*(1-eEarth^2)/(1+eEarth*cos(fEA(i))),'bo')
            polar(fMH(i)+Theta,aMars*(1-eMars^2)/(1+eMars*cos(fMH(i))),'go')
            polar(fEOMA,aEarth*(1-eEarth^2)/(1+eEarth*cos(fEOMA)),'yo')
            polar(f+fEA(i),r,'g')
            polar(f,R,'b')
            polar(0,0,'kx')
            polar(f+Theta,RR,'r')
            rlim = aH*(1+eH)*1.25;
            axis([-1 1 -1 1]*rlim);
            axis square
        end
    end
    
% ************************ Determination of All Outbound Trajectories to Mars *****************************************************
  if t > TEarth*0
    disp(i)
    eT = 0.999;
    RAH = aMars*(1-eMars^2)/(1+eMars*cos(fEA(i)-Theta+pi));
    RPH = aEarth*(1-eEarth^2)/(1+eEarth*cos(fEA(i)));
    eH = (RAH-RPH)/(RAH+RPH);
    eStep = (eT-eH)/500; % Step size of Transfer Eccentricity
    n = 0; % Initiate Options Per Day Counter
    while eT >= eH
        
        aT = REH/(1-eT);
        Y = aT*(1-eT^2)/(aMars*(1-eMars^2));
        fT0 = pi/3; % Initilize True Anomaly for Loop
        fT1 = -pi/2; % Initilize True Anomaly for Loop
        
        while abs(fT1-fT0) > 0.0005 % Condition On Which to Run Loop
         
            if fT1 > pi
                fT1 = fT1-pi;
            elseif fT1 < 0
                fT1 = fT1+pi;
            else
                fT0 = fT1;
                fMOA = fT0+fEA(i)-Theta;
                fT1 = fT0-(1+eT*cos(fT0)-Y-Y*eMars*cos(fMOA))/(Y*eMars*sin(fMOA)-eT*sin(fT0)); % Newton's Root Finding Method
            end
        end  
         
        
        % Determination of Transfer Time to Mars's Orbit
        
        ET = 2*atan(((1-eT)/(1+eT))^(1/2)*tan(fT1/2)); % Determines Transfer Orbit Eccentric Anomaly
        
        if ET < 0
            ET = ET + 2*pi;
        end
        
        MT = ET-eT*sin(ET);
        tT = MT*(aT^3/muSun)^(1/2);
        
        % Determination of Mars's Time On Arrival
        
        EMOA = 2*atan(((1-eMars)/(1+eMars))^(1/2)*tan(fMOA/2)); % Determines Mars Eccentric Anomaly On Arrival
        
        if EMOA < 0
            EMOA = EMOA + 2*pi;
        end
        
        MMOA = EMOA-eMars*sin(EMOA);
        TMOA = MMOA*(aMars^3/muSun)^(1/2);
        
        % Determination of Mars's Time On Arrival
        
        EMA = 2*atan(((1-eMars)/(1+eMars))^(1/2)*tan(fMA(i)/2)); % Determines Mars Eccentric Anomaly On Earth Departure
        
        if EMA < 0
            EMA = EMA + 2*pi;
        end
        
        MMA = EMA-eMars*sin(EMA);
        TMA = MMA*(aMars^3/muSun)^(1/2);
        
        % Mars Period Advance From True Position and Hypothetical Transfer Orbit Intersection
        
        dTM = (TMOA-TMA);
        
        if dTM < 0
            dTM = dTM + TMars;
        end
        
   % Return Path Transfer Options to Mars
        
        fTRP0 = 3*pi/2; % Initilize True Anomaly for Loop
        fTRP1 = 2*pi/3; % Initilize True Anomaly for Loop
    
        while abs(fTRP1-fTRP0) > 0.0005 % Condition On Which to Run Loop

            if fTRP1 > 2*pi
                fTRP1 = fTRP1-pi;
            elseif fTRP1 < pi
                fTRP1 = fTRP1+pi;
            else
                fTRP0 = fTRP1;
                fMOARP = fTRP0+fEA(i)-Theta;
                fTRP1 = fTRP0-(1+eT*cos(fTRP0)-Y-Y*eMars*cos(fMOARP))/(Y*eMars*sin(fMOARP)-eT*sin(fTRP0)); % Newton's Root Finding Method
            end
        end  
        
        % Determination of Transfer Time to Mars's Orbit
        
        ETRP = 2*atan(((1-eT)/(1+eT))^(1/2)*tan(fTRP1/2)); % Determines Transfer Orbit Eccentric Anomaly
        if ETRP < 0
            ETRP = ETRP + 2*pi;
        end
        
        MTRP = ETRP-eT*sin(ETRP);
        tTRP = MTRP*(aT^3/muSun)^(1/2);
        
        % Determination of Mars's Time On Arrival
        
        EMOARP = 2*atan(((1-eMars)/(1+eMars))^(1/2)*tan(fMOARP/2)); % Determines Mars Eccentric Anomaly On Arrival
        
        if EMOARP < 0
            EMOARP = EMOARP + 2*pi;
        end
        
        MMOARP = EMOARP-eMars*sin(EMOARP);
        TMOARP = MMOARP*(aMars^3/muSun)^(1/2);
        
        % Determination of Mars's Time On Arrival
        
        EMA = 2*atan(((1-eMars)/(1+eMars))^(1/2)*tan(fMA(i)/2)); % Determines Mars Eccentric Anomaly On Earth Departure
        
        if EMA < 0
            EMA = EMA + 2*pi;
        end
        
        MMA = EMA-eMars*sin(EMA);
        TMA = MMA*(aMars^3/muSun)^(1/2);
        
        % Mars Period Advance From True Position and Hypothetical Transfer Orbit Intersection
        
        dTMRP = (TMOARP-TMA);
        
        if dTMRP < 0
            dTMRP = dTMRP + TMars;
        end
        
        % Test to See if Mars is in the Right Location for Transfer
        
        if abs(tT-dTM) < 0.5*86400
            
            z = 1;
            TrajNum = TrajNum+1; % Trajectory Number
            n = n+1; % Options Per Day Counter
            JDN = (t+JDN0*86400)/86400; % Computes Julian Day Number of Earth Departure
            JDNA = (t+tT+JDN0*86400)/86400; % Computes Julian Day Number of Mars Arrival
            
            fEOA=1;
            
            % The Following Uses the Method of Patched Conics to Determine Delta V
            
            DVTMI = (((1+eT)/(1-eT))*muSun/aT)^(1/2)-(muSun*(2/REH-1/aEarth))^(1/2); % Delta V for Trans-Martian Injection [km/s]
            aPark = 330; % Altitude of Parking Orbit Above Earth [km]
            PerHyper = aPark + rEarth; % Periapsis Distance From Center of Earth [km]
            ah = muEarth/DVTMI^2; % Determines Hyperbolic Semi-Major Axis
            eh = PerHyper/ah+1; % Determines Hyperbolic Eccentricity
            DVTMIactual = ((eh+1)/(eh-1))^(1/2)*DVTMI-(muEarth/PerHyper)^(1/2); % Actual Delta V for Trans-Martian Injection [km/s]

            % Delta V to Low Mars Orbit With Patched Conics
            
            RT = aT*(1-eT^2)/(1+eT*cos(fT1)); % Determines Orbit Radius of Transfer Orbit
            VTOA = (muSun*(2/RT-1/aT))^(1/2); % Determines Velocity On Arrival for Transfer Orbit
            VMOA = (muSun*(2/RT-1/aMars))^(1/2); % Determines Velocity On Arrival for Mars
            phiT = asin((aT^2*(1-eT^2)/(RT*(2*aT-RT)))^(1/2));
            phiM = asin((aMars^2*(1-eMars^2)/(RT*(2*aMars-RT)))^(1/2));
            DVMOI = (VTOA^2+VMOA^2-2*VTOA*VMOA*cos(phiM-phiT))^(1/2); % Delta V for Mars Orbit Insertion
            aPark = 150; % Altitude of Parking Orbit Above Mars [km]
            PerHyperM = aPark + rMars; % Periapsis Distance From Center of Mars
            ahM = muMars/DVMOI^2; % Determines Hyperbolic Semi-Major Axis
            ehM = PerHyperM/ahM+1; % Determines Hyperbolic Eccentricity
            DVMOIactual = ((ehM+1)/(ehM-1))^(1/2)*DVMOI-(muMars/PerHyperM)^(1/2); % Actual Delta V for Mars Orbit Insertion

            % Total DeltaV for Transit to Mars

            DVtotal = DVTMIactual+DVMOIactual;
            
            format long
            TrajectoryData(TrajNum,1) = JDN;
            TrajectoryData(TrajNum,2) = JDNA;
            TrajectoryData(TrajNum,3) = tT/86400;
            TrajectoryData(TrajNum,4) = DVTMIactual;
            TrajectoryData(TrajNum,5) = DVMOIactual;
            TrajectoryData(TrajNum,6) = DVtotal;
            TrajectoryData(TrajNum,7) = fEA(i);
            TrajectoryData(TrajNum,8) = fMA(i);
            TrajectoryData(TrajNum,9) = eT;
            TrajectoryData(TrajNum,10) = aT/aMars;
            %TrajectoryData(TrajNum,11) = fMM(i);
            %TrajectoryData(TrajNum,11) = REH;
            TrajectoryData(TrajNum,11) = i;
            TrajectoryData(TrajNum,12) = fMOA;
        end
        
        if abs(tTRP-dTMRP) < 0.5*86400
            
            z = 1;
            TrajNum = TrajNum+1; % Trajectory Number
            n = n+1; % Options Per Day Counter
            JDN = (t+JDN0*86400)/86400; % Computes Julian Day Number of Earth Departure
            JDNARP = (t+tTRP+JDN0*86400)/86400; % Computes Julian Day Number of Mars Arrival
            
            % The Following Uses the Method of Patched Conics to Determine Delta V
            
            DVTMI = (((1+eT)/(1-eT))*muSun/aT)^(1/2)-(muSun*(2/REH-1/aEarth))^(1/2); % Delta V for Trans-Martian Insertion [km/s]
            aPark = 330; % Altitude of Parking Orbit Above Earth [km]
            PerHyper = aPark + rEarth; % Periapsis Distance From Center of Earth [km]
            ah = muEarth/DVTMI^2; % Determines Hyperbolic Semi-Major Axis
            eh = PerHyper/ah+1; % Determines Hyperbolic Eccentricity
            DVTMIactual = ((eh+1)/(eh-1))^(1/2)*DVTMI-(muEarth/PerHyper)^(1/2); % Actual Delta V for Trans-Martian Insertion [km/s]

            % Delta V to Low Mars Orbit With Patched Conics
            
            RTRP = aT*(1-eT^2)/(1+eT*cos(fTRP1)); % Determines Orbit Radius of Transfer Orbit
            VTOARP = (muSun*(2/RTRP-1/aT))^(1/2); % Determines Velocity On Arrival for Transfer Orbit
            VMOARP = (muSun*(2/RTRP-1/aMars))^(1/2); % Determines Velocity On Arrival for Mars
            phiTRP = asin((aT^2*(1-eT^2)/(RTRP*(2*aT-RTRP)))^(1/2));
            phiMRP = asin((aMars^2*(1-eMars^2)/(RTRP*(2*aMars-RTRP)))^(1/2));
            DVMOIRP = (VTOARP^2+VMOARP^2-2*VTOARP*VMOARP*cos(phiMRP-phiTRP))^(1/2); % Delta V for Mars Orbit Injection
            aPark = 150; % Altitude of Parking Orbit Above Mars [km]
            PerHyperM = aPark + rMars; % Periapsis Distance From Center of Mars
            ahMRP = muMars/DVMOIRP^2; % Determines Hyperbolic Semi-Major Axis
            ehMRP = PerHyperM/ahMRP+1; % Determines Hyperbolic Eccentricity
            DVMOIactualRP = ((ehMRP+1)/(ehMRP-1))^(1/2)*DVMOIRP-(muMars/PerHyperM)^(1/2); % Actual Delta V for Mars Orbit Injection

            % Total DeltaV for Transit to Mars

            DVtotalRP = DVTMIactual+DVMOIactualRP;
            
            format long
            TrajectoryData(TrajNum,1) = JDN;
            TrajectoryData(TrajNum,2) = JDNARP;
            TrajectoryData(TrajNum,3) = tTRP/86400;
            TrajectoryData(TrajNum,4) = DVTMIactual;
            TrajectoryData(TrajNum,5) = DVMOIactualRP;
            TrajectoryData(TrajNum,6) = DVtotalRP;
            TrajectoryData(TrajNum,7) = fEA(i);
            TrajectoryData(TrajNum,8) = fMA(i);
            TrajectoryData(TrajNum,9) = eT;
            TrajectoryData(TrajNum,10) = aT/aMars;
            %TrajectoryData(TrajNum,11) = fMM(i);
            %TrajectoryData(TrajNum,11) = REH;
            TrajectoryData(TrajNum,11) = i;
            TrajectoryData(TrajNum,12) = fMOARP;
        end
        
        eT = eT-eStep;
        
    end
    
%     if TrajNum > TrajNumOld
%         C = C+1;
%         sum = zeros(n,12);
%         avg = zeros(n,12);
%         
%        for j = 1:12 
%            b = 0;
%            sumup = 0;
%         for ii = (TrajNum-(n-1)):TrajNum
%             b = b+1;
%             sumup = sumup+TrajectoryData(ii,j);
%             sum(b,j) = sumup; 
%             AvgTrajectoryData(C,j) = sum(b,j)/b;
%         end
%         
%         TrajNumOld = TrajNum;
%         
%        end
%     end

% Routine that Reduces the Amount of Trajectory Data to Three Trajectories Per Day

    if TrajNum > TrajNumOld
        if n > 3
            
            C = C+3;
            for j = 1:12 
                ReducedTrajectoryData(C-2,j) = TrajectoryData(TrajNum-(n-1),j);
                ReducedTrajectoryData(C-1,j) = TrajectoryData(round((TrajNum+TrajNum-(n-1))/2),j);
                ReducedTrajectoryData(C,j) = TrajectoryData(TrajNum,j);
            end       
         elseif n <= 3
               for ii = (TrajNum-(n-1)):TrajNum
                   C = C+1;
                   for j = 1:12
                       ReducedTrajectoryData(C,j) = TrajectoryData(ii,j);
                   end
               end

               TrajNumOld = TrajNum;
        end
     end
% ************************ Determination of All Inbound Trajectories to Earth *****************************************************
    
    fEH = fMA(i)+Theta-pi;
    RMHER = aMars*(1-eMars^2)/(1+eMars*cos(fMA(i)));
    REHER = aEarth*(1-eEarth^2)/(1+eEarth*cos(fEH));
    aHER = (RMHER+REHER)/2;
    eHER = (RMHER-REHER)/(RMHER+REHER);
    
    eTER = 0.5;
    eStepER = (eTER-eHER)/500; % Step size of Transfer Eccentricity Earth Return
    nER = 0; % Initiate Options Per Day Counter
    
    while eTER >= eHER
        
        aTER = RMHER/(1+eTER);
        Y = aTER*(1-eTER^2)/(aEarth*(1-eEarth^2));
        fT0 = 3*pi/2; % Initilize True Anomaly for Loop
        fT1 = -pi; % Initilize True Anomaly for Loop
        
        while abs(fT1-fT0) > 0.0005 % Condition On Which to Run Loop
         
            if fT1 > 2*pi
                fT1 = fT1-pi;
            elseif fT1 < pi
                fT1 = fT1+pi;
            else
                fT0 = fT1;
                fEER = fT0+fMA(i)+Theta-pi;
                fT1 = fT0-(1+eTER*cos(fT0)-Y-Y*eEarth*cos(fEER))/(Y*eEarth*sin(fEER)-eTER*sin(fT0)); % Newton's Root Finding Method
            end
        end  
         
        
        % Determination of Transfer Time to Earth's Orbit
        
        ETER = 2*atan(((1-eTER)/(1+eTER))^(1/2)*tan(fT1/2)); % Determines Transfer Orbit Eccentric Anomaly
        
        if ETER < 0
            ETER = ETER + 2*pi;
        end
        
        MTER = ETER-eTER*sin(ETER);
        tTER = (aTER^3/muSun)^(1/2)*(MTER-pi);
        
        % Determination of Earth's Time On Arrival
        
        EER = 2*atan(((1-eEarth)/(1+eEarth))^(1/2)*tan(fEER/2)); % Determines Earth Eccentric Anomaly On Arrival
        
        if EER < 0
            EER = EER + 2*pi;
        end
        
        MER = EER-eEarth*sin(EER);
        TER = MER*(aEarth^3/muSun)^(1/2);
        
        % Determination of Earth's Time On Mars Departure
        
        EEA = 2*atan(((1-eEarth)/(1+eEarth))^(1/2)*tan(fEA(i)/2)); % Determines Earth Eccentric Anomaly On Mars Departure
        
        if EEA < 0
            EEA = EEA + 2*pi;
        end
        
        MEA = EEA-eEarth*sin(EEA);
        TEA = MEA*(aEarth^3/muSun)^(1/2);
        
        % Earth Period Advance From True Position and Hypothetical Transfer Orbit Intersection
        
        dTMER = (TER-TEA);
        
        if dTMER < 0
            dTMER = dTMER + TEarth;
        end
        
    % Return Path Transfer Options to Earth
        
        fTRP0 = pi/3; % Initilize True Anomaly for Loop
        fTRP1 = -pi/2; % Initilize True Anomaly for Loop
        
        while abs(fTRP1-fTRP0) > 0.0005 % Condition On Which to Run Loop
         
            if fTRP1 > pi
                fTRP1 = fTRP1-pi;
            elseif fTRP1 < 0
                fTRP1 = fTRP1+pi;
            else
                fTRP0 = fTRP1;
                fEERRP = fTRP0+fMA(i)+Theta-pi;
                fTRP1 = fTRP0-(1+eTER*cos(fTRP0)-Y-Y*eEarth*cos(fEERRP))/(Y*eEarth*sin(fEERRP)-eTER*sin(fTRP0)); % Newton's Root Finding Method
            end
        end  
         
        
        % Determination of Transfer Time to Earth's Orbit
        
        ETERRP = 2*atan(((1-eTER)/(1+eTER))^(1/2)*tan(fTRP1/2)); % Determines Transfer Orbit Eccentric Anomaly
        
        if ETERRP < 0
            ETERRP = ETERRP + 2*pi;
        end
        
        MTERRP = ETERRP-eTER*sin(ETERRP);
        tTERRP = (MTERRP+pi)*(aTER^3/muSun)^(1/2);
        
        % Determination of Earth's Time On Arrival
        
        EERRP = 2*atan(((1-eEarth)/(1+eEarth))^(1/2)*tan(fEERRP/2)); % Determines Eartth Eccentric Anomaly On Arrival
        
        if EERRP < 0
            EERRP = EERRP + 2*pi;
        end
        
        MERRP = EERRP-eEarth*sin(EERRP);
        TERRP = MERRP*(aEarth^3/muSun)^(1/2);
        
        % Determination of Earth's Time On Mars Departure
        
        EEARP = 2*atan(((1-eEarth)/(1+eEarth))^(1/2)*tan(fEA(i)/2)); % Determines Earth Eccentric Anomaly On Mars Departure
        
        if EEARP < 0
            EEARP = EEARP + 2*pi;
        end
        
        MEARP = EEARP-eEarth*sin(EEARP);
        TEARP = MEARP*(aEarth^3/muSun)^(1/2);
        
        % Earth Period Advance From True Position and Hypothetical Transfer Orbit Intersection
        
        dTMERRP = (TERRP-TEARP);
        
        if dTMERRP < 0
            dTMERRP = dTMERRP + TEarth;
        end
    
        % Test to See if Mars is in the Right Location for Transfer
        
        if abs(tTER-dTMER) < 0.5*86400
            
            TrajNumER = TrajNumER+1; % Trajectory Number
            nER = nER+1; % Options Per Day Counter
            JDN = (t+JDN0*86400)/86400; % Computes Julian Day Number of Earth Departure
            JDNAER = (t+tTER+JDN0*86400)/86400; % Computes Julian Day Number of Mars Arrival
            
            fEOA=1;
            
            % The Following Uses the Method of Patched Conics to Determine Delta V
            
            DVTEI = (muSun*(2/RMHER-1/aMars))^(1/2)-(((1-eTER)/(1+eTER))*muSun/aTER)^(1/2); % Delta V for Trans-Earth Injection [km/s]
            aPark = 150; % Altitude of Parking Orbit Above Mars [km]
            PerHyper = aPark + rMars; % Periapsis Distance From Center of Mars [km]
            ah = muMars/DVTEI^2; % Determines Hyperbolic Semi-Major Axis
            eh = PerHyper/ah+1; % Determines Hyperbolic Eccentricity
            DVTEIactual = ((eh+1)/(eh-1))^(1/2)*DVTEI-(muMars/PerHyper)^(1/2); % Actual Delta V for Trans-Earth Injection [km/s]

            % Delta V to Low Earth Orbit With Patched Conics
            
            RTER = aTER*(1-eTER^2)/(1+eTER*cos(fT1)); % Determines Orbit Radius of Transfer Orbit
            VTOAER = (muSun*(2/RTER-1/aTER))^(1/2); % Determines Velocity On Arrival for Transfer Orbit
            VEOA = (muSun*(2/RTER-1/aEarth))^(1/2); % Determines Velocity On Arrival for Earth
            phiTER = asin((aTER^2*(1-eTER^2)/(RTER*(2*aTER-RTER)))^(1/2));
            phiEER = asin((aEarth^2*(1-eEarth^2)/(RTER*(2*aEarth-RTER)))^(1/2));
            DVEOI = (VTOAER^2+VEOA^2-2*VTOAER*VEOA*cos(phiTER-phiEER))^(1/2); % Delta V for Earth Orbit Insertion
            aPark = 330; % Altitude of Parking Orbit Above Mars [km]
            PerHyperE = aPark + rEarth; % Periapsis Distance From Center of Earth
            ahE = muEarth/DVEOI^2; % Determines Hyperbolic Semi-Major Axis
            ehE = PerHyperE/ahE+1; % Determines Hyperbolic Eccentricity
            DVEOIactual = ((ehE+1)/(ehE-1))^(1/2)*DVEOI-(muEarth/PerHyperE)^(1/2); % Actual Delta V for Earth Orbit Insertion

            % Total DeltaV for Transit to Earth

            DVtotalER = DVTEIactual+DVEOIactual;
            
            format long
            TrajectoryDataER(TrajNumER,1) = JDN;
            TrajectoryDataER(TrajNumER,2) = JDNAER;
            TrajectoryDataER(TrajNumER,3) = tTER/86400;
            TrajectoryDataER(TrajNumER,4) = DVTEIactual;
            TrajectoryDataER(TrajNumER,5) = DVEOIactual;
            TrajectoryDataER(TrajNumER,6) = DVtotalER;
            TrajectoryDataER(TrajNumER,7) = fEA(i);
            TrajectoryDataER(TrajNumER,8) = fMA(i);
            TrajectoryDataER(TrajNumER,9) = eTER;
            TrajectoryDataER(TrajNumER,10) = aTER/aEarth;
            TrajectoryDataER(TrajNumER,11) = fMM(i);
            TrajectoryDataER(TrajNumER,12) = fEER;
            
        end

        if abs(tTERRP-dTMERRP) < 0.5*86400
            
            TrajNumER = TrajNumER+1; % Trajectory Number
            nER = nER+1; % Options Per Day Counter
            JDN = (t+JDN0*86400)/86400; % Computes Julian Day Number of Earth Departure
            JDNAERRP = (t+tTERRP+JDN0*86400)/86400; % Computes Julian Day Number of Mars Arrival
            
            fEOA=1;
            
            % The Following Uses the Method of Patched Conics to Determine Delta V
            
            DVTEI = (muSun*(2/RMHER-1/aMars))^(1/2)-(((1-eTER)/(1+eTER))*muSun/aTER)^(1/2); % Delta V for Trans-Earth Insertion [km/s]
            aPark = 150; % Altitude of Parking Orbit Above Mars [km]
            PerHyper = aPark + rMars; % Periapsis Distance From Center of Mars [km]
            ah = muMars/DVTEI^2; % Determines Hyperbolic Semi-Major Axis
            eh = PerHyper/ah+1; % Determines Hyperbolic Eccentricity
            DVTEIactual = ((eh+1)/(eh-1))^(1/2)*DVTEI-(muMars/PerHyper)^(1/2); % Actual Delta V for Trans-Earth Insertion [km/s]

            % Delta V to Low Earth Orbit With Patched Conics
            
            RTERRP = aTER*(1-eTER^2)/(1+eTER*cos(fTRP1)); % Determines Orbit Radius of Transfer Orbit
            VTOAERRP = (muSun*(2/RTERRP-1/aTER))^(1/2); % Determines Velocity On Arrival for Transfer Orbit
            VEOARP = (muSun*(2/RTERRP-1/aEarth))^(1/2); % Determines Velocity On Arrival for Earth
            phiTERRP = asin((aTER^2*(1-eTER^2)/(RTERRP*(2*aTER-RTERRP)))^(1/2));
            phiEERRP = asin((aEarth^2*(1-eEarth^2)/(RTERRP*(2*aEarth-RTERRP)))^(1/2));
            DVEOIRP = (VTOAERRP^2+VEOARP^2-2*VTOAERRP*VEOARP*cos(phiTERRP-phiEERRP))^(1/2); % Delta V for Earth Orbit Injection
            aPark = 330; % Altitude of Parking Orbit Above Mars [km]
            PerHyperE = aPark + rEarth; % Periapsis Distance From Center of Earth
            ahERP = muEarth/DVEOIRP^2; % Determines Hyperbolic Semi-Major Axis
            ehERP = PerHyperE/ahERP+1; % Determines Hyperbolic Eccentricity
            DVEOIactualRP = ((ehERP+1)/(ehERP-1))^(1/2)*DVEOIRP-(muEarth/PerHyperE)^(1/2); % Actual Delta V for Earth Orbit Injection

            % Total DeltaV for Transit to Earth

            DVtotalERRP = DVTEIactual+DVEOIactualRP;
            
            format long
            TrajectoryDataER(TrajNumER,1) = JDN;
            TrajectoryDataER(TrajNumER,2) = JDNAERRP;
            TrajectoryDataER(TrajNumER,3) = tTERRP/86400;
            TrajectoryDataER(TrajNumER,4) = DVTEIactual;
            TrajectoryDataER(TrajNumER,5) = DVEOIactualRP;
            TrajectoryDataER(TrajNumER,6) = DVtotalERRP;
            TrajectoryDataER(TrajNumER,7) = fEA(i);
            TrajectoryDataER(TrajNumER,8) = fMA(i);
            TrajectoryDataER(TrajNumER,9) = eTER;
            TrajectoryDataER(TrajNumER,10) = aTER/aEarth;
            TrajectoryDataER(TrajNumER,11) = fMM(i);
            TrajectoryDataER(TrajNumER,12) = fEERRP;
        end        
        
        eTER = eTER-eStep;
    end    
   
%     
%     if TrajNumER > TrajNumOldER
%         CER = CER+1;
%         sum = zeros(nER,12);
%         avg = zeros(nER,12);
%         
%        for j = 1:12 
%            b = 0;
%            sumup = 0;
%         for ii = (TrajNumER-(nER-1)):TrajNumER
%             b = b+1;
%             sumup = sumup+TrajectoryDataER(ii,j);
%             sum(b,j) = sumup; 
%             AvgTrajectoryDataER(CER,j) = sum(b,j)/b; % Average Daily Trajectory Data for Earth Return from Mars
%         end
%         TrajNumOldER = TrajNumER;
%        end
%     end

% Routine that Reduces the Amount of Trajectory Data to Three Trajectories Per Day

     if TrajNumER > TrajNumOldER
        if nER > 3

            CER = CER+3;
            for j = 1:12 
                ReducedTrajectoryDataER(CER-2,j) = TrajectoryDataER(TrajNumER-(nER-1),j);
                ReducedTrajectoryDataER(CER-1,j) = TrajectoryDataER(round((TrajNumER+TrajNumER-(nER-1))/2),j);
                ReducedTrajectoryDataER(CER,j) = TrajectoryDataER(TrajNumER,j);
            end       
         elseif nER <= 3
               for ii = (TrajNumER-(nER-1)):TrajNumER
                   CER = CER+1;
                   for j = 1:12
                       ReducedTrajectoryDataER(CER,j) = TrajectoryDataER(ii,j);
                   end
               end

               TrajNumOldER = TrajNumER;
        end
     end
      
   end 
end

% % Creates Plot of Transfer
% 
%     subplot(1,3,3)
%     for i = 1:2
%         f = 0:(2*pi/998):2*pi;
%         eH =(RMH-REH)/(RMH+REH);
%         b = i;
%         
%         for k=1:999
%             r(k) = AvgTrajectoryData(b,10)*aMars*(1-AvgTrajectoryData(b,9)^2)/(1+AvgTrajectoryData(b,9)*cos(f(k)));
%         end
%         
%         for k=1:999
%             R(k) = aEarth*(1-eEarth^2)/(1+eEarth*cos(f(k)));
%         end
%         
%         for k=1:999
%             RR(k) = aMars*(1-eMars^2)/(1+eMars*cos(f(k)));
%         end
%         
%         hold on
%         polar(AvgTrajectoryData(b,11),aMars*(1-eMars^2)/(1+eMars*cos(AvgTrajectoryData(b,8))),'ro')
%         polar(AvgTrajectoryData(b,7),aEarth*(1-eEarth^2)/(1+eEarth*cos(AvgTrajectoryData(b,7))),'bo')
%         polar(AvgTrajectoryData(b,12)+Theta,aMars*(1-eMars^2)/(1+eMars*cos(AvgTrajectoryData(b,12))),'go')
%         %polar(fEOMA,aEarth*(1-eEarth^2)/(1+eEarth*cos(fEOMA)),'yo')
%         polar(f+AvgTrajectoryData(b,7),r,'g')
%         polar(f,R,'b')
%         polar(f+Theta,RR,'r')
%         %rlim = AvgTrajectoryData(100,10)*aMars*(1+AvgTrajectoryData(100,9))*1.25;
%         rlim = aH*(1+eH)*1.25;
%         axis([-1 1 -1 1]*rlim);
%         axis square
%     end
% 

%      
% Creates Animation of The Actual Relative Positions of Earth and Mars on the Plane of the Eccliptic

% for i = 1 : length(fMM)
%     hold on
%     p1 = polar(fMM(i),aMars*(1-eMars^2)/(1+eMars*cos(fMA(i))),'ro'); % Note: Mars point of Perihelion is Adjusted, but is Normal True Anomaly is Used for the Radius
%     p2 = polar(fEA(i),aEarth*(1-eEarth^2)/(1+eEarth*cos(fEA(i))),'bo');
%     %set(p1,'erasemode','xor')
%     axis([-ApMars ApMars -ApMars ApMars]) % Sets the Axis Max and Min to Mars Aphelion
%     axis square % Makes Each Axis Equal in Length
%     hold off
%     drawnow
% end
% hold off


k = 0; % Roud Trip Counter
C = 0; % Trajectory Options Counter
for i= 1:length(ReducedTrajectoryDataER)
    for j = 1:length(ReducedTrajectoryData)
        if ReducedTrajectoryDataER(i,1)-ReducedTrajectoryData(j,2) > 0
            k = k+1;
            RoundTrip(k) = ReducedTrajectoryDataER(i,2)-ReducedTrajectoryData(j,1);
            if RoundTrip(k) < 720 && (RoundTrip(k)-ReducedTrajectoryData(j,3)-ReducedTrajectoryDataER(i,3)) > 30 && (RoundTrip(k)-ReducedTrajectoryData(j,3)-ReducedTrajectoryDataER(i,3)) < 100 && (ReducedTrajectoryData(j,6)+ReducedTrajectoryDataER(i,6)) < 40
                C = C+1;
                format long
                StayTime = RoundTrip(k)-ReducedTrajectoryData(j,3)-ReducedTrajectoryDataER(i,3);
                for l = 3:14
                    TrajectoryOptions(C,1) = RoundTrip(k);
                    TrajectoryOptions(C,2) = StayTime;
                    TrajectoryOptions(C,l) = ReducedTrajectoryData(j,l-2);
                    TrajectoryOptions(C,l+12) = ReducedTrajectoryDataER(i,l-2);
                end
%                 disp(RoundTrip(k))
%                 disp(StayTime)
%                 format long
%                 disp(AvgTrajectoryData(j,1))
%                 disp(AvgTrajectoryDataER(i,2))
            end
        end
    end
end

% Mass Ratios

Isp = 1000; % Specific Impusle [s]
MassRatios = zeros(length(TrajectoryOptions(:,1)),4);

for i = 1:length(TrajectoryOptions(:,1))
    MassRatios(i,1) = exp(TrajectoryOptions(i,6)/(Isp*gEarth));
    MassRatios(i,2) = exp(TrajectoryOptions(i,7)/(Isp*gEarth));
    MassRatios(i,3) = exp(TrajectoryOptions(i,18)/(Isp*gEarth));
    MassRatios(i,4) = exp(TrajectoryOptions(i,19)/(Isp*gEarth));
    MassRatios(i,5) = MassRatios(i,1)*MassRatios(i,2)*MassRatios(i,3)*MassRatios(i,4);
end

% Trajectory Optimization Mass Based

MassInitial = 100; % Percent Mass
drop1 = 1.6; % Percent Mass Drop
drop2 = 0.6 + 6.4; % Percent Mass Drop
drop3 = 0.5; % Percent Mass Drop 
FinalMass = zeros(1,length(TrajectoryOptions(:,1))); % Initilize Array

for i = 1:length(TrajectoryOptions(:,1))
    FinalMass(i)=(((MassInitial/MassRatios(i,1)-drop1)/MassRatios(i,2)-drop2)/MassRatios(i,3)-drop3)/MassRatios(i,4);
end

[Final,OptimalTrajectory2]=max(FinalMass);

% Trajectory Optimization Delta-V Based

minDVtotal = min(TrajectoryOptions(:,8)+TrajectoryOptions(:,20));
minStayTime = min(TrajectoryOptions(:,2));
minRoundTrip = min(TrajectoryOptions(:,1));
minDVnotEarth = min((TrajectoryOptions(:,7)+TrajectoryOptions(:,18)+TrajectoryOptions(:,19)));
Optimize = zeros(1,length(TrajectoryOptions(:,1)));

for i = 1:length(TrajectoryOptions(:,1))
    Optimize(i) = TrajectoryOptions(i,1)/minRoundTrip+(TrajectoryOptions(i,8)+TrajectoryOptions(i,20))/minDVtotal+TrajectoryOptions(i,2)/minStayTime+(TrajectoryOptions(i,19)+TrajectoryOptions(i,7)+TrajectoryOptions(i,18))/minDVnotEarth;
end

% for i = 1:length(TrajectoryOptions(:,1))
%     Optimize(i) = TrajectoryOptions(i,1)/minRoundTrip+(TrajectoryOptions(i,8)+TrajectoryOptions(i,20))/minDVtotal+TrajectoryOptions(i,2)/minStayTime;
% end

[OptimalParameter,OptimalTrajectory]=min(Optimize);

disp('Optimal Round-Trip Mission Total Delta V [km/s]')
disp(TrajectoryOptions(OptimalTrajectory,8)+TrajectoryOptions(OptimalTrajectory,20))
disp('Optimal Round-Trip Duration [Days]')
disp(TrajectoryOptions(OptimalTrajectory,1))
disp('Optimal Mars Stay-Time [Days]')
disp(TrajectoryOptions(OptimalTrajectory,2))
disp('Optimal Earth Departure Date [JDN]')
disp(TrajectoryOptions(OptimalTrajectory,3))
disp('Optimal Mars Departure Date [JDN]')
disp(TrajectoryOptions(OptimalTrajectory,15))
disp('Trajectory Number')
disp(OptimalTrajectory)
disp('Trajectory Number Mass Based')
disp(OptimalTrajectory2)

subplot(1,3,2)

f = 0:(2*pi/998):2*pi;
eH =(RMHER-REHER)/(RMHER+REHER);

for k=1:999
    r(k) = TrajectoryOptions(OptimalTrajectory,24)*aEarth*(1-TrajectoryOptions(OptimalTrajectory,23)^2)/(1+TrajectoryOptions(OptimalTrajectory,23)*cos(f(k)));
end
for k=1:999
    rr(k) = TrajectoryOptions(OptimalTrajectory,12)*aMars*(1-TrajectoryOptions(OptimalTrajectory,11)^2)/(1+TrajectoryOptions(OptimalTrajectory,11)*cos(f(k)));
end
for k=1:999
    R(k) = aEarth*(1-eEarth^2)/(1+eEarth*cos(f(k)));
end

for k=1:999
    RR(k) = aMars*(1-eMars^2)/(1+eMars*cos(f(k)));
end

hold on
polar(TrajectoryOptions(OptimalTrajectory,22)+Theta,aMars*(1-eMars^2)/(1+eMars*cos(TrajectoryOptions(OptimalTrajectory,22))),'ko')
polar(TrajectoryOptions(OptimalTrajectory,26),aEarth*(1-eEarth^2)/(1+eEarth*cos(TrajectoryOptions(OptimalTrajectory,26))),'ko')
polar(TrajectoryOptions(OptimalTrajectory,9),aEarth*(1-eEarth^2)/(1+eEarth*cos(TrajectoryOptions(OptimalTrajectory,9))),'go')
polar(TrajectoryOptions(OptimalTrajectory,14)+Theta,aMars*(1-eMars^2)/(1+eMars*cos(TrajectoryOptions(OptimalTrajectory,14))),'go')
polar(f+TrajectoryOptions(OptimalTrajectory,22)-pi+Theta,r,'y'),polar(f+TrajectoryOptions(OptimalTrajectory,9),rr,'g')
polar(f,R,'b')
polar(f+Theta,RR,'r')
rlim = aH*(1+eH)*1.25;
axis([-1 1 -1 1]*rlim);
axis square

% figure
% for i = 7:7
%         f = 0:(2*pi/998):2*pi;
%         eH =(RMHER-REHER)/(RMHER+REHER);
%         b = i;
% 
%         for k=1:999
%             r(k) = TrajectoryOptions(b,24)*aEarth*(1-TrajectoryOptions(b,23)^2)/(1+TrajectoryOptions(b,23)*cos(f(k)));
%         end
%         for k=1:999
%             rr(k) = TrajectoryOptions(b,12)*aMars*(1-TrajectoryOptions(b,11)^2)/(1+TrajectoryOptions(b,11)*cos(f(k)));
%         end
%         for k=1:999
%             R(k) = aEarth*(1-eEarth^2)/(1+eEarth*cos(f(k)));
%         end
% 
%         for k=1:999
%             RR(k) = aMars*(1-eMars^2)/(1+eMars*cos(f(k)));
%         end
% 
%         hold on
%         polar(TrajectoryOptions(b,22)+Theta,aMars*(1-eMars^2)/(1+eMars*cos(TrajectoryOptions(b,22))),'ko')
%         polar(TrajectoryOptions(b,26),aEarth*(1-eEarth^2)/(1+eEarth*cos(TrajectoryOptions(b,26))),'ko')
%         polar(TrajectoryOptions(b,9),aEarth*(1-eEarth^2)/(1+eEarth*cos(TrajectoryOptions(b,9))),'go')
%         polar(TrajectoryOptions(b,14)+Theta,aMars*(1-eMars^2)/(1+eMars*cos(TrajectoryOptions(b,14))),'go')
%         polar(f+TrajectoryOptions(b,22)-pi+Theta,r,'y'),polar(f+TrajectoryOptions(b,9),rr,'g')
%         polar(f,R,'b')
%         polar(f+Theta,RR,'r')
%         rlim = aH*(1+eH)*1.25;
%         axis([-1 1 -1 1]*rlim);
%         axis square
% end
    
subplot(1,3,1)
X=1:length(TrajectoryOptions(:,2));
hold on
scatter(TrajectoryOptions(X,2),(TrajectoryOptions(X,8)+TrajectoryOptions(X,20)),'bo');
scatter(TrajectoryOptions(OptimalTrajectory,2),(TrajectoryOptions(OptimalTrajectory,8)+TrajectoryOptions(OptimalTrajectory,20)),'ro');
xlabel('Mars Stay-Time [Days]')
ylabel('Round-Trip Total Delta V [km/s]')
hold off
subplot(1,3,3)
hold on
scatter(TrajectoryOptions(X,2),(TrajectoryOptions(X,1)+TrajectoryOptions(X,20)),'bo');
scatter(TrajectoryOptions(OptimalTrajectory,2),(TrajectoryOptions(OptimalTrajectory,1)),'ro');
xlabel('Mars Stay-Time [Days]')
ylabel('Round-Trip Duration [Days]')
hold off

% Launch Window for Earth Departure

i=1;
while TrajectoryData(i,1) < TrajectoryOptions(OptimalTrajectory,3)
    i = i+1;
end

j = i+1;
k = 0;

while (TrajectoryData(j,1)-TrajectoryData(j-1,1)) == 1 || (TrajectoryData(j,1)-TrajectoryData(j-1,1)) == 0
    k = k+1;
    for ii = 1:length(TrajectoryData(1,:))
        LaunchWindow(k,ii) = TrajectoryData(j-1,ii);
    end
    j = j+1;
end

figure
hold on
scatter(LaunchWindow(:,1),LaunchWindow(:,5),'ro')
scatter(LaunchWindow(:,1),LaunchWindow(:,4),'bo')
scatter(LaunchWindow(1,1),LaunchWindow(1,5),'ko')
scatter(LaunchWindow(1,1),LaunchWindow(1,4),'ko')
xlabel('Date of Earth Departure [JD]')
ylabel('Impusle Delta-V [km/s] (at Earth: blue, at Mars: red)')

% Launch Window for Mars Departure

i=1;
while TrajectoryDataER(i,1) < TrajectoryOptions(OptimalTrajectory,15)
    i = i+1;
end

j = i+1;
jm = j-1;
k = 0;

while (TrajectoryDataER(j,1)-TrajectoryDataER(jm,1)) == 1 || (TrajectoryDataER(j,1)-TrajectoryDataER(jm,1)) == 0
    k = k+1;
    for ii = 1:length(TrajectoryDataER(1,:))
        LWER(k,ii) = TrajectoryDataER(jm,ii);
    end
    j = j+1;
    jm = j-1;
    if j > length(TrajectoryDataER(:,1))
        j = 2;
        jm = length(TrajectoryDataER(:,1));
    end
end
kk = k;
for ii = 1:k % Loop Reverses Row Order of Matrix
    for jj = 1:length(LWER(1,:))
        LaunchWindowER(ii,jj) = LWER(k+1-ii,jj);
    end
end

j = i-1;

while j > 0 && ((TrajectoryDataER(j+1,1)-TrajectoryDataER(j,1)) == 1 || (TrajectoryDataER(j+1,1)-TrajectoryDataER(j,1)) == 0)
    k = k+1;
    for ii = 1:length(TrajectoryDataER(1,:))
        LaunchWindowER(k,ii) = TrajectoryDataER(j+1,ii);
    end
    j = j-1;
end


figure
hold on
scatter(LaunchWindowER(:,1),LaunchWindowER(:,5),'bo')
scatter(LaunchWindowER(:,1),LaunchWindowER(:,4),'ro')
scatter(LaunchWindowER(kk,1),LaunchWindowER(kk,5),'ko')
scatter(LaunchWindowER(kk,1),LaunchWindowER(kk,4),'ko')
xlabel('Date of Mars Departure [JD]')
ylabel('Impusle Delta-V [km/s] (at Earth: blue, at Mars: red)')

% Plots Selected Trajectory From Launch Windows

% figure 
% for i = 1122:1122
%         f = 0:(2*pi/998):2*pi;
%         eH =(RMHER-REHER)/(RMHER+REHER);
%         b = i;
% 
% %         for k=1:999
% %            r(k) = LaunchWindowER(b,10)*aEarth*(1-LaunchWindowER(b,9)^2)/(1+LaunchWindowER(b,9)*cos(f(k)));
% %         end
%         for k=1:999
%             rr(k) = LaunchWindow(b,10)*aMars*(1-LaunchWindow(b,9)^2)/(1+LaunchWindow(b,9)*cos(f(k)));
%         end
%         for k=1:999
%             R(k) = aEarth*(1-eEarth^2)/(1+eEarth*cos(f(k)));
%         end
% 
%         for k=1:999
%             RR(k) = aMars*(1-eMars^2)/(1+eMars*cos(f(k)));
%         end
% 
%         hold on
%         %polar(LaunchWindowER(b,8)+Theta,aMars*(1-eMars^2)/(1+eMars*cos(LaunchWindowER(b,8))),'ko')
%         %polar(LaunchWindowER(b,12),aEarth*(1-eEarth^2)/(1+eEarth*cos(LaunchWindowER(b,12))),'ko')
%         polar(LaunchWindow(b,7),aEarth*(1-eEarth^2)/(1+eEarth*cos(LaunchWindow(b,7))),'go')
%         polar(LaunchWindow(b,12)+Theta,aMars*(1-eMars^2)/(1+eMars*cos(LaunchWindow(b,12))),'go')
%         %polar(f+LaunchWindowER(b,8)-pi+Theta,r,'y')
%         polar(f+LaunchWindow(b,7),rr,'g')
%         polar(f,R,'b')
%         polar(f+Theta,RR,'r')
%         rlim = aH*(1+eH)*1.25;
%         axis([-1 1 -1 1]*rlim);
%         axis square
% end

% figure 
% for i = 1122:1122
%         f = 0:(2*pi/998):2*pi;
%         eH =(RMHER-REHER)/(RMHER+REHER);
%         b = i;
% 
%         for k=1:999
%            r(k) = LaunchWindowER(b,10)*aEarth*(1-LaunchWindowER(b,9)^2)/(1+LaunchWindowER(b,9)*cos(f(k)));
%         end
%  
%         for k=1:999
%             R(k) = aEarth*(1-eEarth^2)/(1+eEarth*cos(f(k)));
%         end
% 
%         for k=1:999
%             RR(k) = aMars*(1-eMars^2)/(1+eMars*cos(f(k)));
%         end
% 
%         hold on
%         polar(LaunchWindowER(b,8)+Theta,aMars*(1-eMars^2)/(1+eMars*cos(LaunchWindowER(b,8))),'ko')
%         polar(LaunchWindowER(b,12),aEarth*(1-eEarth^2)/(1+eEarth*cos(LaunchWindowER(b,12))),'ko')
%         polar(f+LaunchWindowER(b,8)-pi+Theta,r,'y')
%         polar(f,R,'b')
%         polar(f+Theta,RR,'r')
%         rlim = aH*(1+eH)*1.25;
%         axis([-1 1 -1 1]*rlim);
%         axis square
% end

% Output .txt File 

fid = fopen('TrajectoryOptions.txt', 'w');

fprintf(fid, 'Trajectory Otions\r\n\r\n');
fprintf(fid, 'Number  Round Trip  Stay Time  Departure Date  Arrival Date  Transit Time  DeltaV Earth  DeltaV Mars  Total DeltaV  Earth f    Mars f     Transfer e  Transfer a   Perihelion      fMOA       Departure Date  Arrival Date  Transit Time  DeltaV Mars  DeltaV Earth  Total DeltaV  Earth f    Mars f     Transfer e  Transfer a \r\n\r\n');
fprintf(fid, '        [Days]      [Days]     [JDN]           [JDN]         [Days]        [km/s]        [km/s]       [km/s]        [radians]  [radians]  [unitless]  [Mars-norm]  [km]            [radians]  [JDN]           [JDN]         [Days]        [km/s]        [km/s]       [km/s]        [radians]  [radians]  [unitless]  [Mars-norm]\r\n\r\n'); 

j = 0;
for i = 1:length(TrajectoryOptions(:,1))
    j = j+1;
fprintf(fid, '%6.0f  %10.4f  %9.4f  %14.4f  %12.4f  %12.4f  %12.4f  %11.4f  %12.4f  %9.4f  %9.4f  %10.4f  %11.4f  %14.4f  %9.4f  %14.4f  %12.4f  %12.4f  %11.4f  %12.4f  %12.4f  %9.4f  %9.4f  %10.4f  %11.4f  %.4f  %.4f\r\n',j,TrajectoryOptions(i,:));
end
fclose(fid);    

%
fid = fopen('LaunchWindowEarth.txt', 'w');

fprintf(fid, 'Launch Window Earth\r\n\r\n');
fprintf(fid, 'Number  Departure Date  Arrival Date  Transit Time  DeltaV Earth  DeltaV Mars  Total DeltaV  Earth f    Mars f     Transfer e  Transfer a   Perihelion      fMOA     \r\n\r\n');
fprintf(fid, '        [JDN]           [JDN]         [Days]        [km/s]        [km/s]       [km/s]        [radians]  [radians]  [unitless]  [Mars-norm]  [km]            [radians]\r\n\r\n'); 

j = 0;
for i = 1:length(LaunchWindow(:,1))
    j = j+1;
fprintf(fid, '%6.0f  %14.4f  %12.4f  %12.4f  %12.4f  %11.4f  %12.4f  %9.4f  %9.4f  %10.4f  %11.4f  %14.4f  %9.4f\r\n',j,LaunchWindow(i,:));
end
fclose(fid);

%
fid = fopen('LaunchWindowMars.txt', 'w');

fprintf(fid, 'Launch Window Mars\r\n\r\n');
fprintf(fid, 'Number  Departure Date  Arrival Date  Transit Time  DeltaV Mars  DeltaV Earth  Total DeltaV  Earth f    Mars f     Transfer e  Transfer a  \r\n\r\n');
fprintf(fid, '        [JDN]           [JDN]         [Days]        [km/s]        [km/s]       [km/s]        [radians]  [radians]  [unitless]  [Earth-norm]\r\n\r\n'); 

j = 0;
for i = 1:length(LaunchWindowER(:,1))
    j = j+1;
fprintf(fid, '%6.0f  %14.4f  %12.4f  %12.4f  %11.4f  %12.4f  %12.4f  %9.4f  %9.4f  %10.4f  %11.4f  %.4f  %.4f\r\n',j,LaunchWindowER(i,:));
end
fclose(fid);

toc