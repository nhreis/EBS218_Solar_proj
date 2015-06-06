%%% Matlab program to solve transient flat-plate air heater model which is S shaped with fin
clear all
clc
%% Find EoT
% constants
rho_g = 0.20;   % diffuse reflectance of ground
tilt = 12;    % tilt angle of solar collector (deg)
lat = 38;   % latitude (deg North)
time = [1:24];  % time of day (24 hr clock)
gamma = 0;  % surface azimuth angle (zero due south)
pos=[lat, 121];
day=148;
SurAz=0;
%local meridian for Davis is 120 degrees
merid = 120;
%Calculate Solar Declination Angle
SolDec = 23.45*sind(360*(284+day)/365);
%Calculate solar time (solT)
%find Equation of Time (EoT)
B=(day -1)*(360/365);
EoT = 229.2*(0.000075+0.001868*cosd(B) - 0.032077*sind(B) - 0.014615*cosd(2*B) - 0.04089*sind(2*B));

%% Read data to get Gb , Gd , Wind speed & ambient air temperature
% read variables from file
davis = csvread('Solar2014csv.csv');    % reads pyranometer data from CSV file
%weather = xlsread('WeatherStationData2.csv' ); % reads weather data

%wind = weather(1:length(weather),8);    % wind speed column
%wind_avg = mean(wind);  % average wind speed
wind_avg=1.06;
Gt_vert = davis(3456:3479,10);  % south vertical pyranometer (W/m2)
G_horz = davis(3456:3479,7);    % horizontal pyranometer (W/m2)
air_Temp = davis(3456:3479,4);  % air temp (deg C)
day = davis(3456,2);    % day of year


% calculate zenith and incidence angle for pyranometers
beta = 90;  % slope angle (vertical surface)
omega = 15*(time-12);   % hour angle (deg)
decl = 23.45*sind(360*(284+day)/365);    % declination angle (deg)
theta_z = acosd(cosd(lat)*cosd(decl)*cosd(omega)+sind(lat)*sind(decl));   % zenith angle
theta = acosd(sind(decl)*sind(lat)*cosd(beta)...
    -sind(decl)*cosd(lat)*sind(beta)*cosd(gamma)...
    +cosd(decl)*cosd(lat)*cosd(beta)*cosd(omega)...
    +cosd(decl)*sind(lat)*sind(beta)*cosd(gamma)*cosd(omega)...
    +cosd(decl)*sind(beta)*sind(gamma)*sind(omega));    % incidence angle

% calculate G_d and G_b
p1 = (cosd(theta)./cosd(theta_z))+rho_g.*((1-cosd(beta))/2);
p2 = ((1+cosd(beta))./2)-(cosd(theta)./cosd(theta_z));

G_d = (Gt_vert - G_horz.*(p1'))./p2'; % diffuse radiation on horizontal surface
G_b = G_horz - G_d; % beam radiation on horizontal surface
G_d(G_d<0)=0;   % replaces negative values with zero
G_b(G_b<0)=0;   % "                                 "


%% Define Variables
% Solves equations from ebs 218 lecture notes using explicit finite
% difference method.

Starttime=1; %When the model begins
sigma=5.67e-8;
DuctHeight=.058;
DuctWidth=0.33;
DuctLength=3;
Plate_emit=.98;
Back_emit=.98;

Initial_AirT = air_Temp(Starttime); %

%taualfa=.85;

Inlet_T=air_Temp(Starttime); % Inlet temperature equals to ambient air temp at start time

vw = wind_avg; % Wind speed

CoverN=1;
Cover_emit=.88;
Back_loss=0;
% Absorber thermal properties for Alumin Alloy 195, Table A-3,
% pg 948, Heat Transfer: A Practical Approach by Y.A. Cengel
% McGraw-Hill Book Co., New York.
Plate_rho=2700;
Plate_specificH=904;
Plate_thick=.0015875;
Plate_Conduct=235;
Back_rho=30;
Back_specificH=2000;
Back_thick=0.0254;
Back_Conduct=0.00001;


FlowRate=0.061;

%% parameters for finite difference solution
dt=3;
Nseg=60;
tstop=23*3600;
ksav=1;
isav=0;
isv=1;
tyme=0;
% calculate absolute temperatures
Air_Tk=Initial_AirT+273.15;
% set inital values for air flow and solar radiation
amdot=FlowRate;

% get air property values from curve fits of data tables
[denf,therm_diff,amu,kin_vis,akf,FlowSpec,prf]=fluid_props_air(Inlet_T+273.15);
% collector areas and hydraulic diameter
ac=DuctWidth*DuctLength;
ad=DuctHeight*DuctWidth;
dh=2*DuctHeight;
% variables needed for finite difference calculations
np=Nseg+1;
npm1=np-1;
dx=DuctLength/Nseg;
x=linspace(0,DuctLength,np);
c1=dt/(Plate_rho*Plate_specificH*Plate_thick);
c2=Plate_Conduct*Plate_thick/(dx*dx);
c3=dt/(Back_rho*Back_specificH*Back_thick);
c4=Back_Conduct*Back_thick/(dx*dx);
c5=amdot*FlowSpec/(DuctWidth*dx);

%% initialize temperature arrays
for i=1:np
    tp(i,1)=Initial_AirT+.01;
    tp(i,2)=Initial_AirT+.01;
    tb(i,1)=Initial_AirT;
    tb(i,2)=Initial_AirT;
    t(i)=Initial_AirT+.001;
end

Gb=G_b(Starttime); % Get Gb and Gd
Gd=G_d(Starttime);
% Get sun position
[ SolZen , SolAz, SolAlt, hrAng, indAng ] = SunPosition( pos, Starttime, SolDec, tilt, SurAz );
    
%Calculate solar energy absorbed by the plate
SolarEnergy=Transimission(Gb,Gd,tilt,indAng);

% initial value for solar radiation
% forced convection due to wind  eqn 3.15.3, page 174 d&b
hw=2.8 + 3*vw;
% parameters for top loss coefficient equation (6.4.9), page 260 d&b
if(tilt>=70)
    tilt=70;
end
c=520*(1-0.000051*tilt^2);
f=(1+0.089*hw-0.1166*hw*Plate_emit)*(1+0.07866*CoverN);
term1=1/(Plate_emit+0.00591*CoverN*hw);
term2=(2*CoverN+f-1+0.133*Plate_emit)/Cover_emit;
t4=term1+term2-CoverN;
% set initial values
tpave=(tp(1,2)+tp(np,2))/2;
tbave=(tb(1,2)+tb(np,2))/2;
tfave=(t(1)+t(np))/2;
savtim(ksav)=tyme;
savtp(:,ksav)=tp(:,2);
savt(:,ksav)=t(:);
savt(ksav)=t(np);
savgt(ksav)=SolarEnergy;
saveff(ksav)=0;
savamd(ksav)=amdot;
ih=0;
im=0;

xlswrite(FinalData,savtime,1)
xlswrite(FinalData,savtp,2)





%% time loop
for tyme=dt:dt:tstop
    isav=isav+1;
    
    clocktime=Starttime+tyme/3600; % Calculate the corresponding clock time
    clock=floor(clocktime);
    solT = 4*(merid-pos(2))+ EoT + clocktime*60;
    Gb=G_b(clock); % Get Gb and Gd
    Gd=G_d(clock);
    % Get sun position
    [ SolZen , SolAz, SolAlt, hrAng, indAng ] = SunPosition( pos, solT, SolDec, tilt, SurAz );
    
    %Calculate solar energy absorbed by the plate
    SolarEnergy=Transimission(Gb,Gd,tilt,indAng);
    %SolarEnergy=gt*taualfa;
    % calculate hrpb, radiation heat trassfer coefficient
    tpk=tpave+273.15;
    tbk=tbave+273.15;
    hrpb=sigma*(tpk+tbk)*(tpk*tpk+tbk*tbk)/(1/Plate_emit+1/Back_emit-1);
    % get air property values from curve fits of data tables
    tfk=tfave+273.15;
    [denf,therm_diff,amu,kin_vis,akf,FlowSpec,prf]=fluid_props_air(tfk);
    % reynolds number in duct (see example 3.14.2, p173 d&b)
    re=2*amdot/(DuctWidth*amu);
    if(re<2100)
        % laminar flow  convective heat transfer coefficient (eqn 3.14.7, p172, d&b)
        fac=re*prf*dh/DuctLength;
        anud=4.9 + .0606*fac^1.2/(1+0.0909*(fac^0.7)*(prf^0.17));
    else
        % turbulent convective heat transfer coefficient (eqn 3.14.6, p171, d&b)
        anud=0.0158*(re^0.8);
    end
    hcpf=anud*akf/dh;
    hcbf=hcpf;
    % collector top heat loss coefficient, ut, use empirical correlation
    % eqn (6.4.9) on page 260 of d&b.
    et=0.430*(1-100/tpk);
    t1=CoverN/((c/tpk)*((tpk-Air_Tk)/(CoverN+f))^et);
    t2=1/(t1+1/hw);
    t3=sigma*(tpk+Air_Tk)*(tpk*tpk+Air_Tk*Air_Tk);
    if tpk>Air_Tk
        utop=t2+t3/t4
    else
        utop=0;
    end
    ul=utop+Back_loss;
    
    %Fin
    FinThickness = Plate_thick;
    Ac = FinThickness * dx;
    P = 2*dx + 2*FinThickness;
    m=(hcpf * P / ( Plate_Conduct * Ac) )^(1/2);
    qcoef=(hcpf * P * Plate_Conduct * Ac )^(1/2) * tand( m * DuctHeight );
    
    
    % finite difference equations to calculate temperatures
    % distance loop
    Inlet_T=air_Temp(clock)
    t(1)=Inlet_T;
    for i=2:npm1
        if (dx*i <=0.67) || (dx*i>=1 && dx*i<=1.33)||(dx*i>=2.33) % One side fin
            N_fin=1;
        elseif (dx*i>=1.33&&dx*i<=1.67) % Two side fin
            N_fin=2;
        else % No side fin
            N_fin=0;
        end
        
        tp(i,2)=tp(i,1) + c1*(SolarEnergy+c2*(tp(i+1,1)-2.*tp(i,1)...
            +tp(i-1,1))-hcpf*(tp(i,1)-t(i))-hrpb*(tp(i,1)...
            -tb(i,1))-utop*(tp(i,1)-Initial_AirT)...
            - N_fin * qcoef * ( tp(i,1) - t(i)) );
        
        tb(i,2)=tb(i,1) + c3*(c4*(tb(i+1,1)-2.*tb(i,1)...
            +tb(i-1,1))-hcbf*(tb(i,1)-t(i))+hrpb*(tp(i,1)...
            -tb(i,1))-Back_loss*(tb(i,1)-Initial_AirT));
        
        t(i)=(hcpf*tp(i,2)+hcbf*tb(i,2)+c5*t(i-1) ...
            + N_fin * qcoef * ( tp(i,1) - t(i) ))/(hcpf+hcbf+c5);
    end  % of distance loop
    % boundary conditions
    tp(1,2)=tp(2,2);
    tp(np,2)=tp(npm1,2);
    tb(1,2)=tb(2,2);
    tb(np,2)=tb(npm1,2);
    t(np)=(hcpf*tp(np,2)+hcbf*tb(np,2)+c5*t(npm1))/(hcpf+hcbf+c5);
    % average temperatures
    tpave=(tp(1,2)+tp(np,2))/2;
    tbave=(tb(1,2)+tb(np,2))/2;
    tfave=(t(1)+t(np))/2;
    % useful heat and efficiency
    qu=FlowRate*FlowSpec*(t(np)-t(1));
    eff=100.*qu/(SolarEnergy*ac);
    % save some data for plotting
    if(isav>=isv)
        isav=0;
        ksav=ksav+1;
        savtim(ksav)=tyme;
        savtb(:,ksav)=tb(:,2);
        savtp(:,ksav)=tp(:,2);
        savt(:,ksav)=t(:);
        savtfout(ksav)=t(np);
        savgt(ksav)=SolarEnergy;
        saveff(ksav)=eff;
        savamd(ksav)=FlowRate;
        savqu(ksav)=qu;
    end
    % reset temperatures
    for i=1:np
        tp(i,1)=tp(i,2);
        tb(i,1)=tb(i,2);
    end
end  % of tyme loop

%% Plot results
figure(1)
subplot(3,1,1); plot(savtim/3600,savtfout,'k-')
xlabel('Time (hour)')
ylabel('Temperature (C)')
titl=['Model Response, dt = ',num2str(dt),' (sec)',...
    ', dx=',num2str(dx,3),'(m)'];
title(titl)
grid
subplot(3,1,2); plot(savtim/3600,savgt,'k-')
 xlabel('Time (hour)')
ylabel('Gt (W/m2)')
grid
subplot(3,1,3); plot(savtim/3600,savamd,'k-')
xlabel('Time (hour)')
ylabel('Flow Rate (kg/s)')
axis([0 tstop/3600 0 0.1])
grid
figure(2)
subplot(1,2,1); plot(x,savtp(:,1),'k-',x,savtp(:,fix(end/8)),'k:'...
    ,x,savtp(:,fix(end/2)),'k--',x,savtp(:,end),'k-.')
xlabel('Distance from Air Inlet (m)')
ylabel('Temperature (C)')
titl=['Tp for dt=',num2str(dt),...
    '(sec)',', dx=',num2str(dx),'(m)'];
title(titl)
t1=['t=',num2str(savtim(1)/3600,3),' (hour)'];
t2=['t=',num2str(savtim(fix(end/8))/3600,3),' (hour)'];
t3=['t=',num2str(savtim(fix(end/2))/3600,3),' (hour)'];
t4=['t=',num2str(savtim(fix(3*end/4))/3600,3),' (hour)'];
t5=['t=',num2str(savtim(end)/3600,3),' (hour)'];
legend(t1,t2,t3,t5)
grid
subplot(1,2,2); plot(x,savt(:,1),'k-',x,savt(:,fix(end/8)),'k:'...
    ,x,savt(:,fix(end/2)),'k--',x,savt(:,end),'k-.')
xlabel('Distance from Air Inlet (m)')
% ylabel('Temperature (C)')
titl=['Tf for dt=',num2str(dt),...
    '(sec)',', dx=',num2str(dx),'(m)'];
title(titl)
legend(t1,t2,t3,t5)
grid
% 3D surface plots
% pick out selected data for surface plotting
% the mesh graph gets too cluttered if you don't
% do this step
xs=x(1:4:length(x));
ts=savtim(1:10:length(savtim));
temp_air=savt(1:4:length(x),1:10:length(savtim));
temp_abs=savtp(1:4:length(x),1:10:length(savtim));
[X,T]=meshgrid(xs,ts);
figure(3);
mesh(X,T/3600,temp_air');
title('Collector Air Temperature');
xlabel('Distance (m)')
ylabel('Time (hour)')

zlabel('Temperature (C)')
figure(4);
mesh(X,T/3060,temp_abs');
title('Collector Absorber Temperature');
xlabel('Distance (m)')
ylabel('Time (hour)')
zlabel('Temperature (C)')
