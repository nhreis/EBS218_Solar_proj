Variables
% Matlab program to solve transient flat-plate air heater model
% Solves equations from ebs 218 lecture notes using explicit finite
% difference method.

% variable definitions
% d = air duct height (m)
% w = air duct width (m)
% al = collector length (m)
% ep = absorber emittance
% eb = back plate emittance
% ta = ambient air temperature (c)
% gt = solar radiation incident on collector (w/m2)
% taualfa = collector transmittance-absorbtance product
% amdot = mass flow rate of air (kg/s)
% tfin = inlet air temperature (c)
% tilt = collector tilt (deg)
% vw = wind speed (m/s)
% anc = number of acrylic covers
% eg = emittance of acrylic
% ub = back heat loss coefficient (w/m2-c)
% pcset = maximum percent change in iteration on temperatures (%)

clc
%clear all

% Establish specifications for single glazed air-heating solar energy
% collector
sigma = 5.67e-8;    % Stefan-Boltzmann constant (W/m^2*K^4)
d = 0.0508;    % air duct height (m)
w = 0.33;  % collector width (m)
al = 3; % collector length (m)
ep = 0.95;  % emissivity of underside of absorber plate
eb = 0.001;  % emissivity of rear plate
ta = 20; % ambient air temperature in collector (degrees C)
taualfa = 0.85; % collector transmittance-absorbtance product
tfin = 30;  % air inlet temperature (degrees C)
beta = 30;
vw = 25;    % ambient wind speed
anc = 1;    % number of acrylic covers
eg = 0.94;  % emissivity of acrylic cover
ub = 0; % back heat loss coefficient (W/m^2-C)

pos = [38.55, 121.74];
SurAz = 0;

% Absorber thermal properties for Alumin Alloy 195, Table A-3,
% pg 948, Heat Transfer: A Practical Approach by Y.A. Cengel
% McGraw-Hill Book Co., New York.
rhop=2700;  % density of aluminum absorber (kg/m^3 = g/L)
cpp=883;    % specific heat of aluminum plate (J/kg/K)
delp=0.001; % plate thickness (m)
akp=235;    % thermal conductivity of aluminum plate (W/m-K)
rhob=30;  % density of pink insulation (kg/m^3 = g/L)
cpb=5000;    % specific heat of pink insulation (J/kg/K)
delb=0.25; % insulation thickness (m)
akb=0.05;    % thermal conductivity of pink insulation (W/m-K)

fprintf('\nTransient Flat Plate Air Heater Model \n')
fprintf('d (m) = %5.3f\n',d)
fprintf('w (m) = %5.1f\n',w)
fprintf('al (m) = %5.1f\n',al)
fprintf('beta (deg) = %5.1f\n',beta)
fprintf('N = %5.1f\n',anc)
fprintf('tau_alfa (m) = %5.1f\n',taualfa)
fprintf('ep = %5.3f\n',ep)
fprintf('eb = %5.3f\n',eb)
fprintf('eg = %5.3f\n',eg)
fprintf('Ta (C) = %5.1f\n',ta)
fprintf('Vw (m/s) = %5.1f\n',vw)
fprintf('Tfin (C) = %5.1f\n',tfin)

% information for step changes in Gt and air flow rate
gtz=850; % initial value for solar radiation
gtstep=-400; % value for step change in solar radiation
tgtstep=10*3600; % time that step change occurs
amdotz=0.061; % initial value for air flow rate
amdstep=-0.031; % value for step change in air flow rate
tamdstep=0.5*3600; % time that step change occurs
fprintf('Gt(0) (W/m2) = %5.1f\n',gtz)
fprintf('Gt step (W/m2) = %5.1f\n',gtstep)
fprintf('time for Gt step (min) = %5.1f\n',tgtstep/60)
fprintf('mdot(0) (kg/s) = %6.4f\n',amdotz)
fprintf('mdot step (kg/s) = %6.4f\n',amdstep)
fprintf('time for mdot step (min) = %5.1f\n',tamdstep/60)

% parameters for finite difference solution
dt=1;  % time step (sec)
nseg=60;
tstop=3600;
ksav=1;
isav=0;
isv=1;
tyme=0;

% calculate absolute temperatures
tfink=tfin+273.15;
tak=ta+273.15;

% set inital values for air flow and solar radiation
amdot=amdotz;
gt=gtz;

% get air property values from curve fits of data tables
[denf,therm_diff,amu,kin_vis,akf,cpf,prf]=fluid_props_air(tfink);

% collector areas and hydraulic diameter
ac=w*al;
ad=d*w;
dh=2*d;

% variables needed for finite difference calculations
np=nseg+1;
npm1=np-1;
dx=al/nseg;
x=linspace(0,al,np);
c1=dt/(rhop*cpp*delp);
c2=akp*delp/(dx*dx);
c3=dt/(rhob*cpb*delb);
c4=akb*delb/(dx*dx);
c5=amdot*cpf/(w*dx);

% initialize temperature arrays
for i=1:np
    tp(i,1)=ta+.01;
    tp(i,2)=ta+.01;
    tb(i,1)=ta;
    tb(i,2)=ta;
    t(i)=ta+.001;
end

% forced convection due to wind eqn 3.15.3, page 174 d&b
hw=2.8 + 3*vw;

% parameters for top loss coefficient equation (6.4.9), page 260 d&b
if(beta>=70)
    beta=70;
end
c=520*(1-0.000051*beta^2);
f=(1+0.089*hw-0.1166*hw*ep)*(1+0.07866*anc);
term1=1/(ep+0.00591*anc*hw);
term2=(2*anc+f-1+0.133*ep)/eg;
t4=term1+term2-anc;

% set initial values
tpave=(tp(1,2)+tp(np,2))/2;
tbave=(tb(1,2)+tb(np,2))/2;
tfave=(t(1)+t(np))/2;
savtim(ksav)=tyme;
savtp(:,ksav)=tp(:,2);
savt(:,ksav)=t(:);
savt(ksav)=t(np);
savgt(ksav)=gt;
saveff(ksav)=0;
savamd(ksav)=amdot;
ih=0;
im=0;

% time loop
for tyme=dt:dt:tstop
    isav=isav+1;
    tyme
    Clocktime = 6 + tyme/3600;
    [Amb.Gb, Amb.Gd, Amb.Tair, Gt] = irradianceFunc(Clocktime);
    % Where is the sun
    [ SolZen , SolAz, SolAlt, SolDec, solT, hrAng, indAng ] = WhereIsSun( pos, [147, Clocktime * 60], Collect.tilt, SurAz );
    % Calculate solar energy
    s = 1000*RadiationGt(Amb,Ap,Ply,Collect,Cp,indAng,Clocktime)
    %s=gt*taualfa;
    % calculate hrpb, radiation heat transfer coefficient
    tpk=tpave+273.15;
    tbk=tbave+273.15;
    %hrpb=sigma*(tpk+tbk)*(tpk*tpk+tbk*tbk)/(1/ep+1/eb-1)
    hrpb=0;
    % get air property values from curve fits of data tables
    tfk=tfave+273.15;
    [denf,therm_diff,amu,kin_vis,akf,cpf,prf]=fluid_props_air(tfk);
    % reynolds number in duct (see example 3.14.2, p173 d&b)
    re=2*amdot/(w*amu);
    if(re<2100)
        % laminar flow convective heat transfer coefficient (eqn 3.14.7, p172, d&b)
        fac=re*prf*dh/al;
        anud=4.9 + .0606*fac^1.2/(1+0.0909*(fac^0.7)*(prf^0.17));
    else
        % turbulent convective heat transfer coefficient (eqn 3.14.6, p171, d&b)
        anud=0.0158*(re^0.8);
    end
    hcpf=anud*akf/dh;
    hcbf=hcpf
    
    % collector top heat loss coefficient, ut, use empirical correlation
    % eqn (6.4.9) on page 260 of d&b.
    if tpk > tak
        et=0.430*(1-100/tpk);
        t1=anc/((c/tpk)*((tpk-tak)/(anc+f))^et);
        t2=1/(t1+1/hw);
        t3=sigma*(tpk+tak)*(tpk*tpk+tak*tak);
        utop=t2+t3/t4;
        ul=utop+ub
    else
        utop=0;
    end
    
    %Fin
    FinThickness = delp;
    Ac = FinThickness * dx;
    P = 2*dx + 2*FinThickness;
    m=(hcpf * P / ( akp * Ac) )^(1/2);
    qcoef=(hcpf * P * akp * Ac )^(1/2) * tand( m * d );
    
    % finite difference equations to calculate temperatures
    % distance loop
    t(1)=tfin;
    for i=2:npm1
        if (dx*i <=0.67) || (dx*i>=1 && dx*i<=1.33)||(dx*i>=2.33) % One side fin
            N_fin=1;
        elseif (dx*i>=1.33&&dx*i<=1.67) % Two side fin
            N_fin=2;
        else % No side fin
            N_fin=0;
        end
        tp(i,2)=tp(i,1) + c1*(s+c2*(tp(i+1,1)-2.*tp(i,1)...
            +tp(i-1,1))-hcpf*(tp(i,1)-t(i))-hrpb*(tp(i,1)...
            -tb(i,1))-utop*(tp(i,1)-ta)...
            - N_fin * qcoef * ( tp(i,1) - t(i) ));
        tb(i,2)=tb(i,1) + c3*(c4*(tb(i+1,1)-2.*tb(i,1)...
            +tb(i-1,1))-hcbf*(tb(i,1)-t(i))+hrpb*(tp(i,1)...
            -tb(i,1))-ub*(tb(i,1)-ta));
        t(i)=(hcpf*tp(i,2)+hcbf*tb(i,2)+c5*t(i-1)...
            + N_fin * qcoef * ( tp(i,1) - t(i) ))/(hcpf+hcbf+c5);
    end % of distance loop
    
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
    qu=amdot*cpf*(t(np)-t(1));
    eff=100.*qu/(gt*ac);
    
    % save some data for plotting
    if(isav>=isv)
        isav=0;
        ksav=ksav+1;
        savtim(ksav)=tyme;
        savtb(:,ksav)=tb(:,2);
        savtp(:,ksav)=tp(:,2);
        savt(:,ksav)=t(:);
        savtfout(ksav)=t(np);
        savgt(ksav)=gt;
        saveff(ksav)=eff;
        savamd(ksav)=amdot;
        savqu(ksav)=qu;
    end
    
    % reset temperatures
    for i=1:np
        tp(i,1)=tp(i,2);
        tb(i,1)=tb(i,2);
    end
end % of tyme loop
figure(1)
subplot(3,1,1); plot(savtim/60,savtfout,'k-')
% xlabel('Time (min)')
ylabel('Temperature (C)')
titl=['Model Response, dt = ',num2str(dt),' (sec)',...
    ', dx=',num2str(dx,3),'(m)'];
title(titl)
grid
subplot(3,1,2); plot(savtim/60,savgt,'k-')
% xlabel('Time (min)')
ylabel('Gt (W/m2)')
grid
subplot(3,1,3); plot(savtim/60,savamd,'k-')
xlabel('Time (min)')
ylabel('Flow Rate (kg/s)')
axis([0 tstop/60 0 0.1])
grid
figure(2)
subplot(1,2,1); plot(x,savtp(:,1),'k-',x,savtp(:,fix(end/8)),'k:'...
    ,x,savtp(:,fix(end/2)),'k--',x,savtp(:,end),'k-.')
xlabel('Distance from Air Inlet (m)')
ylabel('Temperature (C)')
titl=['Tp for dt=',num2str(dt),...
    '(sec)',', dx=',num2str(dx),'(m)'];
title(titl)
t1=['t=',num2str(savtim(1)/60,3),' (min)'];
t2=['t=',num2str(savtim(fix(end/8))/60,3),' (min)'];
t3=['t=',num2str(savtim(fix(end/2))/60,3),' (min)'];
t4=['t=',num2str(savtim(fix(3*end/4))/60,3),' (min)'];
t5=['t=',num2str(savtim(end)/60,3),' (min)'];
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
mesh(X,T/60,temp_air','EdgeColor','black');
title('Collector Air Temperature');
xlabel('Distance (m)')
ylabel('Time (min)')
zlabel('Temperature (C)')
figure(4);
mesh(X,T/60,temp_abs','EdgeColor','black');
title('Collector Absorber Temperature');
xlabel('Distance (m)')
ylabel('Time (min)')
zlabel('Temperature (C)')
