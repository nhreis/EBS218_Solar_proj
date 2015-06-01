DuctLength=3*Ap.W;
DuctWidth=1/3*Ap.W
Nseg=30; %Number of segments
x=linspace(0,DuctLength,Nseg+1);
dx=DuctLength/Nseg;
DuctHeight=0.0508;

Bp_k=0.003;
Bp_thick=0.25;

Tstop=3600*8; %Calculate this value by compute the day length (s)
dt=1; %time step (s)

%The coeficients in finite difference equations
%For absorber plate
c1p=dt/(Ap.thick * Ap.C * Ap.rho);
c2p=Ap.k * Ap.thick/(dx * dx);

pos = [38.55, 121.74];
SurAz = 0;
FlowRate = Collect.mDot;
ta=20;
Vw=Amb.wind;
sigma = 5.6697e-8;
Bp_epsilon=0;

%For insulation plate
c1b=dt/(Bp_thick * 10000 * Foam.rho);
c2b=Bp_k * Bp_thick/(dx*dx);



% initialize temperature arrays
for i=1:(Nseg+1)
    tp(i,1)=ta+.01;
    tp(i,2)=ta+.01;
    tb(i,1)=ta;
    tb(i,2)=ta;
    t(i,1)=ta+.001;
    t(i,2)=ta+.001;
end

tpk = ( tp(1,2) + tp(Nseg+1,2) ) / 2 + 273.15;
tbk = ( tb(1,2) + tb(Nseg+1,2) ) / 2 + 273.15;
tk = ( t(1,2) + t(Nseg+1,2) ) / 2 + 273.15;

%The coeficients in finite difference equations
%For the fluid
[denf,therm_diff,amu,kin_vis,akf,cpf,prf]=fluid_props_air(tk);
cair=FlowRate*cpf/(DuctWidth*dx);



%time loop
for time=dt:dt:Tstop
    k=time/dt
    
    Clocktime = 6 + time/3600;
    [Amb.Gb, Amb.Gd, Amb.Tair, Gt] = irradianceFunc(Clocktime);
    % Where is the sun
    [ SolZen , SolAz, SolAlt, SolDec, solT, hrAng, indAng ] = WhereIsSun( pos, [147, Clocktime * 60], Collect.tilt, SurAz );
    % Calculate solar energy
    SolarEnergy = 1000*RadiationGt(Amb,Ap,Ply,Collect,Cp,indAng,Clocktime);
    % Calculate top loss coef
    Utop=TopLoss(tpk, Amb.Tair, Vw, Collect.tilt , Ap.epsilon , Cp.epsilon);
    
    % Calculate radiation transfer coef
    hrpb=sigma*(tpk+tbk)*(tpk*tpk+tbk*tbk)/(1/Ap.epsilon+1/Bp_epsilon-1)
    
    %Get air properties
    [den,therm_diff,dyn_vis,kin_vis,therm_con,AirC,Pr]=fluid_props_air(tk);
    
    % Convection coef
    [hcpf,hcbf] = ConvectionCoef(FlowRate, tk, dyn_vis, Pr, therm_con);
    
    %The coeficients in finite difference equations
    %for air flow
    cair=FlowRate * AirC / ( DuctWidth * dx );
    
    %Fin
    FinThickness = Ap.thick;
    Ac = FinThickness * dx;
    P = 2*dx + 2*FinThickness;
    m=(hcpf * P / ( Ap.k * Ac) )^(1/2);
    qcoef=(hcpf * P * Ap.k * Ac )^(1/2) * tand( m * DuctHeight );
    
    
    % distance loop
    %t(1,k)=ta;
    
    for i=2:Nseg
        
        if (dx*i <=0.67) || (dx*i>=1 && dx*i<=1.33)||(dx*i>=2.33) % One side fin
        N_fin=1;
        elseif (dx*i>=1.33&&dx*i<=1.67) % Two side fin
          N_fin=2;
        else % No side fin
           N_fin=0;
        end
        
        %%Update rear plate temperature
        tb(i,k+1) = t(i,k);
        
        % Update plate temperature
        tp(i,k+1) = tp(i,k) + c1p * ( SolarEnergy + c2p * ( tp(i+1,k) - 2 .* tp(i,k)...
            + tp(i-1,k) ) - hcpf * ( tp(i,k) - t(i,k) ) - hrpb * ( tp(i,k)...
            - tb(i,k) ) - Utop * ( tp(i,k) - ta ) ...
            - N_fin * qcoef * ( tp(i,k) - t(i,k) ) ) ;
        
        % Update rear plate temperature
        %tb(i,k+1)= tb(i,k) + c1b * ( c2b * ( tb(i+1,k) - 2 .* tb(i,k)...
         %  + tb(i-1,k) ) - hcpf * ( tb(i,k) - t(i,k) ) + hrpb * ( tb(i,k)...
          %- tb(i,k) ));
        
        % Update air flow temperature
        t(i,k+1) = ( hcpf * tp(i,k+1) + hcbf * tb(i,k+1) + cair * t(i-1, k + 1 )...
            + N_fin * qcoef * ( tp(i,k) - t(i,k) ) ) /(hcpf+hcbf+cair);
        %
        
        
        
    end %of distance loop
    
    % boundary conditions
    tp(1,k+1)=tp(2,k+1);
    tp(Nseg+1,k+1)=tp(Nseg,k+1);
    tb(1,k+1)=tb(2,k+1);
    tb(Nseg+1,k+1)=tb(Nseg,k+1);
    t(1,k+2) = t(1,k+1);
    t(Nseg+1,k+1)=(hcpf*tp(Nseg+1,k+1)+hcbf*tb(Nseg+1,k+1)+cair*t(Nseg,k+1) ...
        + qcoef * ( tp(Nseg+1,k) - t(Nseg+1,k) ))/(hcpf+hcbf+cair);
    
    % Kelvin temperature
    tpk = ( tp(1,k+1) + tp(Nseg+1,k+1) ) / 2 + 273.15;
    tbk = ( tb(1,k) + tb(Nseg+1,k) ) / 2 + 273.15;
    tk = ( t(1,k) + t(Nseg+1,k) ) / 2 + 273.15
    
    % useful heat and efficiency
    qu(k) = FlowRate * AirC * ( t(Nseg+1,k)-t(1,k));
    eff(k)=100.*qu(k)/(Gt*Ap.A);
    
    
end % of time loop
