DuctLength=3*Ap.W;
DuctWidth=1/3*Ap.W
Nseg=60; %Number of segments
x=linspace(0,DuctLength,Nseg+1);
dx=DuctLength/Nseg;


Tstop=3600*8; %Calculate this value by compute the day length (s)
dt=10; %time step (s)

%The coeficients in finite difference equations
%For absorber plate
c1p=dt/(Ap.thick * Ap.C * Ap.rho);
c2p=Ap.k * Ap.thick/(dx * dx);

%For insulation plate




% initialize temperature arrays
for i=1:(Nseg+1)
    tp(i,1)=ta+.01;
    tb(i,1)=ta;
    t(i,1)=ta+.001;
end
t(1,2)=ta+.001;


%time loop
for time=dt:dt:Tstop
    k=time/dt;
    % Where is the sun
    [ SolZen , SolAz, SolAlt, SolDec, solT, hrAng, indAng ] = WhereIsSun( pos, t, tilt, SurAz )
    % Calculate solar energy
    [SolarEnergy,Q_ApSky,Q_ApPly,Q_ApGrnd] = Radiation(Amb,Ap,Ply,Colect,Cp,indAng,t);
    % Calculate top loss coef
    Utop=TopLoss(tp(1,k), tp(Nseg+1,k), Ta, Vw, Collect.tilt , Ap.epsilon , Cp.epsilon);
    
    % Kelvin temperature
    tpk = ( tp(1,k) + tp(Nseg+1,k) ) / 2 + 273.15;
    tbk = ( tb(1,k) + tb(Nseg+1,k) ) / 2 + 273.15;
    tk = ( t(1,k) + t(Nseg+1,k) ) / 2 + 273.15;
    
    % Calculate radiation transfer coef
    hrpb=sigma*(tpk+tbk)*(tpk*tpk+tbk*tbk)/(1/Ap.epsilon+1/Cp.epsilon-1);
    
    %Get air properties
    [den,therm_diff,dyn_vis,kin_vis,therm_con,Air.C,Pr]=fluid_props_air(tk);
    
    % Convection coef
    [hcpf,hcbf] = ConvectionCoef(Airflow, tk, dyn_vis, Pr, Air.C)
    
    %The coeficients in finite difference equations
    %for air flow
    cair=Airflow * Air.C / ( DuctWidth * dx );
    
    %Fin
    FinThickness = Ap.thick;
    Ac = FinThickness * dx;
    P = 2w + 2t;
    m=(hcpf * P / ( Ap.k * Ac) )^(1/2)
    qcoef=(hcpf * P * Ap.k * Ac )^(1/2) * tand( m * DuctHeight );
    
    
    % distance loop
    t(1,k)=T_in;
    
    for i=2:Nseg
        
        if (dx*i <=0.67) || (dx*i>=1 && dx*i<=1.33)||(dx*i>=2.33) % One side fin
            N_fin=1;
        elseif (dx*i>=1.33&&dx*i<=1.67) % Two side fin
            N_fin=2;  
        else % No side fin
            N_fin=0;
        end
        
        % Update plate temperature
        tp(i,k+1)= tp(i,k) + c1p * ( SolarEnergy + c2p * ( tp(i+1,k) - 2 .* tp(i,k)...
            + tp(i-1,k) ) - hcpf * ( tp(i,k) - t(i,k) ) - hrpb * ( tp(i,k)...
            - tb(i,k) ) - utop * ( tp(i,k) - Ta ) ...
            - N_fin * qcoef * ( tp(i,k) - t(i,k) ) ) ;
        
        % Update rear plate temperature
        tb(i,k+1)= tb(i,k) + c1b * ( c2b * ( tb(i+1,k) - 2 .* tb(i,k)...
            + tb(i-1,k) ) - hcpf * ( tb(i,k) - t(i,k) ) + hrpb * ( tb(i,k)...
            - tb(i,k) ));
        
        % Update air flow temperature
        t(i,k+1) = ( hcpf * tp(i,k+1) + hcbf * tb(i,k+1) + cair * t(i-1, k + 1 )...
            + N_fin * qcoef * ( tp(i,k) - t(i,k) ) ) /(hcpf+hcbf+cair);
        
        
        
    end %of distance loop
    
    % boundary conditions
    tp(1,k+1)=tp(2,k+1);
    tp(Nseg+1,k+1)=tp(Nseg,k+1);
    tb(1,k+1)=tb(2,k+1);
    tb(Nseg+1,k+1)=tb(Nseg,k+1);
    t(1,k+2) = t(1,k+1);
    t(Nseg+1,k+1)=(hcpf*tp(Nseg+1,k+1)+hcbf*tb(Nseg+1,k+1)+cair*t(Nseg,k+1) ...
        + qcoef * ( tp(Nseg+1,k) - t(Nseg+1,k) ))/(hcpf+hcbf+cair);
    
    % useful heat and efficiency
    qu(k) = FlowRate * Air.C * ( t(Nseg+1,k)-t(1,k));
    eff(k)=100.*qu(k)/(Gt*Ap.A);
    
end % of time loop
