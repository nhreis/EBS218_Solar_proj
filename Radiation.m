function [Solar,Q_ApSky,Q_ApPly,Q_ApGrnd] = Radiation(Amb,Ap,Ply,Colect,Cp,Theta,t)
    %Calculates and returns radiative heat transfers for component of struct/objects(?) in Watts.
    
%% Incoming Solar Radiation
sigma = 5.6697e-8;  %Boltzmann constant as defined in Textbook
%Transmittance of Acrylic cover
Rb  = 1.71;    %Geometric ratio ======GUESS FROM NOTES======
rhog = 0.2;  %Diffuse reflectance off ground
n1 = 1;     %refractive index of air
n2 = Cp.n;  %refractive index of acrylic cover
Theta1 = Theta'; %Incidence angle
Theta_effGround = (90 - 0.5788*Colect.tilt + 0.002693*(Colect.tilt)^2)*ones(length(t),1); %Effective incidence angle for ground reflected radiation
Theta_effDiffuse = (59.7 - 0.1388*Colect.tilt + 0.001497*Colect.tilt^2)*ones(length(t),1); %Effeective incidence angle for diffuse  radiation

ThetaVec = {Theta1 Theta_effDiffuse Theta_effGround};   %Defining the order of Theta in vector                                 %Initialise vector space
I = ones(length(t),1);                            
for i = 1:length(ThetaVec)
    Theta2 = asind((n1/n2)*sind(ThetaVec{i}));          %Refracted angle
    alpha = 0.93.*(1 + 1.364e-3.*ThetaVec{i} - 2.745e-4.*ThetaVec{i}.^2 + 1.767e-5.*ThetaVec{i}.^3 ... 
                    - 5.282e-7.*ThetaVec{i}.^4 + 7.117e-9.*ThetaVec{i}.^5 - 3.604e-11.*ThetaVec{i}.^6); %Abosorptance
    r_perp = (sind(Theta2-ThetaVec{i}).^2)./(sind(Theta2+ThetaVec{i}).^2); %Perpendicular unpolarised reflected radiation 
    r_par = (tand(Theta2-ThetaVec{i}).^2)./(tand(Theta2+ThetaVec{i}).^2);  %parallel unpolarised reflected radiation 
    tau_r = 0.5.*((I-r_par)./(I+r_par)+(I-r_perp)./(I+r_perp));            %Total Transmittance   
    tau_a = exp(-Cp.kextinct*Cp.thick./(cosd(Theta2)));                  %Accounting for absoprtion losses in acrylic cover
    tau = tau_a .* tau_r;            
    TauAlpha{i} = 1.01.*tau.*alpha;
end
G = Amb.Gb + Amb.Gd;
Solar = Amb.Gb*Rb.*TauAlpha{1} + Amb.Gd.*TauAlpha{2}.*((1+cosd(Colect.tilt))/2) + G*rhog.*TauAlpha{3}.*((1-cosd(Colect.tilt)/2));



%% Radiative heat transfer: Absorber to Back plate
%qr = hr(T2-T1) is the basic equation
%Defining the radiative heat transfer coefficient hr
hr_Ap2Ply = (sigma*(Ap.T^2 + Ply.T^2)*(Ap.T + Ply.T)) / (1/Ap.epsilon + 1/Ply.epsilon -1);
Q_ApPly = hr_Ap2Ply*(Ap.T-Ply.T);


%% Radiative heat transfer with Sky
T_dp = 5*ones(length(t),1); %Due point temperature taken as the median of the range given in D&B (-20 to 30) Tsky does not affect heating performance much 
T_sky = (Amb.Tair).*(0.711*ones(length(t),1) + 0.0056.*T_dp + 0.000073.*(T_dp).^2 + 0.013.*cos(15.*t')).^(1/4); %Defining temperature of sky as function of time
hr_sky = Ap.epsilon*sigma.*(Ap.T^2*ones(length(T_sky),1) + T_sky.^2).*(Ap.T + T_sky);        %linearized radiative heat transfer coefficient
Q_ApSky = hr_sky*Ap.A.*(Ap.T*ones(length(T_sky),1) - T_sky);



%% Radiative heat transfer with ground

Tgrnd = Amb.Tair;   %assume same as ambient air temperature
hr_Ap2grnd = (sigma*(Ap.T^2 + Tgrnd^2)*(Ap.T + Tgrnd)) / (1/Ap.epsilon + 1/Ply.epsilon -1);
Q_ApGrnd = hr_Ap2grnd*(Ap.T - Tgrnd);
