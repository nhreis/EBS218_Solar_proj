function Solar = Transimission(Gb,Gd,tilt,Theta,Zen)
    %Calculates and returns radiative heat transfers for component of struct/objects(?) in Watts.
    
%% Incoming Solar Radiation
sigma = 5.6697e-8;  %Boltzmann constant as defined in Textbook
%Transmittance of Acrylic cover
if Theta>90
    Theta=90;
end
if Zen>85
    Zen=85;
end
Rb  = cosd(Theta)/cosd(Zen);    %Geometric ratio ======GUESS FROM NOTES======
rhog = 0.2;  %Diffuse reflectance off ground
n1 = 1;     %refractive index of air
n2 = 1.49;  %refractive index of acrylic cover
%Theta %Incidence angle
kextinct=0;
thick=0.00635;
Theta_effGround = (90 - 0.5788*tilt + 0.002693*(tilt)^2); %Effective incidence angle for ground reflected radiation
Theta_effDiffuse = (59.7 - 0.1388*tilt + 0.001497*tilt^2); %Effeective incidence angle for diffuse  radiation

ThetaVec = {Theta Theta_effDiffuse Theta_effGround};   %Defining the order of Theta in vector                                 %Initialise vector space
I = 1;                            
for i = 1:length(ThetaVec)
    Theta2 = asind((n1/n2)*sind(ThetaVec{i}));          %Refracted angle
    alpha = 0.93.*(1 + 1.364e-3.*ThetaVec{i} - 2.745e-4.*ThetaVec{i}.^2 + 1.767e-5.*ThetaVec{i}.^3 ... 
                    - 5.282e-7.*ThetaVec{i}.^4 + 7.117e-9.*ThetaVec{i}.^5 - 3.604e-11.*ThetaVec{i}.^6); %Abosorptance
    r_perp = (sind(Theta2-ThetaVec{i}).^2)./(sind(Theta2+ThetaVec{i}).^2); %Perpendicular unpolarised reflected radiation 
    r_par = (tand(Theta2-ThetaVec{i}).^2)./(tand(Theta2+ThetaVec{i}).^2);  %parallel unpolarised reflected radiation 
    tau_r = 0.5.*((I-r_par)./(I+r_par)+(I-r_perp)./(I+r_perp));            %Total Transmittance   
    tau_a = exp(-kextinct*thick./(cosd(Theta2)));                  %Accounting for absoprtion losses in acrylic cover
    tau = tau_a .* tau_r;            
    TauAlpha{i} = 1.01.*tau.*alpha;
end
G = Gb +Gd;
Solar = Gb*Rb.*TauAlpha{1} + Gd.*TauAlpha{2}.*((1+cosd(tilt))/2); %+ G*rhog.*TauAlpha{3}.*((1-cosd(tilt)/2));
