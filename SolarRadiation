function Gt= SolarRadiation(time,Gb, Gd)
%Total solar radiation on the tilted surface using Liu-Jordan Model
%Gb and Gd are data collected by pyranometer

rho_g=0.2; %Pre-set diffuse reflectance of ground
beta=15; %Pre-set tilted angle
pos=[38.5539, 121.7381]; % Pre-set position
SurAz=0; % Pre-set Surface Azimuth angle

%Calculate position of the sun
[ SolZen , SolAz, SolAlt, SolDec, solT, hrAng, indAng] = WhereIsSun( pos, time, beta, SurAz )
 
Gt=cosd(indAng)/cosd(SolZen)*Gb + Gd * ((1+cosd(beta))/2) + ...
    (Gb+Gd)* rho *((1-cosd(beta))/2);

end
