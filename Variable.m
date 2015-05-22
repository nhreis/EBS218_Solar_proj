clc
clear all
close all

%% Solar Collector Model, EBS 218 Spring Quarter, 2015
% Initialise Sun-Earth Geometry

n = 198;           %nth Day of year
Lat = 38.5;        %Latitude in degrees North=positive
Long = 121.5;      %Longitude West of Prime meridian 
L_st = 120;        %Standard Meridian West of Prime meridian
t = 6:18;         %Standard Time (can be scalar or vector)

[ThetaZ,GammaS,SolarTime] = SunLocation(n,Lat,Long,L_st,t);
beta = 90;      %Surface tilt angle in degrees (0 faceup to 180 facedown)
Gamma = -45;     %Surface azimuth angle relative to local meridian (Zero due south, East negative) 
Theta = zeros(1,length(t)); %initialise null space vector for theta
for i = 1:length(t);
    Theta(i) = acosd(cosd(ThetaZ(i))*cosd(beta) + sind(ThetaZ(i))*sind(beta)*cosd((GammaS(i)-Gamma)));
end

for i = 1:length(Theta)
    if Theta(i) > 90
       Theta(i) = 90;
    end
end


%% Incoming Measured Data

Gb = 515;   %Incoming beam Radiation
Gd = 98;    %Incoming diffuse Radiation


%% Define Objects (Structs) to be passed to functions

%Collector variables
Colect.tilt = 30;   %Tile angle of the collector, degrees
Colect.Az   = 0;    %Azimuth angle of collector, degrees


%Aluminum Absorber Plate
Ap.T = 100;  %Temperature, Degrees C
Ap.A = 1;    %Area, m^2
Ap.thick = .00124; %Thickness of plate, m (18 gauge)
Ap.alpha = 0.98;    %Absorbptivity of black paint (0< alpha < 1) [Incropera & DeWitt 6th ed, pg 956]
Ap.epsilon = 0.98;  %Emissivity of black paint, [Incropera & DeWitt 6th ed, pg 956]
Ap.C = 904;        %Specific heat, J/kgK [WolframAlpha]
Ap.rho = 2700;    %Density Aluminum, kg/m^3 [WolframAlpha]

%Acrylic Cover Plate
Cp.T = 50;    %Temperature, C
Cp.A = 1;     %Area, m^2
Cp.thick = 0.00635;  %1/4" thick plate, m.
Cp.alpha = 0.05;    %Absorbptance  ===THIS IS A GUESS===
Cp.relect = 0.05;   %Reflectivity  ===ALSO A GUESS===  
Cp.tau = 0.894;     %Lecture 9, Slide 4
Cp.C = 1470;       %Specific Heat, J/kgK
Cp.k = 0.2;       %Thermal Conductivity, W/(m K)[WolframAlpha]
Cp.rho = 1180;   %Density, kg/m^3 [WolframAlpha]
Cp.epsilon = 0;  %Emissivity of acrylic =====NEED VALUE=====
Cp.n = 1.25;    %index of refraction of acrylic =========THIS IS A GUESS=======
Cp.kextinct = 0.1; %coefficient of extinction ========GUESSED VALUE=======.

%2x4 Walls
Wall.T= 40;              %Wall Tempearture, C
Wall.thick = 0.0381;  %Wall thickmess, m
Wall.height = 0.0889; %Wall height, m
Wall.L = 1;           %Wall length, m
Wall.A = 4*Wall.height*Wall.L;  %Total surface area of walls, m^2
Wall.epsilon = 0.87;    %Emissivity of wood, [Incropera & DeWitt 6th ed, pg 955]
Wall.rho = 500;      %Estimation of wood density, kg/m^3 [Incropera & DeWitt 6th ed, pg 940]
Wall.C = 2806;      %Specific heat of Yellow Pine, J/kgK [Incropera & DeWitt 6th ed, pg 940]
Wall.k = 0.13;      %Estimation of thermal conductivity of wood, W/(m K), [Incropera & DeWitt 6th ed, pg 940]
Wall.alpha = 0;     %Absorbance of wood ======NEED VALUE!=========
Wall.R = 4.375*0.176; %R-Value of plywood, Km^2/W http://coloradoenergy.org/procorner/stuff/r-values.htm

%Plywood Base
Ply.T = 20;   %Plywood base temp, 20C  *assumption*
Ply.A = 1;    %Area of plywood base, m^2
Ply.thick = 0.0127; %Thickness of plywood, m (.5in?) ===is this correct?===
Ply.k = 0.125; %Thermal conductibity of plywood, W/(m K)
Ply.R = 0.31*0.176;  %R-Value of Plywood, K m^2/W,  http://coloradoenergy.org/procorner/stuff/r-values.htm
Ply.C = 2500;  % Specific heat plywood, J/kg K   http://www.ehow.co.uk/info_8534542_thermal-properties-plywood.html
Ply.epsilon = 0.82; %http://www.infrared-thermography.com/material-1.htm

%Insulating Pink Foam
Foam.T= 20;  %Foam temperature, C
Foam.A = 1; %Foam area, m^2
Foam.R=1.37;  %R-value of the foam, K m^2/W http://insulation.owenscorning.com/assets/0/428/429/431/fd144b72-c513-45e0-99d8-9ef16c90a65b.pdf
Foam.thick = 0.0254; %thickness of foam, m  (1 in? is this right?)
Foam.rho = 0;  %density of the foam ,kg/m^3 ======need value=====


%Ambient conditions
Amb.Tair = 30; %Ambient air temperature, C
Amb.Tsky = 0;  %Sky Temperature ======need value/function for this!=========
Amb.Gb = 800;  %Beam radiation on horizontal, W/m^2
Amb.Gd = 5;    %Diffuse radiation on horizontal, W/m^2
Amb.wind = 3;  %Wind speed, m/s 
Amb.windDir = 0; %Wind direction, compass degrees
Amb.RH = 0.60; %Relative Humidity, %
Amb.t = 12;     %Time at which data point was taken, h (24)


%% Radiation Function Test

[Solar,Q_ApSky,Q_ApPly,Q_ApGrnd] = Radiation(Amb,Ap,Ply,Colect,Cp,Theta,t);










