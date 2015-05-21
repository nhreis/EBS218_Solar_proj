%% Solar Collector Model, EBS 218 Spring Quarter, 2015

%% Define Objects (Structs) to be passed to functions

%Collector variables
Collect.tilt = 0;   %Tile angle of the collector, degrees
Collect.Az   =0;    %Azimuth angle of collector, degrees
Collect.lat = 38.5;  %Latitude position of collector, degrees (Davis)
Collect.lng = 121.7; %Collector longitude, degrees (Davis)
Collect.incident = 15; %Solar Incidence angle on collector, degrees
Collect.mDot = 0.04;  %Mass flow rate of air through collector, kg/s

%Aluminum Absorber Plate
Ap.T = 100;  %Temperature, Degrees C
Ap.A = 1;    %Area, m^2
Ap.thick = .00124; %Thickness of plate, m (18 gauge)
Ap.alpha = 0.98;    %Absorbptivity of black paint (0< alpha < 1) [Incropera & DeWitt 6th ed, pg 956]
Ap.epsilon = 0.98;  %Emissivity of black paint, [Incropera & DeWitt 6th ed, pg 956]
Ap.C = 904;        %Specific heat, J/kgK [WolframAlpha]
Ap.rho = 2700;    %Density Aluminum, kg/m^3 [WolframAlpha]
Ap.k = 235;     %Thermal Conductivity of Al, W/(m K)  [WolframAlpha]
Ap.W = 1;       %Width of collector, m

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
Cp.n = 1.49;    %Index of refraction of Acrylic [http://en.wikipedia.org/wiki/List_of_refractive_indices]
Cp.kextinct =0; %Coefficient of extinction   [http://www.filmetrics.com/refractive-index-database/Acrylic/Acrylate-Lucite-Plexiglass]

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

%Insulating Pink Foam
Foam.T= 20;  %Foam temperature, C
Foam.A = 1; %Foam area, m^2
Foam.R=1.37;  %R-value of the foam, K m^2/W http://insulation.owenscorning.com/assets/0/428/429/431/fd144b72-c513-45e0-99d8-9ef16c90a65b.pdf
Foam.thick = 0.0254; %thickness of foam, m  (1 in? is this right?)
Foam.rho = 0;  %density of the foam ,kg/m^3 ======need value=====


%Ambient conditions
%Includes properties of air
Amb.Tair = 30; %Ambient air temperature, C
Amb.Tsky = 0;  %Sky Temperature ======need value/function for this!=========
Amb.Gb = 800;  %Beam radiation on Pyronometer (parallel with collector), W/m^2
Amb.Gd = 5;    %Diffuse radiation on Pyronometer (parallel with collector), W/m^2
Amb.Gr = 1;    %Reflected solar radiation on Collector, W/m^2
Amb.wind = 3;  %Wind speed, m/s 
Amb.windDir = 0; %Wind direction, compass degrees
Amb.RH = 0.60; %Relative Humidity, %
Amb.day = 140; %Day on which test was run, n^th day of the year
Amb.t = 12*60;     %Time at which data point was taken, min since midnight
%These will be inputs from the model we did for HW, which should be
%converted to a function.
Amb.dec = 20;  %Solar declanation, degrees
Amb.hrAng = 15; %Hour Angle, degrees
Amb.Az = 0; %Solar Azimuth angle, degrees
Amb.zen = 45;  %Solar zenith angle
Amb.CpAir = 1007; %Specific Heat constant pressure of air, J/(kg K) [WolframAlpha]
Amb.CvAir = 721; %Specific heat constant volume air, J/(kg K), [http://www.ohio.edu/mechanical/thermo/property_tables/air/air_Cp_Cv.html]
Amb.kAir = 0.02587;  %thermal conductivity of air, W/(m K) [WolframAlpha]
Amb.rhoAir = 1.164; %Density air at 30C, kg/m^3 [WolframAlpha]
Amb.vAir = 1.87E-5; %Dynamic viscosity of air at 30C, Pa s, [WolframAlpha]
Amb.SolT = 0;       %Solar time
Amb.SolAlt =0;      %Solar Altitude Angle, degrees


%% Set values to locate the sun based on time of day, tilt of collector, direction of collector.
%Also finds solar incident angle on collector based on tilt angle and time
%of day
[Amb.Zen , Amb.Az, Amb.SolAlt, Amb.dec, Amb.solT, Amb.hrAng, Collect.incident] = WhereIsSun([Collect.lat, Collect.lng], [Amb.day, Amb.t], Collect.tilt, Collect.Az);

