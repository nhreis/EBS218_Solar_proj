% EBS 218 S 2015, Professor Jenkins
% Matlab script to calculate solar radiation on a tilted surface
% using the Liu and Jordan model. 

% data is in an excel file called Solar2014csv.csv
% data is organized as such: Year, Julian Day of Year, Hour of Day (PST),
% Air Temperature (°C), PV Module Temperature (°C), 
% Pyranometer—South at Latitude (W/m2), Pyranometer—Horizontal (W/m2)
% Pyranometer—South at Latitude minus 16°  (W/m2), 
% Pyranometer—South at Latitude plus 16°  (W/m2), 
% Pyranometer— South Vertical  (W/m2)

clc 
clear all

% constants
rho_g = 0.20;   % diffuse reflectance of ground
tilt = 15;    % tilt angle of solar collector (deg)
lat = 38;   % latitude (deg North)
time = [0:23];  % time of day (24 hr clock)
gamma = 0;  % surface azimuth angle (zero due south)

% read variables from file
davis = csvread('Solar2014csv.csv');    % reads data from CSV file
Gt_vert = davis(3456:3479,10);  % south vertical pyranometer (W/m2)
G_horz = davis(3456:3479,7);    % horizontal pyranometer (W/m2)
day = davis(3456,2);    % day of year


% calculate zenith and incidence angle for pyranometers
beta = 90;  % slope angle (vertical surface)
omega = 15*(time-12);   % hour angle (deg)
decl = 23.45*sin(360*(284+day)/365);    % declination angle (deg)
theta_z = acos(cos(lat)*cos(decl)*cos(omega)+sin(lat)*sin(decl));   % zenith angle
theta = acos(sin(decl)*sin(lat)*cos(beta)...
        -sin(decl)*cos(lat)*sin(beta)*cos(gamma)...
        +cos(decl)*cos(lat)*cos(beta)*cos(omega)...
        +cos(decl)*sin(lat)*sin(beta)*cos(gamma)*cos(omega)...
        +cos(decl)*sin(beta)*sin(gamma)*sin(omega));    % incidence angle 

% calculate G_d and G_b 
G_d = (Gt_vert - G_horz*(rho_g*((1-cos(beta))/2)+cos(theta)/cos(theta_z)))...
      /(((1+cos(beta))/2)-cos(theta)/cos(theta_z)); % diffuse radiation on horizontal surface
G_b = G_horz - G_d; % beam radiation on horizontal surface

% calculate G_T
G_T = G_b*(cos(theta)/cos(theta_z))+G_d*((1+cos(180-tilt))/2)...
      +G_horz*rho_g*((1-cos(180-tilt))/2);  % irradiance on solar collector
G_T = cat(2,time',G_T); % irradiance vs. time
 
 

  
        
