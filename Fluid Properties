function [den,therm_diff,dyn_vis,kin_vis,therm_con,cp,Pr]=fluid_props_air(Tk)
% Matlab function to evaluate properties of air at one atmosphere
% Data taken from Table A-11, page 961 of Cengel, Y.A.  1998.
% Heat Transfer: A Practical Approach.  WCB McGraw-Hill, Boston.
% Temerpature range of data extracted from Table is from 200 to 400 K
% Curve fits done in Matlab script: air_props_fit.m, using Matlab Polyfit.
% Tk = temperature (K)
% den = density (kg/m3)
% cp = specific heat (J/kg C)
% cond = thermal conductivity (W/m C)
% diff = thermal diffusivity (m2/s)
% mu = dynamic viscosity (kg/m s)
% nu = kinematic viscosity (m2/s)
% Pr = Prandtl number
den=1.4162e-5*Tk*Tk-1.2798e-2*Tk+3.7474;
therm_diff=1.4659e-10*Tk*Tk+4.6488e-8*Tk-5.0077e-6;
dyn_vis=-3.3899e-11*Tk*Tk+6.7723e-8*Tk+1.2295e-6;
kin_vis=1.0568e-10*Tk*Tk+2.8644e-8*Tk-2.3601e-6;
therm_con=-5.1906e-8*Tk*Tk+1.0627e-4*Tk-1.1287e-3;
cp=2.8678e-4*Tk*Tk-1.2256e-1*Tk+1.0160e3;
Pr=9.0437e-7*Tk*Tk-7.2540e-4*Tk+8.4878e-1;
end
