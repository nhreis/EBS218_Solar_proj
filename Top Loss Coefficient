function Utop=TopLoss(Tp1, Tplast, Ta, Vw)
%Tp1 the temperature of the first segment of plate
%Tp1 the temperature of the last segment of plate
%Ta the temperature of the ambient air
%Vw wind speed
sigma=5.67e-8  % Stefan-Boltzmann constant (W/m^2*K^4)

tpk=(Tp1+Tplast)/2 +273.15;
tak=Ta+273.15;
hw=2.8 + 3*vw; % forced convection due to wind  eqn 3.15.3, page 174 d&b

et=0.430*(1-100/tpk);
t1=anc/((c/tpk)*((tpk-tak)/(anc+f))^et);
t2=1/(t1+1/hw);
t3=sigma*(tpk+tak)*(tpk*tpk+tak*tak);
Utop=t2+t3/t4;
end
