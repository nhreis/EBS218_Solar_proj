function Utop=TopLoss(Tp1, Tplast, Ta, Vw, beta , Ap.epsilon , Cp.epsilon)
%Tp1 the temperature of the first segment of plate
%Tp1 the temperature of the last segment of plate
%Ta the temperature of the ambient air
%Vw wind speed
%beta tilted angle
Ncover=1; %Number of covers
sigma=5.67e-8;  % Stefan-Boltzmann constant (W/m^2*K^4)
hw=2.8 + 3*vw; % forced convection due to wind  eqn 3.15.3, page 174 d&b

% parameters for top loss coefficient equation (6.4.9), page 260 d&b
if(beta>=70)
beta=70;
end

c=520*(1-0.000051*beta^2);
f=( 1 + 0.089 * hw - 0.1166 * hw * Ap.epsilon ) * ( 1 + 0.07866 * Ncover);
term1 = 1 / ( Ap.epsilon + 0.00591 * Ncover * hw );
term2 = ( 2 * Ncover + f - 1 + 0.133 * Ap.epsilon ) / Cp.epsilon;
t4 = term1 + term2 - Ncover;

tpk=(Tp1+Tplast)/2 +273.15; % average plate temperature
tak=Ta+273.15; % average air temperature


et = 0.430*(1-100/tpk);
t1 = Ncover / ( ( c / tpk ) * ( ( tpk - tak ) / ( Ncover + f))^et);
t2 = 1/(t1+1/hw);
t3 = sigma * ( tpk + tak ) * ( tpk * tpk + tak * tak );
Utop = t2 + t3 / t4;
end
