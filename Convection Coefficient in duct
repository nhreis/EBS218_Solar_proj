function hc = ConvectionCoef(FlowRate, MeanT_air)
w=0.33 % Pre-set collector width 0.33(m) for one dimensional model
d=0.012 %air duct height (m)
dh=2*d;
l=3; %Collector length for one dimension (m)
%get air property values from curve fits of data tables
[denf,therm_diff,dyn_vis,kin_vis,akf,cpf,prf]=fluid_props_air(MeanT_air);
% Reynolds number in duct (see example 3.14.2, p173 d&b)
Re=2*FlowRate/(w*dyn_vis);
if(Re<2100)
% laminar flow convective heat transfer coefficient (eqn 3.14.7, p172, d&b)
fac=Re*prf*dh/l;
anud=4.9 + .0606*fac^1.2/(1+0.0909*(fac^0.7)*(prf^0.17));
else
% turbulent convective heat transfer coefficient (eqn 3.14.6, p171, d&b)
anud=0.0158*(Re^0.8);
end
hcpf=anud*akf/dh;
hcbf=hcpf;
end
