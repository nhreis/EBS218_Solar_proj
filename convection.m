
function hc = convection(P,T,h,Vw,Amb,Ap)
 
    %Inputs 
   %P,Initial air pressure, Pa
   %T,Initial air temperature, K
   %h,Initial Relative Humidity of Air
   %Vw,Wind Speed, m/s

   %Measuring Perimeters and Areas
    Ad % Air Duct Area, m^2
    Pd % Air Duct Perimeter, m
   
    %Determining Saturation Vapor Pressure, Psv
    A= 1.2378847E-05; %Empirical Coefficient, 1/K^2
    B1=-1.9121316E-02; %Empirical Coefficient, 1/K
    C=33.93711047; %Empirical Coefficient
    D=-6.341645E3; %Empirical Coefficient

    Psv=exp(A*T^2+B1*T+C+D/T); %Equation

    %Determining the Enhancement Factor f
    a=1.00662;
    B2=3.14E-8; %1/Pa
    y=5.6E-7; %1/K^2
    t=T-273.15; %Temperature (C)
    f=a+B2*P+y*t^2; %Equation

    %Mole Fraction of Water, Xv
    Xv=h*f*Psv/P; %Equation

    %Defining the Compressibility Factor
    a0= 1.58123E-6; %K/Pa, CIPM-2007
    a1= -2.9331E-8; %1/Pa, CIPM-2007
    a2=1.1043E-10; %1/(K*Pa), CIPM-2007
    b0=0.0005707;  %K/Pa, CIPM-2007
    b1=-2.051E-8;  %1/Pa, CIPM-2007
    c0= 0.00019898; %1/Pa, CIPM-2007
    d=1.83E-11;  %K^2/P^2, CIPM-2007
    e=-7.65E-9; %K^2/P^2, CIPM-2007

    Z= 1-P/T*(a0+a1*t+a2*t^2+(b0+b1*t)*Xv)+c0+c1*t)*Xv^2+P^2/T^2*(d+e*Xv^2); %Equation, CIPM-2007
    
    %Air Density, rho.air
    Ma= 28.96546E-3; %molar mass of dry air, kg/mol, CIPM-2007
    Mv= 0.018; %molar mass of water vapor, kg/mol, CIPM-2007
    R=8.314472; %Universal Gas Constant J/(mol*K), CIPM-2007

    rho_air=P*Ma/(Z*R*T)*(1-Xv*(1-Mv/Ma)); %Equation, CIPM-2007

    %Dynamic and Kinematic Visocities of Air, mu and nu
    b=1.458E-6; %kg/(m*S*K^0.5), http://www-mdp.eng.cam.ac.uk/web/library/enginfo/aerothermal_dvd_only/aero/fprops/propsoffluids/node5.html
    S=110.4; %K, air (http://www-mdp.eng.cam.ac.uk/web/library/enginfo/aerothermal_dvd_only/aero/fprops/propsoffluids/node5.html)

    mu=(b*T^1.5)/(T+S); %http://www-mdp.eng.cam.ac.uk/web/library/enginfo/aerothermal_dvd_only/aero/fprops/propsoffluids/node5.html
    nu=mu/rho_air; %Lecture 5, EBS 218

    %Characteristic Length of Air, L
    L=4*Ad/Pd; %Lecture 5, EBS 218

    %Reynolds Number of Air, Re
    Re=Vw*L/nu; %Lecture 5, EBS 218

    %Volumetric Coefficient of Expansion of Air, Beta
    Beta=1/T; %Shock Absorber Handbook, 2007

    %Thermal Conductivity of Air, kappa
    kappa=0.00075*T+0.0051; %obtained from curve-fitting tabulated data from Engineering ToolBox

    %Thermal Diffusivity of Air, alpha (Dewitt+Incropera, 4th Edition)
    Cp=1.007; %KJ/kg*K, 
    alpha=kappa/(P*Cp);
    
    %Prandtl Number of Air, Pr
    Pr=nu/alpha; %Lecture 5, EBS 218

    %Grashof Number of Air (Lecture 5, EBS 218)
    g=9.81; %gravitational constant, m/s2
    Gr=G*Beta*(Tp-Ta)*L^3/nu^2 ; 

    %Rayleigh Number of the Air (Lecture 5, EBS 218)
    Lp=Ap.thick; %Characteristic Length of the Plate
    Ra=g*Beta*(Tp-Ta)*Ap.thickness^3/(nu*alpha); 

    %Do we have free, forced, or mixed convection?
    if Gr/Re^2 <0.5
       Nusselt= free(Ra);
    elseif Gr/Re^2 <=5
        Nusselt=mixed(forced,free);
    else Gr/Re^2 >5
        Nusselt=forced(Re,Pr);
    end 

    %Free Convection Function
    function Nusselt= free(Ra)

    %Calculating the Nusselt Number
        if Ra<=10000
            Nusselt=disp('invalid')
        elseif Ra<=10^7
            Nusselt= 0.54*Ra^0.25
        elseif Ra<=10^10
            Nusselt= 0.15*Ra^(1/3)
        else Ra>10^10
            Nusselt= disp('invalid')
        end 
    end 

    %Forced Convection Function (Lecture 5, EBS 218)
    function Nusselt= forced(Re,Pr)

        if Re<20000
            Nusselt= disp('invalid');
        elseif Re<=90000
            Nusselt=0.86*Re^0.5*Pr^(1/3);
        else Re>90000
            Nusselt=disp('invalid');
        end 
    end 

    %Mixed Convection Function
    function Nusselt=mixed(forced,free)
        Nusselt=(forced(Ra))^3+(free(Re,Pr))^3)^(1/3); 
    end 

    
    
    %Calculating the Convective Heat Transfer Coefficient (Lecture 5, EBS
    %218)
    hc= Nusselt*Kappa/L
end
    
    











