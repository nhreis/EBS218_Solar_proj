function [ SolZen , SolAz, SolAlt, SolDec, solT, hrAng, indAng ] = WhereIsSun( pos, t, tilt, SurAz )
%Input formats for position and time
%pos = [0, 0]; %position of Earth's surface in [deg Latitude, deg Longitude]
%t = [0,0]; % time in [day of year (of 365), mintues (standard clock time)]

%local meridian for Davis is 120 degrees
merid = 120;

%Calculate Solar Declination Angle
SolDec = 23.45*sind(360*(284+t(1))/365);


%Calculate solar time (solT)

    %find Equation of Time (EoT)
    B=(t(1) -1)*(360/365);
    EoT = 229.2*(0.000075+0.001868*cosd(B) - 0.032077*sind(B) - 0.014615*cosd(2*B) - 0.04089*sind(2*B));

    
solT = 4*(merid-pos(2))+ EoT + t(2);



%Calculate hour angle
hrAng = 15*(solT/60 -12);

%Calculate Solar Zenith Angle
SolZen = cosd(pos(1))*cosd(SolDec)*cosd(hrAng) + sind(pos(1))*sind(SolDec);
SolZen = acosd(SolZen); 
%Calculate Solar Altitude Angle
SolAlt = 90-SolZen;


%Calculate Solar Azimuth Angle 
A = asind((sind(hrAng)*cosd(SolDec))/sind(SolZen));

if cosd(SolZen)>= sind(SolDec)/sind(pos(1))
    SolAz = A;
else if cosd(SolZen)< sind(SolDec)/sind(pos(1))
        if A<0
            SolAz = -A-180;
            
        else if A>=0
                SolAz = 180 -A;
            end
        end
    end
end


%%b) Find Direction to sun from any surface

indAng = sind(SolDec)*sind(pos(1))*cosd(tilt) - sind(SolDec)*cosd(pos(1))*sind(tilt)*cosd(SurAz) + cosd(SolDec)*cosd(pos(1))*cosd(tilt)*cosd(hrAng) + cosd(SolDec)*sind(pos(1))*sind(tilt)*cosd(SurAz)*cosd(hrAng) + cosd(SolDec)*sind(tilt)*sind(SurAz)*sind(hrAng);
%indAng = cosd(SolZen)*cosd(tilt)+sind(SolAz)*sind(tilt)*cosd(SolAz-SurAz);
%The simplified expression.
indAng = acosd(indAng);




end
