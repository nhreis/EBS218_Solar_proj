function [ SolZen , SolAz, SolAlt, hrAng, indAng ] = SunPosition( pos, solT, SolDec, tilt, SurAz )
%Input formats for position and time
%pos = [0, 0]; %position of Earth's surface in [deg Latitude, deg Longitude]


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
