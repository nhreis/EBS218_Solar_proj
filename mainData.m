clear all
data = xlsread('mainData.xlsx'); % reads collector data from Excel file
start = 292;    % row # where useful data begins

% Day 150 Data
fin150 = 454;   % row # of last day150 data 
day150.day = data(start:fin150,2);  % day #
day150.day(day150.day==151)=150;    % set 151s to 150

day150.hour = data(start:fin150,3); % hour
day150.min = data(start:fin150,4);  % minute
day150.T_air_In = data(start:fin150,8);  % air temp at inlet (deg C)
day150.T_air_Out = data(start:fin150,9);   % air temp at outlet (deg C)
day150.T_back = data(start:fin150,10); % back plate temp (deg C)

day150.Gt = data(start:fin150,11); % solar radiation on collector (W/m2)
day150.Gt(day150.Gt==-6999)=NaN;    % remove erroneous data
day150.Gt(day150.Gt<0)=0;   % sets negative values to zero

day150.m_Dot = data(start:fin150,15);  % mass flow rate of air (kg/s)
day150.airT = data(start:fin150,16);    % air temp (deg C)

day150.T_air_In(day150.T_air_In<-6000)=NaN;
day150.Gt(1)=900;

% Day 151 Data
start151 = fin150+1;  % day #
fin151 = 747;   % row # of last day151 data
day151.day = data(start151:fin151,2);   % day #
day151.day(day151.day==152)=151;    % set 152s to 151

day151.hour = data(start151:fin151,3);  % hour
day151.min = data(start151:fin151,4);   % minute
day151.T_air_In = data(start151:fin151,8);  % air temp at inlet (deg C)
day151.T_air_Out = data(start151:fin151,9);   % air temp at outlet (deg C)
day151.T_back = data(start151:fin151,10); % back plate temp (deg C)

day151.Gt = data(start151:fin151,11); % solar radiation on collector (W/m2)
day151.Gt(day151.Gt==-6999)=NaN;    % remove erroneous data
day151.Gt(day151.Gt<0)=0;   % sets negative values to zero

day151.m_Dot = data(start151:fin151,15);  % mass flow rate of air (kg/s)
day151.airT = data(start151:fin151,16);    % air temp (deg C)
day151.T_air_In(day151.T_air_In<-6000)=NaN;
day151.T_air_In(147)=day151.T_air_In(147)+14;
day151.T_air_In(145)=day151.T_air_In(144);
day151.T_air_In(146)=day151.T_air_In(145);
day151.T_back(day151.T_back<-5000)=0;



% Day 152 Data
start152 = fin151+1;  % day #
fin152 = 1035;   % row # of last day152 data
day152.day = data(start152:fin152,2);   % day #
day152.day(day152.day==153)=152;    % set 153s to 152

day152.hour = data(start152:fin152,3);  % hour
day152.min = data(start152:fin152,4);   % minute
day152.T_air_In = data(start152:fin152,8);  % air temp at inlet (deg C)
day152.T_air_Out = data(start152:fin152,9);   % air temp at outlet (deg C)
day152.T_back = data(start152:fin152,10); % back plate temp (deg C)

day152.Gt = data(start152:fin152,11); % solar radiation on collector (W/m2)
day152.Gt(day152.Gt==-6999)=NaN;    % remove erroneous data
day152.Gt(day152.Gt<0)=0;   % sets negative values to zero

day152.m_Dot = data(start152:fin152,15);  % mass flow rate of air (kg/s)
day152.airT = data(start152:fin152,16);    % air temp (deg C)


% Day 153 Data
start153 = fin152+1;  % day #
fin153 = 1323;   % row # of last day152 data
day153.day = data(start153:fin153,2);   % day #
day153.day(day153.day==154)=153;    % set 154s to 153

day153.hour = data(start153:fin153,3);  % hour
day153.min = data(start153:fin153,4);   % minute
day153.T_air_In = data(start153:fin153,8);  % air temp at inlet (deg C)
day153.T_air_Out = data(start153:fin153,9);   % air temp at outlet (deg C)
day153.T_back = data(start153:fin153,10); % back plate temp (deg C)

day153.Gt = data(start153:fin153,11); % solar radiation on collector (W/m2)
day153.Gt(day153.Gt==-6999)=NaN;    % remove erroneous data
day153.Gt(day153.Gt<0)=0;   % sets negative values to zero

day153.m_Dot = data(start153:fin153,15);  % mass flow rate of air (kg/s)
day153.airT = data(start153:fin153,16);    % air temp (deg C)


% Day 154 Data
start154 = fin153+1;  % day #
fin154 = 1611;   % row # of last day152 data
day154.day = data(start154:fin154,2);   % day #
day154.day(day154.day==155)=154;    % set 155s to 154

day154.hour = data(start154:fin154,3);  % hour
day154.min = data(start154:fin154,4);   % minute
day154.T_air_In = data(start154:fin154,8);  % air temp at inlet (deg C)
day154.T_air_Out = data(start154:fin154,9);   % air temp at outlet (deg C)
day154.T_back = data(start154:fin154,10); % back plate temp (deg C)

day154.Gt = data(start154:fin154,11); % solar radiation on collector (W/m2)
day154.Gt(day154.Gt==-6999)=NaN;    % remove erroneous data
day154.Gt(day154.Gt<0)=0;   % sets negative values to zero

day154.m_Dot = data(start154:fin154,15);  % mass flow rate of air (kg/s)
day154.airT = data(start154:fin154,16);    % air temp (deg C)


% Day 155 Data
start155 = fin154+1;  % day #
fin155 = 1899;   % row # of last day152 data
day155.day = data(start155:fin155,2);   % day #
day155.day(day155.day==156)=155;    % set 156s to 155

day155.hour = data(start155:fin155,3);  % hour
day155.min = data(start155:fin155,4);   % minute
day155.T_air_In = data(start155:fin155,8);  % air temp at inlet (deg C)
day155.T_air_Out = data(start155:fin155,9);   % air temp at outlet (deg C)
day155.T_back = data(start155:fin155,10); % back plate temp (deg C)

day155.Gt = data(start155:fin155,11); % solar radiation on collector (W/m2)
day155.Gt(day155.Gt==-6999)=NaN;    % remove erroneous data
day155.Gt(day155.Gt<0)=0;   % sets negative values to zero

day155.m_Dot = data(start155:fin155,15);  % mass flow rate of air (kg/s)
day155.airT = data(start155:fin155,16);    % air temp (deg C)


% Day 156 Data
start156 = fin155+1;  % day #
fin156 = 2123;   % row # of last day152 data
day156.day = data(start156:fin156,2);   % day #

day156.hour = data(start156:fin156,3);  % hour
day156.min = data(start156:fin156,4);   % minute
day156.T_air_In = data(start156:fin156,8);  % air temp at inlet (deg C)
day156.T_air_Out = data(start156:fin156,9);   % air temp at outlet (deg C)
day156.T_back = data(start156:fin156,10); % back plate temp (deg C)

day156.Gt = data(start156:fin156,11); % solar radiation on collector (W/m2)
day156.Gt(day156.Gt==-6999)=NaN;    % remove erroneous data
day156.Gt(day156.Gt<0)=0;   % sets negative values to zero

day156.m_Dot = data(start156:fin156,15);  % mass flow rate of air (kg/s)
day156.airT = data(start156:fin156,16);    % air temp (deg C)
