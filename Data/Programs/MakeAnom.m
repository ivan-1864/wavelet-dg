% ---------------ititalization--------------------------------------------
cfg;
D              = importdata(GPSFile);
Time10Hz       = D.data(:, 1);
Data10Hz      = D.data(:, 2:end); % frequency  = 100 Hz
clear D

fid  = fopen(AnomalyFile2100); 
data2100 = textscan(fid, '%f %f %f %f %f %f %f','HeaderLines', 34);
data2100 = cell2mat(data2100);
fclose(fid);

clear fid

fid  = fopen(AnomalyFile400); 
data400 = textscan(fid, '%f %f %f %f %f %f %f','HeaderLines', 34);
data400 = cell2mat(data400);
fclose(fid);

clear fid

fid = fopen(FlightLinesFile);
C = textscan(fid,'%s %f %c %f','HeaderLines',3); 
fclose(fid);

TimeGalsBeg = C{2};
TimeGalsEnd = C{4};
clear fid C

DeltaTime = mean(diff(Time10Hz));
Freq = round(1/DeltaTime);
TimeGals = TimeGalsEnd(1) - TimeGalsBeg(1) + (TimeGalsBeg(2) - TimeGalsEnd(1));

StartInter = 0;
switch AllGals
    case 7
        EndInter = TimeGalsEnd(7) + (TimeGalsBeg(7) - TimeGalsEnd(6)) / 2 ;
    case 6
        EndInter = (TimeGalsEnd(6) + TimeGalsBeg(7)) / 2; 
    case 5
        EndInter = (TimeGalsEnd(5) + TimeGalsBeg(6)) / 2; 
    case 4
        EndInter = (TimeGalsEnd(4) + TimeGalsBeg(5)) / 2; 
    case 3
        EndInter = (TimeGalsEnd(3) + TimeGalsBeg(4)) / 2; 
    case 2
        EndInter = (TimeGalsEnd(2) + TimeGalsBeg(3)) / 2; 
    case 1
        EndInter = (TimeGalsEnd(1) + TimeGalsBeg(2)) / 2;
    otherwise
        disp('Error: check "AllGals" parametr')
end
IndexInterval = (Time10Hz(:, 1) >= StartInter) & (Time10Hz(:, 1) < EndInter);
IndexGPS = find(IndexInterval);

clear IndexInterval

LenInter = length(IndexGPS);


AnomalyAll      = zeros(LenInter, 3);
AnomalyAll_x    = zeros(LenInter,4);

AnomalyLow=zeros(LenInter, 3);
AnomalyLow_x    = zeros(LenInter,4);

% choose the interval from 2 line (from graph) 
Finish1stGals = 16130.05;
Finish2ndGals = 18549.85;
Index2ndGals=find(data2100(:, 1) >= Finish1stGals & data2100(:, 1) < Finish2ndGals);

Anomaly1Gals    = [data2100(Index2ndGals, 6),data2100(Index2ndGals, 7),data2100(Index2ndGals, 5)] ; % xi [arcsec], eta [arcsec], \Delta g [mGal]
Anomaly1Gals    =  Anomaly1Gals-Anomaly1Gals(1, :); % 1st record equal to zero
NumRec1Gals = TimeGals * Freq; % number of records in 1 line 

NewIndexAnomal=linspace(data2100(Index2ndGals(1), 1), data2100(Index2ndGals(end), 1), NumRec1Gals)';
InterpAnomalAll = interp1(data2100(Index2ndGals, 1), Anomaly1Gals, NewIndexAnomal, 'spline');

% для низкочастотной составляющей
Index2ndGals=find(data400(:, 1) >= Finish1stGals & data400(:, 1) < Finish2ndGals);

Anomaly1Gals    = [data400(Index2ndGals, 6),data400(Index2ndGals, 7),data400(Index2ndGals, 5)] ; % xi [arcsec], eta [arcsec], \Delta g [mGal]
Anomaly1Gals    =  Anomaly1Gals-Anomaly1Gals(1, :); % 1st record equal to zero
NumRec1Gals = TimeGals * Freq; % number of records in 1 line 

NewIndexAnomal=linspace(data400(Index2ndGals(1), 1), data400(Index2ndGals(end), 1), NumRec1Gals)';
InterpAnomalLow = interp1(data400(Index2ndGals, 1), Anomaly1Gals, NewIndexAnomal, 'spline');


% expand anomaly to all lines

for i = 1:AllGals

    StartImu  = TimeGalsBeg(i) - (TimeGalsBeg(i + 1) - TimeGalsEnd(i)) / 2;
    FinishImu = TimeGalsEnd(i) + (TimeGalsBeg(i + 1) - TimeGalsEnd(i)) / 2;
    
    TimeImu = find(Time10Hz(:, 1) >= StartImu & Time10Hz(:, 1) < FinishImu);
    
    if i / 2 - floor(i / 2) == 0
        AnomalyAll(TimeImu, :) = flipud(InterpAnomalAll); %  inverse array 
        AnomalyLow(TimeImu, :) = flipud(InterpAnomalLow); %  inverse array 
    else
        AnomalyAll(TimeImu, :) = InterpAnomalAll;
        AnomalyLow(TimeImu, :) = InterpAnomalLow;
    end
end


for i = 1:LenInter
   L_zx         = Make_L_zx_Matrix(deg2rad(Data10Hz(i, 9)), deg2rad(Data10Hz(i, 8)), deg2rad(Data10Hz(i, 7)));
   AnomalyAll(i, 1)   = -deg2rad(AnomalyAll(i, 1) / 3600) * Geodesy_NormalGravity(Data10Hz(i, 2), Data10Hz(i, 3)); %arcsec->m/s/s
   AnomalyAll(i, 2)   = -deg2rad(AnomalyAll(i, 2) / 3600) * Geodesy_NormalGravity(Data10Hz(i, 2), Data10Hz(i, 3)); %arcsec->m/s/s
   AnomalyAll(i, 3)   = 10^(-5) * AnomalyAll(i, 3); %mGal->m/s/s
   AnomalyAll_x(i, :) = [Time10Hz(i, 1), AnomalyAll(i, :)];
   AnomalyAll(i, :)   = AnomalyAll(i, :) * L_zx'; %Mx->Mz
          
   AnomalyLow(i, 1)   = -deg2rad(AnomalyLow(i, 1) / 3600) * Geodesy_NormalGravity(Data10Hz(i, 2), Data10Hz(i, 3)); %arcsec->m/s/s
   AnomalyLow(i, 2)   = -deg2rad(AnomalyLow(i, 2) / 3600) * Geodesy_NormalGravity(Data10Hz(i, 2), Data10Hz(i, 3)); %arcsec->m/s/s
   AnomalyLow(i, 3)   = 10^(-5) * AnomalyLow(i, 3); %mGal->m/s/s
   AnomalyLow_x(i, :) = [Time10Hz(i, 1), AnomalyLow(i, :)];
   AnomalyLow(i, :)   = AnomalyLow(i, :) * L_zx'; %Mx->Mz
end

