TrajectoryFile = '../Input_data/Trajectory.txt';
FlightLinesFile = '../Input_data/FlightLines.txt';
AnomalyFile2100 = '../Input_data/anomaly/XGM2019e_2159_2100.dat';
AnomalyFile400 = '../Input_data/anomaly/XGM2019e_2159_400.dat';
OuputAnomalyFile = '../Output_data/anomaly_data.txt';
OuputAnomalyLowFile = '../Output_data/anomaly_low_data.txt';
OuputAnomalySatFile = '../Output_data/anomaly_sat_data.txt';
IMUFile = '../Output_data/IMU_data.txt';
GPSFile = '../Output_data/GPS_data.txt';

AllGals         = 6; %кол-во нужных галсов (1, 3, 5, 7)
Add_Noise       = 1; %добавить  шум
Add_Bias        = 1; %добавить смещения и дрейфы
Add_Anomal      = 1; %добавить аномалию

StartTime       = 0; %время начала записи в секунах


AccBias         = [30, -40, 0] * 10^(-5); % [m/s/s]
GyroDrift       = [-3, 3, 1] * 10^(-3); % [deg/h]

StdVel10Hz      = 0.05; % [m/s] 10Hz
StdGyro100Hz    = 3; %[deg/h]  100Hz
StdAcc100Hz     = 30 * 10^(-5); % [m/s/s] 100Hz

TStep100Hz      = 0.01; %timestep



