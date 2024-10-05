clearvars;close all;clc;
targetDate = datetime(2004,01,01):days(1):datetime(2024,01,01);
resolution = 0.025; % resolution for latitude and longitude in degrees (< 3 km)
targetLat = [58.15]; % target latitude
targetLon = [7.85];  % target longitude
N = numel(targetDate);
Nz = 7; % heights levels
Nlat = numel(targetLat);
Nlon = numel(targetLon);
U = NaN(N*24, Nlat, Nlon, Nz); % Mean wind speed profile
D = NaN(N*24, Nlat, Nlon, Nz); % Mean wind direction profile
P0 = NaN(N*24, Nlat, Nlon); % Surface air pressure
RH0 = NaN(N*24, Nlat, Nlon); % Relative humidity
rain = NaN(N*24, Nlat, Nlon); % Precipitation
T0 = NaN(N*24, Nlat, Nlon); % Surface temperature
time = NaT(1, N*24);
for ii = 1:N
    targetYear = num2str(year(targetDate(ii)));
    targetMonth = num2str(month(targetDate(ii)));
    targetDay = num2str(day(targetDate(ii)));
    try

        [mydata] = getNORA3_subset(targetLat, targetLon, targetYear, targetMonth, targetDay, resolution, 'optPara', 'atm_1h');
        [dz,indZ]=min(abs(mydata.z-150));
        if numel(mydata.z)==8 & dz==0,
            mydata.D(:,:,indZ,:)=[];
            mydata.U(:,:,indZ,:)=[];
            mydata.Un(:,:,indZ,:)=[];
            mydata.Ue(:,:,indZ,:)=[];
            mydata.z(indZ)=[];
        end


        startIndex = (ii-1) * 24 + 1;
        endIndex = ii * 24;

        P0(startIndex:endIndex, :, :) = squeeze(mydata.P0);
        T0(startIndex:endIndex, :, :) = squeeze(mydata.T0);
        RH0(startIndex:endIndex, :, :) = squeeze(mydata.RH0);
        rain(startIndex:endIndex, :, :) = squeeze(mydata.precipitation);
        U(startIndex:endIndex, :, :, :) = squeeze(mydata.U)';
        D(startIndex:endIndex, :, :, :) = squeeze(mydata.D)';
        time(startIndex:endIndex) = mydata.time;


    catch exception

        try
        warning('error, trying again...')
        pause(20)
        
        [mydata] = getNORA3_subset(targetLat, targetLon, targetYear, targetMonth, targetDay, resolution, 'optPara', 'atm_1h');
        [dz,indZ]=min(abs(mydata.z-150));
        if numel(mydata.z)==8 & dz==0,
            mydata.D(:,:,indZ,:)=[];
            mydata.U(:,:,indZ,:)=[];
            mydata.Un(:,:,indZ,:)=[];
            mydata.Ue(:,:,indZ,:)=[];
            mydata.z(indZ)=[];
        end


        startIndex = (ii-1) * 24 + 1;
        endIndex = ii * 24;

        P0(startIndex:endIndex, :, :) = squeeze(mydata.P0);
        T0(startIndex:endIndex, :, :) = squeeze(mydata.T0);
        RH0(startIndex:endIndex, :, :) = squeeze(mydata.RH0);
        rain(startIndex:endIndex, :, :) = squeeze(mydata.precipitation);
        U(startIndex:endIndex, :, :, :) = squeeze(mydata.U)';
        D(startIndex:endIndex, :, :, :) = squeeze(mydata.D)';
        time(startIndex:endIndex) = mydata.time;

        catch exception
            warning('Data cannot be read')
            fprintf([exception.message,' \n'])
        end
    end


end


save('NORA3_Kristiansand_2004_2024.mat')

T2m = squeeze(T0);
time_datenum = datenum(time);
u10=squeeze(U(:,1,1,1));
save('NORA3_Kristiansand_2023_2024_python.mat','T2m','time_datenum','u10','targetLat','targetLon','resolution')

%% Plot Mean Wind Speed at the first height level
clearvars;close all;clc;
load('NORA3_Kristiansand_2004_2024.mat')
T2m = squeeze(T0);
time_datenum = datenum(time);
u10=squeeze(U(:,1,1,1));

figure;
plot(time, u10);  % Assuming Nlat = 1 and Nlon = 1 for simplicity
title('Mean Wind Speed at First Height Level');
xlabel('Time');
ylabel('Mean Wind Speed (m/s)');
grid on;

%% Plot Air Temperature at the surface
figure;
plot(time, T2m);
title('Air Temperature at Surface');
xlabel('Time');
ylabel('Air Temperature (K)');
grid on;

%% Plot Precipitation
figure;
y = squeeze(rain(:,1,1));
y(y<0) = 0;
plot(time, y);
title('Precipitation');
xlabel('Time');
ylabel('Precipitation (mm)');
grid on;

%% Plot Relative Humidity at the surface
figure;
plot(time, squeeze(RH0(:,1,1)));
title('Relative Humidity at Surface');
xlabel('Time');
ylabel('Relative Humidity (%)');
grid on;

%% Plot Mean Wind Direction at the first height level
figure;
plot(time, squeeze(D(:,1,1,1)));  % Again, assuming Nlat = 1 and Nlon = 1
title('Mean Wind Direction at First Height Level');
xlabel('Time');
ylabel('Mean Wind Direction (degrees)');
grid on;

%% Convert the data into nc file

% Load the .mat file
data=load('NORA3_Kristiansand_2023_2024_python.mat');

% Define the NetCDF file name
ncfile = 'NORA3ForFirda_Kristiansand_2004_2024.nc';

% Convert time_datenum to seconds since 1970-01-01
origin = datenum('1970-01-01 00:00:00');
time_seconds = (data.time_datenum - origin) * 86400; % Convert days to seconds

% Create dimensions
nccreate(ncfile, 'time_datenum', 'Dimensions', {'time', length(data.time_datenum)}, 'Datatype', 'double');
nccreate(ncfile, 'time_seconds', 'Dimensions', {'time', length(data.time_datenum)}, 'Datatype', 'double');
nccreate(ncfile, 'u10', 'Dimensions', {'time', length(data.time_datenum)}, 'Datatype', 'double');
nccreate(ncfile, 'T2m', 'Dimensions', {'time', length(data.time_datenum)}, 'Datatype', 'double');
nccreate(ncfile, 'targetLat', 'Dimensions', {'scalar', 1}, 'Datatype', 'double');
nccreate(ncfile, 'targetLon', 'Dimensions', {'scalar', 1}, 'Datatype', 'double');
nccreate(ncfile, 'resolution', 'Dimensions', {'scalar', 1}, 'Datatype', 'double');

% Write variables to the NetCDF file
ncwrite(ncfile, 'time_datenum', data.time_datenum);
ncwrite(ncfile, 'time_seconds', time_seconds);
ncwrite(ncfile, 'u10', data.u10);
ncwrite(ncfile, 'T2m', data.T2m);
ncwrite(ncfile, 'targetLat', data.targetLat);
ncwrite(ncfile, 'targetLon', data.targetLon);
ncwrite(ncfile, 'resolution', data.resolution);

% Add attributes to variables
ncwriteatt(ncfile, 'time_datenum', 'units', 'days since 0000-01-01 00:00:00');
ncwriteatt(ncfile, 'time_seconds', 'units', 'seconds since 1970-01-01 00:00:00');
ncwriteatt(ncfile, 'u10', 'long_name', 'Mean wind speed at 10m above surface');
ncwriteatt(ncfile, 'u10', 'units', 'm/s');
ncwriteatt(ncfile, 'T2m', 'long_name', 'Air temperature at 2m above surface');
ncwriteatt(ncfile, 'T2m', 'units', 'K'); % Assuming temperature is in Kelvin
ncwriteatt(ncfile, 'targetLat', 'units', 'degrees_north');
ncwriteatt(ncfile, 'targetLon', 'units', 'degrees_east');
ncwriteatt(ncfile, 'resolution', 'long_name', 'Spatial resolution of the data interpolation');
ncwriteatt(ncfile, 'resolution', 'units', 'degrees');

disp(['NetCDF file "' ncfile '" created successfully using nccreate and ncwrite.']);

%% read the .nc file
clearvars;
ncfile = 'NORA3ForFirda_Kristiansand_2004_2024.nc';
% Read the variables from the NetCDF file
time_datenum = ncread(ncfile, 'time_datenum');
time = seconds(ncread(ncfile, 'time_seconds'))+datetime(1970,1,1);
u10 = ncread(ncfile, 'u10');
T2m = ncread(ncfile, 'T2m');
targetLat = ncread(ncfile, 'targetLat');
targetLon = ncread(ncfile, 'targetLon');
resolution = ncread(ncfile, 'resolution');
