function [data] = getNORA3(targetLat,targetLon,targetYear,targetMonth,targetDay,targetHour,resolution,varargin)
%
% [data] = getNORA3(targetLat,targetLon,targetYear,targetMonth,targetDay,targetHour,resolution)
% reads the data from the NORA3 atmosphere hindcast and resample them as
% gridded data. The data source is available on https://thredds.met.no/.
%
% Input:
%    * targetLat: [1x1] or [1xNlat] double: target latitude where the data
%         should be extracted
%    * targetLon: [1x1] or [1xNlon] double: target longitude where the data
%         should be extracted
%    * targetYear: [1x1] double: year(s) from which the data are
%         taken
%    * targetMonth: [1x1]  double: month(s) from which the data are
%         taken
%    * targetDay: [1x1] double: day(s) from which the data are
%         taken
%    * targetHour: [1x1] double: hour(s) from which the data are
%         taken
%    * resolution: [1x1]:double: resolution, in degree
%    *  varargin: 
%               -'optPara'cell of string: Name of the variables to read and
%               extract the Temperature data from the netcdf file.
%               - 'speedup': double: 1 or anything else. If 1 is chosen,
%               use a faster method, but I have no guarantee that it works
%               100% of the time
% Outputs:
%   * data: structure with the following fields
%       - time: [1x1] datetime
%       - lat:  [Nlon x Nlat] double: gridded latitude
%       - lon:  [Nlon x Nlat] double: gridded longitude
%       - U:  [Nlon x Nlat x Nz] double: wind_speed at Nz heights above the surface in
%       m/s
%       - DD:  [Nlon x Nlat x Nz] double: wind_direction at Nz heights above the surface
%       in deg
%       - D10:  [Nlon x Nlat x Nz] double: reference wind_direction at 10m above the surface
%       in deg
%       - Ue:  [Nlon x Nlat x Nz] double: eastern mean wind speed component
%       - Un:  [Nlon x Nlat x Nz] double: northern mean wind speed component
%       - T0:  [Nlon x Nlat] double: Surface temperature (if asked by user)
%       - P0:  [Nlon x Nlat] double: Surface pressure (if asked by user)
%       - T:  [Nlon x Nlat x Np] double: Temperture profiles (if asked by user)
%       - zT:  [Nlon x Nlat x Np] double: height for the temperture profiles (if asked by user)
%       - p:  [Np x 1] double: pressure levels for the temeprature (if asked by user)
%       - z:  [1 x Nz] double: height above the surface
% Author: E. Cheynet - UiB, Norway - last modified: 07-09-2021


%% Optional aprameters


p = inputParser();
p.CaseSensitive = false;
p.addOptional('optPara',{}); % optional aprameters
p.addOptional('speedup',1); % optional aprameters
p.parse(varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%
optPara = p.Results.optPara ;
speedup= p.Results.speedup ;


if ~any(contains(optPara,'T')) && ~isempty(optPara)
    warning('At the moment, only the temperature profiles are implemented as ''optional parameters'' ')
end
%% Check month and day number + transform date into string
[myYear,myMonth,myDay,myFolder,leapTime,targetDate] = getMyDate(targetYear,targetMonth,targetDay,targetHour);
%% Preallocation and initalisation
data = struct('time',[],'U',[],'D',[],'Un',[],'Ue',[],'lon',[],'lat',[],'D10',[]);
urldat= ['https://thredds.met.no/thredds/dodsC/nora3/',myYear,'/',myMonth,'/',myDay,'/',myFolder,'/fc',myYear,myMonth,myDay,myFolder,'_00',leapTime,'_fp.nc'];
time0 = ncread(urldat,'time')./86400+datenum('1970-01-01 00:00:00');
data.time = datetime(datestr(double(time0)));
if data.time ~= targetDate,    error('Failure to recover the target date');end


lon00 = ncread(urldat,'longitude');
lat00 = ncread(urldat,'latitude');
[N1,N2]=size(lon00);

% height_above_msl= ncread(urldat,'height_above_msl');
% top_of_atmosphere= ncread(urldat,'top_of_atmosphere');
% blh= ncread(urldat,'atmosphere_boundary_layer_thickness');
%% Definition of the grid
[lon,lat] = meshgrid(targetLon(1):resolution:targetLon(end),targetLat(1):resolution:targetLat(end));
[Nlat,Nlon]=size(lat);

%% To speed up the data extraction in one location
if numel(targetLat)==1 && numel(targetLon)==1 && speedup ==1
    
    [~,indStart] = min(sqrt((lat00(:)-(targetLat(:)-0.5)).^2 + abs(lon00(:)-(targetLon(:)-0.5)).^2));
    [~,indEnd] = min(sqrt((lat00(:)-(targetLat(:)+0.5)).^2 + abs(lon00(:)-(targetLon(:)+0.5)).^2));
    [row1, col1] = ind2sub(size(lat00), indStart);
    [row2, col2] = ind2sub(size(lat00), indEnd);
    
    r1 = min(row1,row2); % row start
    cr = abs(row2-row1);% row count
    
    c1 = min(col1,col2); % column start
    cc = abs(col2-col1); % column count
    
    
    lon00 = ncread(urldat,'longitude',[r1,c1],[cr,cc]);
    lat00 = ncread(urldat,'latitude',[r1,c1],[cr,cc]);
    [N1,N2]=size(lon00);
    height2 = ncread(urldat,'height2');
    height4 = ncread(urldat,'height4');
    z = [height4;height2];
    Nz = numel(z);
    ux = zeros(N1,N2,Nz);
    uy = zeros(N1,N2,Nz);
    ux(:,:,2:end) = ncread(urldat,'x_wind_z',[r1,c1,1,1],[cr,cc,Nz-1,1]); % zonal
    uy(:,:,2:end) = ncread(urldat,'y_wind_z',[r1,c1,1,1],[cr,cc,Nz-1,1]); % meridional
    ux(:,:,1) = ncread(urldat,'x_wind_10m',[r1,c1,1,1],[cr,cc,1,1]);
    uy(:,:,1) = ncread(urldat,'y_wind_10m',[r1,c1,1,1],[cr,cc,1,1]);
    dir10 = ncread(urldat,'wind_direction',[r1,c1,1,1],[cr,cc,1,1]);
    
    if any(contains(optPara,'T'))
        Np = 16;
        p = ncread(urldat,'pressure0');
        
        T = ncread(urldat,'air_temperature_pl',[r1,c1,1,1],[cr,cc,Np,1]);
        T0 = ncread(urldat,'air_temperature_0m',[r1,c1,1,1],[cr,cc,1,1]); % surface temeprature
        P0 = ncread(urldat,'surface_air_pressure',[r1,c1,1,1],[cr,cc,1,1]); % surface_air_pressure

        
    end
else
    
    h_10m = ncread(urldat,'height4');
    h = ncread(urldat,'height2');
    z = [h_10m;h];
    
    Nz = numel(z);
    
    ux = zeros(N1,N2,Nz);
    uy = zeros(N1,N2,Nz);
    ux(:,:,2:end) = ncread(urldat,'x_wind_z'); % zonal
    uy(:,:,2:end) = ncread(urldat,'y_wind_z'); % meridional
    ux(:,:,1) = ncread(urldat,'x_wind_10m');
    uy(:,:,1) = ncread(urldat,'y_wind_10m');
    dir10 = ncread(urldat,'wind_direction');
    
    if any(contains(optPara,'T'))
        Np = 16;
        p = ncread(urldat,'pressure0'); % zonal
        T = ncread(urldat,'air_temperature_pl'); % zonal
        T0 = ncread(urldat,'air_temperature_0m'); % surface temeprature
        P0 = ncread(urldat,'surface_air_pressure'); % surface_air_pressure
    end
end
%% Read the data in a for loop for each selected output
%  Data are resampled spatially as gridded data
dummyLat = lat00(:);
dummyLon = lon00(:);
ind = find(dummyLat>=min(targetLat(:)-0.1) & dummyLat <= max(targetLat(:)+0.1) & dummyLon>=min(targetLon(:)-0.1) & dummyLon <= max(targetLon(:)+0.1));

newUx = zeros(Nlat,Nlon,Nz);
newUy = zeros(Nlat,Nlon,Nz);

if any(contains(optPara,'T'))
    newT = zeros(Nlat,Nlon,Np);
    for ii=1:Np
        dummyT = T(:,:,ii);
        dummyT = dummyT(:);
        F_T = scatteredInterpolant(dummyLat(ind),dummyLon(ind),dummyT(ind),'linear','none');
        newT(:,:,ii) = F_T(lat,lon);
    end
    
    dummyT0 = T0(:);
    F_T0 = scatteredInterpolant(dummyLat(ind),dummyLon(ind),dummyT0(ind),'linear','none');
    newT0 = F_T0(lat,lon);
    
    dummyP0 = P0(:);
    F_P0 = scatteredInterpolant(dummyLat(ind),dummyLon(ind),dummyP0(ind),'linear','none');
    newP0 = F_P0(lat,lon);
    
end

clear F_*

for ii=1:Nz
    if ii==1
        H = sqrt(ux(:,:,ii).^2 + uy(:,:,ii).^2);
        Un = H.*cosd(dir10); Un = Un(:);
        Ue = H.*sind(dir10); Ue = Ue(:);
        F_un = scatteredInterpolant(dummyLat(ind),dummyLon(ind),Un(ind),'linear','none');
        F_ue = scatteredInterpolant(dummyLat(ind),dummyLon(ind),Ue(ind),'linear','none');
        newUn = F_un(lat,lon);
        newUe = F_ue(lat,lon);
        newDir = atan2d(newUe,newUn);
    end
    
    dummyUx = ux(:,:,ii); dummyUx = dummyUx(:);
    dummyUy = uy(:,:,ii); dummyUy = dummyUy(:);
    F_ux = scatteredInterpolant(dummyLat(ind),dummyLon(ind),dummyUx(ind),'linear','none');
    F_uy = scatteredInterpolant(dummyLat(ind),dummyLon(ind),dummyUy(ind),'linear','none');
    newUx(:,:,ii) = F_ux(lat,lon);
    newUy(:,:,ii) = F_uy(lat,lon);
end

newDir(newDir>360)= newDir(newDir>360) -360;
newDir(newDir<0)= newDir(newDir<0)  + 360;

data.D10 = newDir;
data.D = atan2d(-newUx,-newUy);
data.D(data.D<0)= 360 + data.D(data.D<0);


offset = (data.D(:,:,1) - data.D10);
data.D = data.D - offset;

data.U = sqrt(newUx.^2 + newUy.^2);
data.lon = lon;
data.lat = lat;

if any(contains(optPara,'T'))
    data.T = flip(newT,3);
    data.T0 = newT0;
    data.P0 = newP0;
    data.p = double(flipud(p));
    zT = getZ_hydrostatic(data.p,newP0/100,newT0);
    data.zT = zT; % height for temperature
end



for ii=1:Nz
    data.Un(:,:,ii) = data.U(:,:,ii).*cosd(data.D(:,:,ii));
    data.Ue(:,:,ii) = data.U(:,:,ii).*sind(data.D(:,:,ii));
end
data.z = double(z);

    function [z] = getZ_hydrostatic(P,P0,T0)
        g = 9.81 ;% Earth-surface gravitational acceleration
        L = -6.5e-3;% Lapse rate
        Rp = 287.053; %  specific gas constant = 287.053 J/(kg K)
        coeff = -L*Rp./g;
        A = (P./P0).^coeff-1;
        z = T0/L.*A;
    end

    function [myYear,myMonth,myDay,myFolder,leapTime,targetDate] = getMyDate(targetYear,targetMonth,targetDay,targetHour)
        % folder 00: 3,4,5,6,7,8,9
        % folder 06: 9,10,11,12,13,14,15
        % folder 12: 15,16,17,18,19,20,21
        % folder 18: 21,22,23,00,01,02,03
        if ischar(targetYear),  targetYear = str2double(targetYear); end
        if ischar(targetMonth),  targetMonth = str2double(targetMonth); end
        if ischar(targetDay),  targetDay = str2double(targetDay); end
        if ischar(targetHour),  targetHour = str2double(targetHour); end
        
        if targetHour<3
            myFolder = '18';
            myDate = datetime(targetYear,targetMonth,targetDay)-days(1);
            
            targetYear0 = targetYear;
            targetMonth0 = targetMonth;
            targetDay0 = targetDay;
            targetYear = year(myDate);
            targetMonth = month(myDate);
            targetDay = day(myDate);
            leapTime = num2str(mod(targetHour-18,24));
            
            targetDate = datetime(targetYear0,targetMonth0,targetDay0,targetHour,0,0);
            
        elseif targetHour>=3 &&  targetHour<9
            myFolder = '00';
            leapTime = num2str(targetHour);
            targetDate = datetime(targetYear,targetMonth,targetDay,targetHour,0,0);
        elseif targetHour>=9 &&  targetHour<15
            myFolder = '06';
            leapTime = num2str(targetHour-6);
            targetDate = datetime(targetYear,targetMonth,targetDay,targetHour,0,0);
        elseif targetHour>=15 &&  targetHour<21
            myFolder = '12';
            leapTime = num2str(targetHour-12);
            targetDate = datetime(targetYear,targetMonth,targetDay,targetHour,0,0);
        elseif targetHour>=21
            myFolder = '18';
            leapTime = num2str(targetHour-18);
            targetDate = datetime(targetYear,targetMonth,targetDay,targetHour,0,0);
        end
        
        if targetMonth<10
            myMonth = ['0',num2str(targetMonth)];
        else
            myMonth = num2str(targetMonth);
        end
        if targetDay<10
            myDay = ['0',num2str(targetDay)];
        else
            myDay = num2str(targetDay);
        end
        myYear = num2str(targetYear);
        
    end



end
