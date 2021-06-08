function [data] = getNORA3(targetLat,targetLon,targetYear,targetMonth,targetDay,targetHour,resolution)
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
%
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
%       - z:  [1 x Nz] double: height above the surface
%
% Author: E. Cheynet - UiB, Norway - last modified: 08-06-2021

%% Definition of the grid
[lon,lat] = meshgrid(targetLon(1):resolution:targetLon(end),targetLat(1):resolution:targetLat(end));
[Nlat,Nlon]=size(lat);
%% Check month and day number + transform date into string
% folder 00: 3,4,5,6,7,8,9
% folder 06: 9,10,11,12,13,14,15
% folder 12: 15,16,17,18,19,20,21
% folder 18: 21,22,23,00,01,02,03

if ischar(targetMonth),  targetMonth = str2double(targetMonth); end
if ischar(targetDay),  targetDay = str2double(targetDay); end

if targetHour<=3
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
    
elseif targetHour>3 &&  targetHour<=9
    myFolder = '00';
    leapTime = num2str(targetHour);
    targetDate = datetime(targetYear,targetMonth,targetDay,targetHour,0,0);
elseif targetHour>9 &&  targetHour<=15
    myFolder = '06';
    leapTime = num2str(targetHour-6);
    targetDate = datetime(targetYear,targetMonth,targetDay,targetHour,0,0);
elseif targetHour>15 &&  targetHour<=21
    myFolder = '12';
    leapTime = num2str(targetHour-12);
    targetDate = datetime(targetYear,targetMonth,targetDay,targetHour,0,0);
elseif targetHour>21
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

data = struct('time',[],'U',[],'D',[],'Un',[],'Ue',[],'lon',[],'lat',[],'D10',[]);
%% Preallocation and initalisation
urldat= ['https://thredds.met.no/thredds/dodsC/nora3/',myYear,'/',myMonth,'/',myDay,'/',myFolder,'/fc',myYear,myMonth,myDay,myFolder,'_00',leapTime,'_fp.nc'];
time0 = ncread(urldat,'time')./86400+datenum('1970-01-01 00:00:00');
data.time = datetime(datestr(double(time0)));

if data.time ~= targetDate
    error('Failure to recover the target date');
end

% height_above_msl= ncread(urldat,'height_above_msl');
% top_of_atmosphere= ncread(urldat,'top_of_atmosphere');
h= ncread(urldat,'atmosphere_boundary_layer_thickness');
lon00 = ncread(urldat,'longitude');
lat00 = ncread(urldat,'latitude');
[N1,N2]=size(lon00);


height2 = ncread(urldat,'height2');
height4 = ncread(urldat,'height4');
z = [height4;height2];

Nz = numel(z);

ux = zeros(N1,N2,Nz);
uy = zeros(N1,N2,Nz);
ux(:,:,1:Nz-1) = ncread(urldat,'x_wind_z'); % zonal
uy(:,:,1:Nz-1) = ncread(urldat,'y_wind_z'); % meridional

ux(:,:,end) = ncread(urldat,'x_wind_10m');
uy(:,:,end) = ncread(urldat,'y_wind_10m');
dir10 = ncread(urldat,'wind_direction');

%% Read the data in a for loop for each selected output
%  Data are resampled spatially as gridded data
dummyLat = lat00(:);
dummyLon = lon00(:);
ind = find(dummyLat>=min(targetLat(:)-0.1) & dummyLat <= max(targetLat(:)+0.1) & dummyLon>=min(targetLon(:)-0.1) & dummyLon <= max(targetLon(:)+0.1));

newUx = zeros(Nlat,Nlon,Nz);
newUy = zeros(Nlat,Nlon,Nz);


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

for ii=1:Nz
    data.Un(:,:,ii) = data.U(:,:,ii).*cosd(data.D(:,:,ii));
    data.Ue(:,:,ii) = data.U(:,:,ii).*sind(data.D(:,:,ii));
end
data.z = double(z);
end
