function [data] = getNORA3_subset(targetLat,targetLon,targetYear,targetMonth,targetDay,resolution,varargin)
%
% [data] = getNORA3_subset(targetLat, targetLon, targetYear, targetMonth, targetDay, resolution, varargin)
% Reads data from the NORA3 atmosphere hindcast, resampling them as
% gridded data. This data source is available on https://thredds.met.no/.
%
% Input:
%    * targetLat: [1x1] or [1xNlat] double: target latitude(s) where the data
%         should be extracted.
%    * targetLon: [1x1] or [1xNlon] double: target longitude(s) where the data
%         should be extracted.
%    * targetYear: [1x1] string: year from which the data are
%         taken.
%    * targetMonth: [1x1] string: month from which the data are
%         taken.
%    * targetDay: [1x1] string: day from which the data are
%         taken.
%    * resolution: [1x1] double: resolution in degrees (latitude and longitude).
%    * varargin: additional parameters including:
%               - 'optPara': cell of strings specifying the variables to be read,
%                 e.g., 'atm_1h' to specify the temporal resolution or type of atmospheric data.
%
% Outputs:
%   * data: structure with the following fields:
%       - time: [N*24 x 1] datetime: times corresponding to the extracted data points.
%       - lat:  [Nlat x Nlon] double: gridded latitude values.
%       - lon:  [Nlat x Nlon] double: gridded longitude values.
%       - U:  [N*24 x Nlat x Nlon x Nz] double: Mean wind speed at Nz heights above the surface in m/s.
%       - D:  [N*24 x Nlat x Nlon x Nz] double: Mean wind direction at Nz heights above the surface in degrees.
%       - P0: [N*24 x Nlat x Nlon] double: Surface air pressure in Pascals.
%       - RH0: [N*24 x Nlat x Nlon] double: Relative humidity at the surface in percentage.
%       - T0:  [N*24 x Nlat x Nlon] double: Surface air temperature in Kelvin.
%       - rain: [N*24 x Nlat x Nlon] double: Precipitation amount in mm.
%       - z:  [1 x Nz] double: Heights above the surface for wind measurements.
%
% Author: E. Cheynet - UiB, Norway - last modified: 2024-05-13


%% Optional aprameters
p = inputParser();
p.CaseSensitive = false;
p.addOptional('optPara',{}); % optional aprameters
p.parse(varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%
optPara = p.Results.optPara ;

if numel(targetDay)==1, targetDay = ['0',targetDay];end
if numel(targetMonth)==1, targetMonth = ['0',targetMonth];end
    

%% Preallocation and initalisation
data = struct('time',[],'U',[],'D',[],'Un',[],'Ue',[],'lon',[],'lat',[]);

[myYear,myMonth,myDay,targetDate] = getMyDate(targetYear,targetMonth,targetDay);
urldat = ['https://thredds.met.no/thredds/dodsC/nora3_subset_atmos/wind_hourly/arome3kmwind_1hr_',myYear,myMonth,'.nc'];

time0 = seconds(ncread(urldat,'time')) +datetime('1970-01-01 00:00:00');

[~,indTime]=min(abs(time0-targetDate));

data.time = time0(indTime:indTime+24-1);

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

if numel(targetLat)==1 &&   numel(targetLon)==1
    offsetLat = 0.3;
    offsetLon = 0.15;
    [~,indStart] = min(sqrt((lat00(:)-(targetLat(1)-offsetLat)).^2 + abs(lon00(:)-(targetLon(1)-offsetLon)).^2));
    [~,indEnd] = min(sqrt((lat00(:)-(targetLat(1)+offsetLat)).^2 + abs(lon00(:)-(targetLon(1)+offsetLon)).^2));
elseif numel(targetLat)==2 &&   numel(targetLon)==2
    offsetLat = 0.5*sqrt(diff(targetLat).^2 + diff(targetLon).^2);
    offsetLon = offsetLat;
    [~,indStart] = min(sqrt((lat00(:)-(targetLat(1)-offsetLat)).^2 + abs(lon00(:)-(targetLon(2)-offsetLon)).^2));
    [~,indEnd] = min(sqrt((lat00(:)-(targetLat(2)+offsetLat)).^2 + abs(lon00(:)-(targetLon(1)+offsetLon)).^2));
else
    error('targetLat and targetLatmust have the same dimensions')
end

[row1, col1] = ind2sub(size(lat00), indStart);
[row2, col2] = ind2sub(size(lat00), indEnd);
r1 = min(row1,row2); % row start
cr = abs(row2-row1+1);% row count
c1 = min(col1,col2); % column start
cc = abs(col2-col1+1); % column count
lon00 = ncread(urldat,'longitude',[r1,c1],[cr,cc]);
lat00 = ncread(urldat,'latitude',[r1,c1],[cr,cc]);
dummyLat = double(lat00(:));
dummyLon = double(lon00(:));
ind = find(dummyLat>=min(targetLat(:)-offsetLat) & dummyLat <= max(targetLat(:)+offsetLat) &...
    dummyLon>=min(targetLon(:)-offsetLon) & dummyLon <= max(targetLon(:)+offsetLon));

z = ncread(urldat,'height');
Nz = numel(z);
meanU = ncread(urldat,'wind_speed',[r1,c1,1,indTime],[cr,cc,Nz,24]); % zonal
windDir = ncread(urldat,'wind_direction',[r1,c1,1,indTime],[cr,cc,Nz,24]);

if any(contains(optPara,'atm_1h'))
    urldat2 = ['https://thredds.met.no/thredds/dodsC/nora3_subset_atmos/atm_hourly/arome3km_1hr_',myYear,myMonth,'.nc'];
    
    rain = ncread(urldat2,'precipitation_amount_hourly',[r1,c1,1,indTime],[cr,cc,1,24]);% precipitation_amount_hourly
    RH0 = ncread(urldat2,'relative_humidity_2m',[r1,c1,1,indTime],[cr,cc,1,24]); % relative humidity at 2 m
    T0 = ncread(urldat2,'air_temperature_2m',[r1,c1,1,indTime],[cr,cc,1,24]); % air temperature at 2 m
    P0 = ncread(urldat2,'air_pressure_at_sea_level',[r1,c1,1,indTime],[cr,cc,1,24]); % surface_air_pressure
end
%% Read the data in a for loop for each selected output
%  Data are resampled spatially as gridded data

N = numel(data.time);
[Nlat, Nlon] = size(lat);

newR = zeros(Nlat, Nlon, N);
newRH = zeros(Nlat, Nlon, N);
newT0 = zeros(Nlat, Nlon, N);
newP0 = zeros(Nlat, Nlon, N);

if any(contains(optPara,'atm_1h'))
    
    % Loop over each time step
    for jj = 1:N
        % Rain
        dummyRain = reshape(rain(:,:,:,jj),[],1);
        F_R = scatteredInterpolant(dummyLat(ind),dummyLon(ind),dummyRain(ind),'linear','none');
        newR(:,:,jj) = F_R(lat,lon);
        
        % Relative Humidity
        dummyRH = reshape(RH0(:,:,jj),[],1);
        F_RH = scatteredInterpolant(dummyLat(ind),dummyLon(ind),dummyRH(ind),'linear','none');
        newRH(:,:,jj) = F_RH(lat, lon);
        
        % Surface Temperature
        dummyT0 = reshape(T0(:,:,jj),[],1);
        F_T0 = scatteredInterpolant(dummyLat(ind),dummyLon(ind),dummyT0(ind),'linear','none');
        newT0(:,:,jj) = F_T0(lat, lon);
        
        % Surface Pressure
        dummyP0 = reshape(P0(:,:,jj),[],1);
        F_P0 = scatteredInterpolant(dummyLat(ind),dummyLon(ind),dummyP0(ind),'linear','none');
        newP0(:,:,jj) = F_P0(lat, lon);
    end
end


clear F_*

% Initialize the matrices for Un, Ue, and wind direction D
data.Un = zeros(Nlat, Nlon, Nz, N);
data.Ue = zeros(Nlat, Nlon, Nz, N);
data.D = zeros(Nlat, Nlon, Nz, N);
data.U = zeros(Nlat, Nlon, Nz, N);

if any(contains(optPara,'atm_1h'))
    
    for ii = 1:Nz
        for jj = 1:N
            % Calculate Un and Ue at each grid point
            H = meanU(:,:,ii,jj);
            Un = H .* cosd(windDir(:,:,ii,jj));
            Un = Un(:);
            Ue = H .* sind(windDir(:,:,ii,jj));
            Ue = Ue(:);
            
            % Create scattered interpolants for Un and Ue
            F_un = scatteredInterpolant(dummyLat(ind),dummyLon(ind),Un(ind),'linear','none');
            F_ue = scatteredInterpolant(dummyLat(ind),dummyLon(ind),Ue(ind),'linear','none');
            
            % Interpolate values at the lat, lon grid
            newUn = F_un(lat,lon);
            newUe = F_ue(lat,lon);
            
            % Calculate new wind direction
            newDir = atan2d(newUe, newUn);
            % Adjust wind directions to be within [0, 360] degrees
            newDir(newDir > 360) = newDir(newDir > 360) - 360;
            newDir(newDir < 0) = newDir(newDir < 0) + 360;
            
            % Store interpolated values and derived direction
            data.Un(:,:,ii,jj) = newUn;
            data.Ue(:,:,ii,jj) = newUe;
            data.D(:,:,ii,jj) = newDir;
            data.U(:,:,ii,jj) = sqrt(newUn.^2 + newUe.^2);
        end
    end
end


data.lon = lon;
data.lat = lat;

if any(contains(optPara,'atm_1h'))
    data.T0 = newT0;
    data.P0 = newP0;
    data.RH0 = newRH;
    data.precipitation = newR;
end

data.z = double(z);
% data.time = time;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [myYear, myMonth, myDay, targetDate] = getMyDate(targetYear, targetMonth, targetDay)
        % Convert input to numeric if they are strings
        if ischar(targetYear), targetYear = str2double(targetYear); end
        if ischar(targetMonth), targetMonth = str2double(targetMonth); end
        if ischar(targetDay), targetDay = str2double(targetDay); end

        % Create the target date as a datetime object
        targetDate = datetime(targetYear, targetMonth, targetDay, 0, 0, 0);

        % Get the date components from the day before the target date
        myYear = num2str(year(targetDate));
        myMonth = sprintf('%02d', month(targetDate));
        myDay = sprintf('%02d', day(targetDate));
    end



end
