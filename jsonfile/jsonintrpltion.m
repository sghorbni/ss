function [lats,lons,times] = jsonintrpltion(filename)

data = jsondecode(fileread(filename));

% Extract trails data
trails = data.trails;
timestamps = fieldnames(trails); % Get all timestamps (keys)

% Initialize arrays for lats, lons, heights, and times
lats = [];
lons = [];
times = [];
% alts = [];

% Loop through each timestamp in trails
timestamps = fieldnames(trails);
for i = 1:length(timestamps)
    trailData = trails.(timestamps{i});
    
    % Extract numeric values from the cell array
    lats = [lats; double(trailData{1})];       % Latitude
    lons = [lons; double(trailData{2})];       % Longitude
    % alts = [alts; double(trailData{3})]; % Height (altitude)
    
    % Remove the 'x' prefix from the timestamp
    timestampStr = timestamps{i};
    if timestampStr(1) == 'x'
        timestampStr = timestampStr(2:end); % Remove the first character ('x')
    end
    
    % Convert timestamp to numeric value
    timestamp = str2double(timestampStr);
    if isnan(timestamp)
        error('Failed to convert timestamp to numeric value: %s', timestamps{i});
    end
    times = [times; timestamp]; % Timestamp in milliseconds
end

% Convert timestamps to seconds (relative to the first timestamp)
times = (times - times(1)) / 1000; % Relative time in seconds

% Ensure lats, lons, and heights are numeric arrays (double)
lats = double(lats);
lons = double(lons);
% alts = double(alts);

% Step 3: Perform interpolation
% Create a finer time grid for interpolation

interpTimes = linspace(min(times), max(times), 1000); % 1000 points for smooth interpolation

% Interpolate latitude and longitude
interpLats = interp1(times, lats, interpTimes, 'spline'); % Spline interpolation
interpLons = interp1(times, lons, interpTimes, 'spline'); % Spline interpolation


figure;
subplot(2, 1, 1);
plot(times, lats, 'bo', 'DisplayName', 'Original Lats'); % Original latitude
hold on;
plot(interpTimes, interpLats, 'r-', 'DisplayName', 'Interpolated Lats'); % Interpolated latitude
xlabel('Time');
ylabel('Latitude');
title('Latitude Interpolation');
legend;
grid on;

subplot(2, 1, 2);
plot(times, lons, 'go', 'DisplayName', 'Original Lons'); % Original longitude
hold on;
plot(interpTimes, interpLons, 'm-', 'DisplayName', 'Interpolated Lons'); % Interpolated longitude
xlabel('Time');
ylabel('Longitude');
title('Longitude Interpolation');
legend;
grid on;

% Step 5: Plot the flight path on a geographic map
figure;
geoplot(lats, lons, 'bo-', 'DisplayName', 'Original Path'); % Original flight path
hold on;
geoplot(interpLats, interpLons, 'r-', 'DisplayName', 'Interpolated Path'); % Interpolated flight path
geobasemap('satellite'); % Use a street map
basemapName = "openstreetmap";
url = "a.tile.openstreetmap.org/${z}/${x}/${y}.png"; 
copyright = char(uint8(169));
attribution = copyright + "OpenStreetMap contributors";
addCustomBasemap(basemapName,url,"Attribution",attribution)
geobasemap('openstreetmap');
title('Flight Path on Geographic Map');
legend;
grid on;

end