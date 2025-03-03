function resample(filename)

%%%%%%%%%%%%%%%%%%%%
%resample data by intrpltion
%the first step in 10 sec%
%%%%%%
% Load JSON file
data = jsondecode(fileread(filename));

% Extract trails data
trails = data.trails;

% Initialize arrays for lats, lons, heights, and times
lats = [];
lons = [];
times = [];

% Loop through each timestamp in trails
timestamps = fieldnames(trails);
for i = 1:length(timestamps)
    trailData = trails.(timestamps{i});
    
    % Extract numeric values from the cell array
    lats = [lats; trailData{1}];       % Latitude
    lons = [lons; trailData{2}];       % Longitude
    
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

% Create a new time vector with 10-second intervals
start_time = times(1);
end_time = times(end);
new_time = (start_time:10:end_time)'; % 10-second intervals (numeric)

% Ensure new_time is within the range of times
new_time = new_time(new_time >= start_time & new_time <= end_time);

% Interpolate latitude and longitude for the new time vector
lat_interp = interp1(times, lats, new_time, 'line'); %spline
lon_interp = interp1(times, lons, new_time, 'line');%spline

% Step 4: Create a new structure for resampled data
resampled_data = struct();
resampled_data.times = new_time; % Keep as numeric (seconds)
resampled_data.lats = lat_interp;
resampled_data.lons = lon_interp;

% Display the resampled data
disp(resampled_data);

% Step 5: Save the resampled data to a new JSON file (optional)
resampled_filename = 'resampled_1995899608.json';
fid = fopen(resampled_filename, 'w');
fprintf(fid, jsonencode(resampled_data));
fclose(fid);

% Step 6: Plot the original and resampled data
figure;
geoplot(lats, lons, 'r-', 'LineWidth', 2); % Original data
hold on;
geoplot(lat_interp, lon_interp, 'b--', 'LineWidth', 2); % Resampled data
geobasemap('satellite');

basemapName = "openstreetmap";
url = "a.tile.openstreetmap.org/${z}/${x}/${y}.png"; 
copyright = char(uint8(169));
attribution = copyright + "OpenStreetMap contributors";
addCustomBasemap(basemapName,url,"Attribution",attribution)

geobasemap('openstreetmap');
legend('Original', 'Resampled (10-second intervals)');
title('Resampled Flight Path');


end