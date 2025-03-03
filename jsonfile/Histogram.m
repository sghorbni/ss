function Histogram(filename)
%%%% COMPUTE HISTOGRAM
% reduce or aggregation path
% %%%%%%%%%%%%%%%

% Load JSON file
data = jsondecode(fileread(filename));

% Extract trails data
trails = data.trails;

% Initialize arrays for lats, lons, heights, and times
lats = [];
lons = [];
alts = [];
times = [];

% Loop through each timestamp in trails
timestamps = fieldnames(trails);
for i = 1:length(timestamps)
    trailData = trails.(timestamps{i});
    
    % Extract numeric values from the cell array
    lats = [lats; trailData{1}];       % Latitude
    lons = [lons; trailData{2}];       % Longitude
    alts = [alts; trailData{3}]; % Height (altitude)
    
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
alts = double(alts);

R = 6371; % Earth's radius in km

distances = zeros(length(lats) - 1, 1);
for i = 1:length(lats) - 1
    dLat = deg2rad(lats(i+1) - lats(i));
    dLon = deg2rad(lons(i+1) - lons(i));
    a = sin(dLat/2)^2 + cos(deg2rad(lats(i))) * cos(deg2rad(lats(i+1))) * sin(dLon/2)^2;
    c = 2 * atan2(sqrt(a), sqrt(1-a));
    distances(i) = R * c; % Distance in km
end

% Compute time differences
original_interpoint_time = diff(times); % Time differences in seconds

% Create a new time vector with 10-second intervals
start_time = times(1);
end_time = times(end);
%%%
%%%%
new_time = (start_time:15:end_time)'; % 10-second intervals (numeric)

% Ensure new_time is within the range of times
new_time = new_time(new_time >= start_time & new_time <= end_time);

% Interpolate latitude and longitude for the new time vector
lat_interp = interp1(times, lats, new_time, 'line');%spline
lon_interp = interp1(times, lons, new_time, 'line');%spline


% Step 4: Create a new structure for resampled data
resampled_data = struct();
resampled_data.times = new_time; % Keep as numeric (seconds)
resampled_data.lats = lat_interp;
resampled_data.lons = lon_interp;

% Display the resampled data
disp(resampled_data);

% Step 5: Save the resampled data to a new JSON file (optional)
resampled_filename = 'resampled_1996545206.json';
fid = fopen(resampled_filename, 'w');
fprintf(fid, jsonencode(resampled_data));
fclose(fid);

% Step 6: Calculate interpoint time for original and resampled data
original_interpoint_time = diff(times); % Original interpoint time
resampled_interpoint_time = diff(new_time); % Resampled interpoint time


% Step 7: Calculate course and speed every 6 points (1 minute apart)
% Initialize arrays for course and speed
course = [];
speed = [];

% Loop through the resampled data, skipping 6 points each time
for i = 1:6:length(lat_interp)-6
    % Get the current and next point (6 points apart)
    lat1 = lat_interp(i);
    lon1 = lon_interp(i);
    lat2 = lat_interp(i+6);
    lon2 = lon_interp(i+6);
    
    % Calculate the course (bearing) between the two points
    dLon = lon2 - lon1;
    y = sind(dLon) * cosd(lat2);
    x = cosd(lat1) * sind(lat2) - sind(lat1) * cosd(lat2) * cosd(dLon);
    bearing = atan2d(y, x);
    bearing = mod(bearing + 360, 360); % Ensure bearing is between 0 and 360 degrees
    course = [course; bearing];
    
    % Calculate the distance between the two points using Haversine formula
    R = 6371; % Earth's radius in kilometers
    dLat = deg2rad(lat2 - lat1);
    dLon = deg2rad(lon2 - lon1);
    a = sin(dLat/2) * sin(dLat/2) + cos(deg2rad(lat1)) * cos(deg2rad(lat2)) * sin(dLon/2) * sin(dLon/2);
    c = 2 * atan2(sqrt(a), sqrt(1-a));
    distance = R * c; % Distance in kilometers
    
    % Calculate the time difference (6 points * 10 seconds = 60 seconds = 1 minute)
    time_diff = 60; % 1 minute in seconds
    
    % Calculate speed in km/h
    speed_kmh = (distance / time_diff) * 3600; % Convert km/s to km/h
    speed = [speed; speed_kmh];
end

% Display the course and speed
disp('Course (degrees):');
disp(course);
disp('Speed (km/h):');
disp(speed);

%%%%%%%%%%%%%%%%%%%%%
% Step 9: Calculate the next point based on course, speed, and initial point
% Initialize arrays for calculated points
calculated_lats = [];
calculated_lons = [];

% Start with the initial point (first point in the resampled data)
initial_lat = lat_interp(1);
initial_lon = lon_interp(1);
calculated_lats = [calculated_lats; initial_lat];
calculated_lons = [calculated_lons; initial_lon];

% Loop through the course and speed vectors to calculate the next points
for i = 1:length(course)
    % Get the current course and speed
    current_course = course(i);
    current_speed = speed(i); % Speed in km/h
    
    % Convert speed to km/s (since time difference is 60 seconds)
    speed_kms = current_speed / 3600; % Convert km/h to km/s
    
    % Calculate the distance traveled in 60 seconds
    distance = speed_kms * 60; % Distance in kilometers
    
    % Calculate the new latitude and longitude using the course and distance
    R = 6371; % Earth's radius in kilometers
    delta = distance / R; % Angular distance in radians
    
    % Convert current latitude and longitude to radians
    lat1 = deg2rad(calculated_lats(end));
    lon1 = deg2rad(calculated_lons(end));
    
    % Calculate the new latitude and longitude
    lat2 = asin(sin(lat1) * cos(delta) + cos(lat1)...
        * sin(delta) * cos(deg2rad(current_course)));
    lon2 = lon1 + atan2(sin(deg2rad(current_course))...
        * sin(delta) * cos(lat1), cos(delta) - sin(lat1) * sin(lat2));
    
    % Convert the new latitude and longitude back to degrees
    lat2 = rad2deg(lat2);
    lon2 = rad2deg(lon2);
    
    % Append the calculated point to the arrays
    calculated_lats = [calculated_lats; lat2];
    calculated_lons = [calculated_lons; lon2];
end

% Step 10: Compare the calculated points with the original points
% Extract the original points corresponding to the calculated points
original_lats = lat_interp(1:6:end);
original_lons = lon_interp(1:6:end);

% Calculate the difference between calculated and original points
lat_diff = calculated_lats - original_lats;
lon_diff = calculated_lons - original_lons;

% Step 12: Reduce the number of aggregation points (optional)
% For example, reduce the number of points by half
reduced_lats = lat_interp(1:12:end);
reduced_lons = lon_interp(1:12:end);

% Step 2: Compute interpoint distance for original and resampled data
% Original interpoint distance
original_interpoint_distance = zeros(length(lats) - 1, 1);
for i = 1:length(lats) - 1
    dLat = deg2rad(lats(i+1) - lats(i));
    dLon = deg2rad(lons(i+1) - lons(i));
    a = sin(dLat/2)^2 + cos(deg2rad(lats(i))) * cos(deg2rad(lats(i+1))) * sin(dLon/2)^2;
    c = 2 * atan2(sqrt(a), sqrt(1-a));
    original_interpoint_distance(i) = R * c; % Distance in km
end

% Resampled interpoint distance
resampled_interpoint_distance = zeros(length(lat_interp) - 1, 1);
for i = 1:length(lat_interp) - 1
    dLat = deg2rad(lat_interp(i+1) - lat_interp(i));
    dLon = deg2rad(lon_interp(i+1) - lon_interp(i));
    a = sin(dLat/2)^2 + cos(deg2rad(lat_interp(i))) * cos(deg2rad(lat_interp(i+1))) * sin(dLon/2)^2;
    c = 2 * atan2(sqrt(a), sqrt(1-a));
    resampled_interpoint_distance(i) = R * c; % Distance in km
end
% Step 3: Plot interpoint time histograms (figure 1)
figure(1);
subplot(2,1,1);
histogram(original_interpoint_time, 'BinWidth', 1, 'FaceColor', 'r');
xlabel('Interpoint Time (seconds)');
ylabel('Frequency');
title('Original Data Interpoint Time Histogram');
grid on;

subplot(2,1,2);
histogram(resampled_interpoint_time, 'BinWidth', 1, 'FaceColor', 'b');
xlabel('Interpoint Time (seconds)');
ylabel('Frequency');
title('Resampled Data Interpoint Time Histogram');
grid on;

% Step 4: Plot interpoint distance histograms (figure 2)
figure(2);
subplot(2,1,1);
histogram(original_interpoint_distance, 'BinWidth', 0.1, 'FaceColor', 'r');
xlabel('Interpoint Distance (km)');
ylabel('Frequency');
title('Original Data Interpoint Distance Histogram');
grid on;

subplot(2,1,2);
histogram(resampled_interpoint_distance, 'BinWidth', 0.1, 'FaceColor', 'b');
xlabel('Interpoint Distance (km)');
ylabel('Frequency');
title('Resampled Data Interpoint Distance Histogram');
grid on;




end