function [lats,lons,alts,times]=jsondistnce2(filename)

% Load JSON file
data = jsondecode(fileread('1996202916.json'));

% Extract trails data
trails = data.trails;

lats = [];
lons = [];
alts = [];
times = [];

% Loop through each timestamp in trails
timestamps = fieldnames(trails);
for i = 1:length(timestamps)
    trailData = trails.(timestamps{i});
    
    % Extract numeric values from the cell array
    lats = [lats; double(trailData{1})];       % Latitude
    lons = [lons; double(trailData{2})];       % Longitude
    alts = [alts; double(trailData{3})]; % Height (altitude)
    
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

% Initialize arrays for distances, speed, last_speeds, and speed_errors
distances = zeros(length(lats) - 1, 1);
speed = zeros(length(lats) - 1, 1);
last_speeds = zeros(length(lats) - 1, 1);
speed_errors = zeros(length(lats) - 1, 1);

% Compute distances, speed, last_speeds, and speed_errors in a loop
R = 6371; % Earth's radius in km
for i = 1:length(lats) - 1
    % Compute distance
    dLat = deg2rad(lats(i+1) - lats(i));
    dLon = deg2rad(lons(i+1) - lons(i));
    a = sin(dLat/2)^2 + cos(deg2rad(lats(i))) * cos(deg2rad(lats(i+1))) * sin(dLon/2)^2;
    c = 2 * atan2(sqrt(a), sqrt(1-a));
    distances(i) = R * c; % Distance in km
    
    % Compute speed (km/h)
    timeDifference = times(i+1) - times(i); % Time difference in seconds
    speed(i) = distances(i) / timeDifference * 3600; % Convert to km/h
    
    % Compute last_speeds and speed_errors
    if i > 1
        last_speeds(i) = speed(i-1);
        speed_errors(i) = abs(speed(i) - last_speeds(i));
    else
        last_speeds(i) = NaN; % No previous speed for the first point
        speed_errors(i) = NaN; % No error for the first point
    end
end

% Smooth data for Figure 3
smoothed_speed = smoothdata(speed, 'movmean', 5); % Smooth speed with a moving average of window size 5
smoothed_last_speeds = smoothdata(last_speeds, 'movmean', 5); % Smooth last speeds
smoothed_speed_errors = smoothdata(speed_errors, 'movmean', 5); % Smooth speed errors

% Plot 1: Flight path on a geographic map
figure(1);
geoplot(lats, lons, '-o', 'LineWidth', 2, 'MarkerSize', 4);
geobasemap('satellite');
basemapName = "openstreetmap";
url = "a.tile.openstreetmap.org/${z}/${x}/${y}.png"; 
copyright = char(uint8(169));
attribution = copyright + "OpenStreetMap contributors";
addCustomBasemap(basemapName,url,"Attribution",attribution)
geobasemap('openstreetmap');
legend('Transformed Path', 'Original Path');
title('Flight Path on Geographic Map');
legend('Flight Path', 'Start', 'End');
grid on;
% Save the figure as PNG using exportgraphics
exportgraphics(gcf, 'FlightPath_GeographicMap.png', 'Resolution', 300);

% Plot 2: Flight path in (x, y) coordinates
figure(2);
plot(lons, lats, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'Color', 'b');
xlabel('Longitude');
ylabel('Latitude');
title('Flight Path in (X, Y) Coordinates');
grid on;
hold on;
plot(lons(1), lats(1), 'go', 'MarkerSize', 10, 'LineWidth', 2); % Start point
plot(lons(end), lats(end), 'ro', 'MarkerSize', 10, 'LineWidth', 2); % End point
legend('Flight Path', 'Start', 'End');
hold off;
% Save the figure as PNG using exportgraphics
exportgraphics(gcf, 'FlightPath_XYCoordinates.png', 'Resolution', 300);

% Plot 3: Speed vs. Time with comparison to last speeds and speed errors (smoothed)
%%%%%%%%%%%%
figure(3);
subplot(3, 1, 1); % First subplot: Current Speed
plot(times(2:end), speed, '-o', 'DisplayName', 'Current Speed');
hold on
plot(times(2:end), smoothed_speed, '-o', 'DisplayName', 'Smoothed Current Speed');
xlabel('Time (seconds)');
ylabel('Speed (km/h)');
title('Current Speed vs. Time');
legend;
grid on;

subplot(3, 1, 2); % Second subplot: Last Speed
plot(times(2:end), last_speeds, '-o', 'DisplayName', 'Last Speed');
hold on
plot(times(2:end), smoothed_last_speeds, '-o', 'DisplayName', 'Smoothed Last Speed');
xlabel('Time (seconds)');
ylabel('Speed (km/h)');
title('Last Speed vs. Time');
legend;
grid on;

subplot(3, 1, 3); % Third subplot: Speed Errors
plot(times(2:end), speed_errors, '-o', 'DisplayName', 'Speed Errors');
hold on
plot(times(2:end), smoothed_speed_errors, '-o', 'DisplayName', 'Smoothed Speed Errors');
xlabel('Time (seconds)');
ylabel('Speed Error (km/h)');
title('Speed Errors vs. Time');
legend;
grid on;
exportgraphics(gcf, 'Smoothed_Speed_vs_Time_and_Errors.png', 'Resolution', 300);

%%%%%%%%%%%
% figure(3);
% plot(times(2:end), smoothed_speed, '-o', 'DisplayName', 'Smoothed Current Speed');
% hold on;
% plot(times(2:end), smoothed_last_speeds, '-o', 'DisplayName', 'Smoothed Last Speed');
% plot(times(2:end), smoothed_speed_errors, '-o', 'DisplayName', 'Smoothed Speed Errors');
% xlabel('Time (seconds)');
% ylabel('Speed (km/h) / Speed Error (km/h)');
% title('Smoothed Speed vs. Time (Current, Last, and Errors)');
% legend;
% grid on;
% hold off;
% % Save the figure as PNG using exportgraphics
% exportgraphics(gcf, 'Smoothed_Speed_vs_Time_and_Errors.png', 'Resolution', 300);

% Plot 4: Height vs. Time
figure(4);
plot(times, alts, '-o');
xlabel('Time (seconds)');
ylabel('Height (feet)');
title('Height vs. Time');
grid on;
% Save the figure as PNG using exportgraphics
exportgraphics(gcf, 'Height_vs_Time.png', 'Resolution', 300);

% Plot 5: Filtered and smoothed flight path on a geographic map
minTimeDifference = 10; % Minimum time difference in seconds
timeDiffs = diff(times);
validIndices = [true; timeDiffs >= minTimeDifference]; % Keep the first point and points with time difference >= 10 seconds

% Ensure validIndices has the correct length
if length(validIndices) ~= length(times)
    validIndices = [true; timeDiffs >= minTimeDifference];
end

% Filtered data
filteredLats = lats(validIndices);
filteredLons = lons(validIndices);
filteredTimes = times(validIndices);
filteredHeights = alts(validIndices);

% Smooth filtered data
smoothed_filteredLats = smoothdata(filteredLats, 'movmean', 5); % Smooth latitudes
smoothed_filteredLons = smoothdata(filteredLons, 'movmean', 5); % Smooth longitudes

% Plot smoothed and filtered flight path on a geographic map
figure(5);
geoplot(smoothed_filteredLats, smoothed_filteredLons, '-o', 'LineWidth', 2, 'MarkerSize', 4);
geobasemap('satellite');
basemapName = "openstreetmap";
url = "a.tile.openstreetmap.org/${z}/${x}/${y}.png"; 
copyright = char(uint8(169));
attribution = copyright + "OpenStreetMap contributors";
addCustomBasemap(basemapName,url,"Attribution",attribution)
geobasemap('openstreetmap');
legend('Transformed Path', 'Original Path');
title('Smoothed and Filtered Flight Path on Geographic Map');
legend('Smoothed Flight Path', 'Start', 'End');
grid on;
% Save the figure as PNG using exportgraphics
exportgraphics(gcf, 'Smoothed_FilteredFlightPath_GeographicMap.png', 'Resolution', 300);


end