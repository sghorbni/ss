function  [lats,lons,alts,times,errors]=jsndistnce1(filename)

% Load JSON file
data = jsondecode(fileread('1996202916.json'));

% Extract trails data
trails = data.trails;

% Initialize arrays for lats, lons, heights, and times
lats = [];
lons = [];
alts = [];
times = [];
errors=[];

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

% Compute azimuth angles
azimuth = zeros(length(lats) - 1, 1);
for i = 1:length(lats) - 1
    azimuth(i) = calculateAzimuth(lats(i), lons(i), lats(i+1), lons(i+1));
end

% Display the first few values of azimuth to verify
disp('First few values of azimuth:');
disp(azimuth(1:10));

% Compute distances between consecutive points
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
timeDifferences = diff(times); % Time differences in seconds

% Compute speed (km/h)
speed = distances ./ timeDifferences * 3600; % Convert to km/h

% Compute last speeds (shifted by one point)
last_speeds = [NaN; speed(1:end-1)]; % Last speeds for comparison
for j=1:length(lats)-1
    errors=[errors;abs(speed(j)-last_speeds(j))];
end

figure(1);
geoplot(lats, lons, '-o', 'LineWidth', 2, 'MarkerSize', 4);
geobasemap('satellite');

basemapName = "openstreetmap";
url = "a.tile.openstreetmap.org/${z}/${x}/${y}.png"; 
copyright = char(uint8(169));
attribution = copyright + "OpenStreetMap contributors";
addCustomBasemap(basemapName,url,"Attribution",attribution)

geobasemap('openstreetmap');
title('Flight Path on Geographic Map');
legend('Flight Path', 'Start', 'End');
grid on;
saveas(gcf, 'FlightPath_GeographicMap.png');

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
hold off
saveas(gcf, 'FlightPath_XYCoordinates.png');
%%%%%%%%%%
% figure(3);
% plot(times(2:end), speed, '-o',times(2:end), last_speeds,...
%     '-o',times(2:end), errors, '-o')
%%%%%%%%
% %Plot 3: Speed vs. Time with comparison to last speeds
% figure(3);
% plot(times(2:end), speed, '-o', 'DisplayName', 'Current Speed');
% hold on;
% plot(times(2:end), last_speeds, '-o', 'DisplayName', 'Last Speed');
% hold on;
% plot(times(2:end), errors, '-o', 'DisplayName', 'errors');
% xlabel('Time (seconds)');
% ylabel('Speed (km/h)');
% title('Speed vs. Time (Current vs. Last Speed)');
% legend;
% grid on;
% hold off
% % % % % % % % % % % % 
figure(3);
subplot(3, 1, 1); % First subplot: Current Speed
plot(times(2:end), speed, '-o', 'DisplayName', 'Current Speed');
xlabel('Time (seconds)');
ylabel('Speed (km/h)');
title('Current Speed vs. Time');
legend;
grid on;

subplot(3, 1, 2); % Second subplot: Last Speed
plot(times(2:end), last_speeds, '-o', 'DisplayName', 'Last Speed');
xlabel('Time (seconds)');
ylabel('Speed (km/h)');
title('Last Speed vs. Time');
legend;
grid on;

subplot(3, 1, 3); % Third subplot: Speed Errors
plot(times(2:end), errors, '-o', 'DisplayName', 'Speed Errors');
xlabel('Time (seconds)');
ylabel('Speed Error (km/h)');
title('Speed Errors vs. Time');
legend;
grid on;
saveas(gcf, 'Speed_vs_Time_and_Errors.png');

% % % % % % % % % % 
% Plot 4: Height vs. Time
figure(4);
plot(times, alts, '-o');
xlabel('Time (seconds)');
ylabel('Height (feet)');
title('Height vs. Time');
grid on;

% Plot 5: Filtered points (remove noisy points) on a geographic map
% Filter condition: Remove points where the time difference is less than 10 seconds
minTimeDifference = 10; % Minimum time difference in seconds
validIndices = [true; timeDifferences >= minTimeDifference]; % Keep the first point and points with time difference >= 10 seconds

% Filtered data
filteredLats = lats(validIndices);
filteredLons = lons(validIndices);
filteredTimes = times(validIndices);
filteredHeights = alts(validIndices);
saveas(gcf, 'Height_vs_Time.png');
% Plot filtered flight path on a geographic map
figure(5);
geoplot(filteredLats, filteredLons, '-o', 'LineWidth', 2, 'MarkerSize', 4);
geobasemap('satellite');

basemapName = "openstreetmap";
url = "a.tile.openstreetmap.org/${z}/${x}/${y}.png"; 
copyright = char(uint8(169));
attribution = copyright + "OpenStreetMap contributors";
addCustomBasemap(basemapName,url,"Attribution",attribution)

geobasemap('openstreetmap');


title('Filtered Flight Path on Geographic Map (Noise Removed)');
legend('Filtered Flight Path', 'Start', 'End');
grid on;
saveas(gcf, 'FilteredFlightPath_GeographicMap.png');

end