function [lats,lons,alts]=jsnpolar1(filename)
data = jsondecode(fileread(filename));

% Extract trails data
trails = data.trails;
timestamps = fieldnames(trails); % Get all timestamps (keys)

% Initialize arrays to store latitude, longitude, and altitude
lats = [];
lons = [];
alts = [];

% Loop through each timestamp and extract the data
for i = 1:numel(timestamps)
    trailData = trails.(timestamps{i}); % Get the data for the current timestamp
    % Ensure the trailData has at least 3 elements (lat, lon, alt)
    if numel(trailData) >= 3
        lats = [lats; trailData{1}]; % Latitude is the first value
        lons = [lons; trailData{2}]; % Longitude is the second value
        alts = [alts; trailData{3}]; % Altitude is the third value
    else
        warning('Skipping timestamp %s: insufficient data.', timestamps{i});
    end
end

% Convert altitude from feet to meters (optional, depending on your use case)
alts = alts * 0.3048; % 1 foot = 0.3048 meters

% Step 2: Convert geographic coordinates to Cartesian coordinates
R = 6371000; % Earth's radius in meters
lat_rad = deg2rad(lats); % Convert latitude to radians
lon_rad = deg2rad(lons); % Convert longitude to radians

% Convert to Cartesian coordinates (x, y, z)
x = (R + alts) .* cos(lat_rad) .* cos(lon_rad);
y = (R + alts) .* cos(lat_rad) .* sin(lon_rad);
z = (R + alts) .* sin(lat_rad);

% Step 3: Convert Cartesian coordinates to polar coordinates
[azimuth, elevation, radius] = cart2sph(x, y, z); % Convert to polar coordinates

% Step 4: Apply rotation to azimuth angle
rotation_angle = deg2rad(25); % Add 10 degrees to the azimuth angle (convert to radians)
azimuth_rotated = azimuth + rotation_angle;

% Step 5: Convert polar coordinates back to Cartesian coordinates
[x_rotated, y_rotated, z_rotated] = sph2cart(azimuth_rotated, elevation, radius);

% Step 6: Convert Cartesian coordinates back to geographic coordinates
lat_rotated = rad2deg(asin(z_rotated ./ (R + alts))); % New latitude
lon_rotated = rad2deg(atan2(y_rotated, x_rotated)); % New longitude

% Step 7: Visualize the results
% Plot original and rotated points in polar coordinates
figure(1);
polarplot(azimuth, radius, 'b-', 'LineWidth', 2); % Original points
hold on;
polarplot(azimuth_rotated, radius, 'r-', 'LineWidth', 2); % Rotated points
title('Polar Plot of Original and Rotated Points');
legend('Original', 'Rotated');
% exportgraphics(gcf, 'polar_plot.png', 'Resolution', 300); % Save as PNG with 300 DPI
% Plot original and rotated points on a geographic map
figure(2);
geoplot(lats, lons, 'b-', 'LineWidth', 2); % Original points
hold on;
geoplot(lat_rotated, lon_rotated, 'r-', 'LineWidth', 2); % Rotated points
geobasemap('satellite'); % Use a street map
basemapName = "openstreetmap";
url = "a.tile.openstreetmap.org/${z}/${x}/${y}.png"; 
copyright = char(uint8(169));
attribution = copyright + "OpenStreetMap contributors";
addCustomBasemap(basemapName,url,"Attribution",attribution)
geobasemap('openstreetmap');
title('Geographic Map of Flight Path');
title('Geographic Map of Flight Path');
legend('Original Path', 'Rotated Path');
grid on;
exportgraphics(gcf, 'geographic_map.png', 'Resolution', 300); % Save as PNG with 300 DPI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end