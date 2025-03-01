function geojsonintrpltion(filename)

data = jsondecode(fileread(filename));
fid = fopen(filename);

if fid == -1
    error('File could not be opened.');
end
raw = fread(fid, inf, 'uint8=>char');
fclose(fid);
data = jsondecode(raw');

% Initialize arrays to store lon and lat values
lons = [];
lats = [];

% Check if the data contains 'features'
if isfield(data, 'features')
    features = data.features;
    
    % Loop through each feature
    for i = 1:length(features)
        feature = features(i);
        
        % Check if the feature contains 'geometry'
        if isfield(feature, 'geometry')
            geometry = feature.geometry;
            
            % Check the geometry type
            if isfield(geometry, 'type') && strcmpi(geometry.type, 'MultiLineString')
                % Extract coordinates for MultiLineString
                coordinates = geometry.coordinates;
                
                % Loop through each line in the MultiLineString
                for j = 1:length(coordinates)
                    line = coordinates(1, j, :); % Get the j-th line
                    lons = [lons; line(:, 1)]; % Extract longitudes
                    lats = [lats; line(:, 2)]; % Extract latitudes
                end
            end
        end
    end
end

% Interpolation
% Create a parameter for interpolation (e.g., distance along the path)
t = 1:length(lons);
tq = 1:0.1:length(lons); % Query points for interpolation

% Perform linear interpolation
lons_interp = interp1(t, lons, tq, 'linear');
lats_interp = interp1(t, lats, tq, 'linear');

% Plot the original and interpolated coordinates
figure;
plot(lons, lats, 'b-', 'LineWidth', 1.5); % Plot original data as blue lines
hold on;
plot(lons_interp, lats_interp, 'r--', 'LineWidth', 1.5); % Plot interpolated data as red dashed lines
xlabel('Longitude');
ylabel('Latitude');
title('GeoJSON Coordinates Plot with Interpolation');
legend('Original', 'Interpolated');
axis equal;
grid on;
hold off;

end