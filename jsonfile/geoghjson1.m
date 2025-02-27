function [lats,lons]=geoghjson1(filename)
data = jsondecode(fileread(filename));
fid = fopen(filename);
if fid == -1
    error('File could not be opened.');
end
raw = fread(fid, inf, 'uint8=>char');
fclose(fid);
data = jsondecode(raw');
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

% Plot the coordinates
figure;
geoplot(lats, lons, 'b-', 'LineWidth', 1.5); % Plot as blue lines

basemapName = "openstreetmap";
url = "a.tile.openstreetmap.org/${z}/${x}/${y}.png"; 
copyright = char(uint8(169));
attribution = copyright + "OpenStreetMap contributors";
addCustomBasemap(basemapName,url,"Attribution",attribution)

geobasemap('openstreetmap');


title('GeoJSON Coordinates Plot on Geographic Map');
grid on;
end