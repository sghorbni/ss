%%%% read lats lons and extract them
function [lats,lons]=extractjsondata2(filename)
data = jsondecode(fileread(filename));
% Extract trails data
trails = data.trails;

keys = fieldnames(trails); % Get all timestamps
lats = [];
lons = [];
% Loop through each timestamp and extract lat/lon
for i = 1:length(keys)
    trailData = trails.(keys{i});% Get data for this timestamp
    
    lat = trailData(1); % Latitude
    lon = trailData(2);
   
    if isnumeric(lat) && isnumeric(lon)
        lats = [lats; lat]; % Append latitude
        lons = [lons; lon]; % Append longitude
    else
       lats = [lats; cell2mat(lat)]; % Convert to numeric
        lons = [lons; cell2mat(lon)]; % Convert to numeric
    end
end
disp('Latitudes (lats):');
disp(lats);
disp('Longitudes (lons):');
disp(lons);

figure;
geoplot(lats, lons, '-o', 'LineWidth', 2, 'MarkerSize', 4); % Plot the flight path
geobasemap('satellite');

basemapName = "openstreetmap";
url = "a.tile.openstreetmap.org/${z}/${x}/${y}.png"; 
copyright = char(uint8(169));
attribution = copyright + "OpenStreetMap contributors";
addCustomBasemap(basemapName,url,"Attribution",attribution)

geobasemap('openstreetmap');

% Add title and labels
title('Flight Path Visualization');
%xlabel('Longitude');
%ylabel('Latitude');

% Set map limits (optional)
geolimits([min(lats) max(lats)], [min(lons) max(lons)]);
% Save the map as an image
saveas(gcf, 'flight_path.png');


end