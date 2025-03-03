function azimuth = calculateAzimuth(lat1, lon1, lat2, lon2)
    % Convert degrees to radians
    lat1 = deg2rad(lat1);
    lon1 = deg2rad(lon1);
    lat2 = deg2rad(lat2);
    lon2 = deg2rad(lon2);

    % Compute azimuth
    dLon = lon2 - lon1;
    y = sin(dLon) * cos(lat2);
    x = cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(dLon);
    azimuth = rad2deg(atan2(y, x));

    % Ensure azimuth is positive (0 to 360 degrees)
    azimuth = mod(azimuth + 360, 360);
end