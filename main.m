close all;
% Read Data
[pressure,temperature,density,airspeed,pitotDynamic,auxDynamic,...
    scanivalve,angle,xpp_scaled,ypp_scaled] = readData();

% Find port 11
p11 = lab5PressureDistribution(scanivalve);

% Add p11 to port data
ports = zeros(9,17,11);
ports(:,1:16,:) = scanivalve;
ports(:,17,:) = p11;

% Plot Pressure Distribution
Cp = pressureDistribution(ports,pitotDynamic,xpp_scaled,angle);