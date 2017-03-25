function [pressure,temperature,density,airspeed,pitotDynamic,auxDynamic,...
    scanivalve,angle,xpp_scaled,ypp_scaled] = readData()
% Requirements
%   - Each Group's data must be 180x28
%   - Filename: 'AirfoilPressure_S013_GXX.csv'
%               XX = Group Number

% Allocation
data = zeros(180,28,11);

% Read data by group (CHANGE GROUP INDICES FOR NUMBER OF FILES)
for group = 1:11
    filename = sprintf('AirfoilPressure_S013_G%02d.csv',group);
    data(:,:,group) = csvread(filename,1,0);
end

% Take averages at each velocity for each angle of attack
avgData = zeros(9,28,11);
for i = 1:9
    for j = 1:28
        for k = 1:11
            avgData(i,j,k) = mean(data((i-1)*20+1:i*20,j,k));
        end
    end
end

% Break out into separate arrays
pressure = avgData(:,1,:);
temperature = avgData(:,2,:);
density = avgData(:,3,:);
airspeed = avgData(:,4,:);
pitotDynamic = avgData(:,5,:);
auxDynamic = avgData(:,6,:);
scanivalve = avgData(:,7:22,:);
angle = avgData(:,23,:);

% Add pressure port locations
xpp = [0, 5, 10, 20, 30, 40, 50, 60, 80, 100, 80, 60, 40, 30, 20, 10, 5];
ypp = [4.19, 9.45, 11.48, 13.6, 14, 13.64, 12.58, 10.95, 6.25, 0, 0, 0, 0, 0, 0.04, 0.5, 1.11];

% Scale the port locations for the chord length given
c = 3.5; %Chord in Inches
xpp_scaled = c*xpp/100;
ypp_scaled = c*ypp/100;

end