function out = lab5PressureDistribution(scanivalve)

% Pressure Distribution Calculations

% Will need the matrix contaning the group data
% Allocation
%[pressure,temperature,density,airspeed,pitotDynamic,auxDynamic,...
    %scanivalve,angle,xpp_scaled,ypp_scaled] = readData();
[row, col, lay] = size(scanivalve);

p11 = zeros(3,1,lay);

x1 = [2.8 2.1];
x2 = [2.1 2.8];
X1 = [ones(1,length(x1)), x1];
X2 = [ones(1,length(x2)), x2];

for i = 1:row
    for k = 1:lay
        y1 = scanivalve(i,8:9,k);
        coef1 = polyfit(x1,y1,1);
        Y1 = coef1(1)*3.5 + coef1(2);
        
        y2 = scanivalve(i,10:11,k);
        coef2 = polyfit(x2,y2,1);
        Y2 = coef2(1)*3.5 + coef2(2);
        
        p11(i,k) = mean([Y1 Y2]);
        
    end
end

out = p11;

end