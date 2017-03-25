% Pressure Distribution Calculations

% Will need the matrix contaning the group data
% Allocation
[pressure,temperature,density,airspeed,pitotDynamic,auxDynamic,...
    scanivalve,angle,xpp_scaled,ypp_scaled] = readData();
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


pressure(:,1:10,:) = scanivalve(:,1:10,:);
pressure(:,11,:) = p11;
pressure(:,12:17,:) = scanivalve(:,11:16,:);
[~, nOfPorts, ~] = size(pressure);

c = 3.5; %cord length in inches
Cp = zeros(9,17,lay);
for i = 1:9
    for j = 1:17
        for k = 1:lay
            Cp(i,j,k) = pressure(i,j,k)/pitotDynamic(i,1,k);
        end
    end
end

n = zeros(row,1,lay);
a = zeros(row,1,lay);
Cn = zeros(row,1,lay);
Ca = zeros(row,1,lay);

for i = 1:row
    for j = 1:(nOfPorts - 1)
        for k = 1:lay
            n(i,k) = n(i,k) - 0.5*(pressure(i,j,k) + pressure(i,j+1,k))...
                *(xpp_scaled(1,j+1) - xpp_scaled(1,j));
            a(i,k) = a(i,k) + 0.5*(pressure(i,j,k) + pressure(i,j+1,k))...
                *(ypp_scaled(1,j+1) - ypp_scaled(1,j));
            
            Cn(i,k) = Cn(i,k) - 0.5*(Cp(i,j,k) + Cp(i,j+1,k))...
                *(xpp_scaled(1,j+1) - xpp_scaled(1,j))/3.5;
            Ca(i,k) = Ca(i,k) + 0.5*(Cp(i,j,k) + Cp(i,j+1,k))...
                *(ypp_scaled(1,j+1) - ypp_scaled(1,j))/3.5;
        end
    end
end

coefLift = Cn.*cosd(angle) - Ca.*sind(angle);
coefDrag = Cn.*sind(angle) + Ca.*cosd(angle);

cL10Layered = coefLift([1 4 7],:,:);
cD10Layered = coefDrag([1 4 7],:,:);
cL20Layered = coefLift([2 5 8],:,:);
cD20Layered = coefDrag([2 5 8],:,:);
cL30Layered = coefLift([3 6 9],:,:);
cD30Layered = coefDrag([3 6 9],:,:);
legends = {'10 m/s', '20 m/s', '30 m/s'};

for k = 1:lay
    cL10(:,k) = cL10Layered(:,:,k);
    cD10(:,k) = cD10Layered(:,:,k);
    angle10(:,k) = angle([1 4 7],1,k);
    
    cL20(:,k) = cL20Layered(:,:,k);
    cD20(:,k) = cD20Layered(:,:,k);
    angle20(:,k) = angle([2 5 8],1,k);
    
    
    cL30(:,k) = cL30Layered(:,:,k);
    cD30(:,k) = cD30Layered(:,:,k);
    angle30(:,k) = angle([3 6 9],1,k);
end



figure
hold on
p1 = plot(angle10,cL10,'bo');
xlabel('Angle of Attack')
ylabel('Coefficient of Lift')
title('Angle of Attack vs. Coefficient of Lift at 10 m/s')

% figure
p2 = plot(angle20,cL20,'ro');
xlabel('Angle of Attack')
ylabel('Coefficient of Lift')
title('Angle of Attack vs. Coefficient of Lift at 20 m/s')

% figure
p3 = plot(angle30,cL30,'ko');
xlabel('Angle of Attack')
ylabel('Lift Coefficient')
title('Angle of Attack vs. Coefficient of Lift at 30 m/s')
legend([p1(1,1) p2(1,1) p3(1,1)], '10 m/s', '20 m/s', '30 m/s')



figure
hold on
p4 = plot(angle10,cD10,'bo');
% xlabel('Angle of Attack')
% ylabel('Coefficient of Drag')
% title('Angle of Attack vs. Coefficient of Drag at 10 m/s')

% figure
p5 = plot(angle20,cD20,'ro');
% xlabel('Angle of Attack')
% ylabel('Coefficient of Drag')
% title('Angle of Attack vs. Coefficient of Drag at 20 m/s')

% figure
p6 = plot(angle30,cD30,'ko');
xlabel('Angle of Attack')
ylabel('Drag Coefficient')
title('Angle of Attack vs. Coefficient of Drag at 30 m/s')
legend([p4(1,1) p5(1,1) p6(1,1)], '10 m/s', '20 m/s', '30 m/s')



figure
hold on 
p7 = plot(angle10,cL10./cD10,'bo');
% xlabel('Angle of Attack')
% ylabel('L/D')
% title('Angle of Attack vs. L/D at 10 m/s')

% figure
p8 = plot(angle20,cL20./cD20,'ro');
% xlabel('Angle of Attack')
% ylabel('L/D')
% title('Angle of Attack vs. L/D at 20 m/s')

% figure
p9 = plot(angle30,cL30./cD30,'ko');
xlabel('Angle of Attack')
ylabel('L/D Ratio')
title('Angle of Attack vs. L/D at 30 m/s')
legend([p7(1,1) p8(1,1) p9(1,1)], '10 m/s', '20 m/s', '30 m/s')





