close all;
clear all;
clc;

% Read in Data Files
VVS013G01 = DataRead('VelocityVoltage_S013_G01.csv'); % This data is not used in the calculations
VVS013G02 = DataRead('VelocityVoltage_S013_G02.csv');
VVS013G03 = DataRead('VelocityVoltage_S013_G03.csv');
VVS013G04 = DataRead('VelocityVoltage_S012_G04.csv'); % The data from Group 4 in Section 13 caused a negative k value so I took the data from section 12
VVS013G05 = DataRead('VelocityVoltage_S013_G05.csv');
VVS013G06 = DataRead('VelocityVoltage_S013_G06.csv');
VVS013G07 = DataRead('VelocityVoltage_S013_G07.csv');
VVS013G08 = DataRead('VelocityVoltage_S013_G08.csv');
VVS013G09 = DataRead('VelocityVoltage_S012_G09.csv'); % Group 9 in Section 13 had a corrupted data file so I took the data from section 12
VVS013G10 = DataRead('VelocityVoltage_S013_G10.csv');
VVS013G11 = DataRead('VelocityVoltage_S013_G11.csv'); % Using data from group 11 in place of group 1 due to the error in the data from group 1

% Calculate k-values and velocity for each Voltage

[V1_1, kV1_1, V2_1, kV2_1, V3_1, kV3_1, V4_1, kV4_1, V5_1, kV5_1, velocity1_1, velocity2_1, velocity3_1, velocity4_1, velocity5_1, PSv1_1, PSv2_1, PSv3_1, PSv4_1, PSv5_1, SRv1_1, SRv2_1, SRv3_1, SRv4_1, SRv5_1] = Calc_v(VVS013G11,VVS013G02);
[V1_2, kV1_2, V2_2, kV2_2, V3_2, kV3_2, V4_2, kV4_2, V5_2, kV5_2, velocity1_2, velocity2_2, velocity3_2, velocity4_2, velocity5_2, PSv1_2, PSv2_2, PSv3_2, PSv4_2, PSv5_2, SRv1_2, SRv2_2, SRv3_2, SRv4_2, SRv5_2] = Calc_v(VVS013G03,VVS013G04);
[V1_3, kV1_3, V2_3, kV2_3, V3_3, kV3_3, V4_3, kV4_3, V5_3, kV5_3, velocity1_3, velocity2_3, velocity3_3, velocity4_3, velocity5_3, PSv1_3, PSv2_3, PSv3_3, PSv4_3, PSv5_3, SRv1_3, SRv2_3, SRv3_3, SRv4_3, SRv5_3] = Calc_v(VVS013G05,VVS013G06);
[V1_4, kV1_4, V2_4, kV2_4, V3_4, kV3_4, V4_4, kV4_4, V5_4, kV5_4, velocity1_4, velocity2_4, velocity3_4, velocity4_4, velocity5_4, PSv1_4, PSv2_4, PSv3_4, PSv4_4, PSv5_4, SRv1_4, SRv2_4, SRv3_4, SRv4_4, SRv5_4] = Calc_v(VVS013G07,VVS013G08);
[V1_5, kV1_5, V2_5, kV2_5, V3_5, kV3_5, V4_5, kV4_5, V5_5, kV5_5, velocity1_5, velocity2_5, velocity3_5, velocity4_5, velocity5_5, PSv1_5, PSv2_5, PSv3_5, PSv4_5, PSv5_5, SRv1_5, SRv2_5, SRv3_5, SRv4_5, SRv5_5] = Calc_v(VVS013G09,VVS013G10);

% Put calculated values into vectors
V = [V1_1 V2_1 V3_1 V4_1 V5_1 V1_2 V2_2 V3_2 V4_2 V5_2 V1_3 V2_3 V3_3 V4_3 V5_3 V1_4 V2_4 V3_4 V4_4 V5_4 V1_5 V2_5 V3_5 V4_5 V5_5];
k = [kV1_1 kV2_1 kV3_1 kV4_1 kV5_1 kV1_2 kV2_2 kV3_2 kV4_2 kV5_2 kV1_3 kV2_3 kV3_3 kV4_3 kV5_3 kV1_4 kV2_4 kV3_4 kV4_4 kV5_4 kV1_5 kV2_5 kV3_5 kV4_5 kV5_5];
velocity = [velocity1_1, velocity2_1, velocity3_1, velocity4_1, velocity5_1, velocity1_2, velocity2_2, velocity3_2, velocity4_2, velocity5_2, velocity1_3, velocity2_3, velocity3_3, velocity4_3, velocity5_3, velocity1_4, velocity2_4, velocity3_4, velocity4_4, velocity5_4, velocity1_5, velocity2_5, velocity3_5, velocity4_5, velocity5_5];
PSv = [PSv1_1, PSv2_1, PSv3_1, PSv4_1, PSv5_1, PSv1_2, PSv2_2, PSv3_2, PSv4_2, PSv5_2, PSv1_3, PSv2_3, PSv3_3, PSv4_3, PSv5_3, PSv1_4, PSv2_4, PSv3_4, PSv4_4, PSv5_4, PSv1_5, PSv2_5, PSv3_5, PSv4_5, PSv5_5];
SRv = [SRv1_1, SRv2_1, SRv3_1, SRv4_1, SRv5_1, SRv1_2, SRv2_2, SRv3_2, SRv4_2, SRv5_2, SRv1_3, SRv2_3, SRv3_3, SRv4_3, SRv5_3, SRv1_4, SRv2_4, SRv3_4, SRv4_4, SRv5_4, SRv1_5, SRv2_5, SRv3_5, SRv4_5, SRv5_5];


%Line of Best Fit
[ypred,sigma_ypred,slope,yintercept,Slopeuncertainty,Yinterceptuncertainty,Q]=linefit(V',k');
figure
hold on
scatter(V,k)
plot(V, ypred)
title('Loss Coefficient vs. Voltage')
xlabel('Voltage [V]')
ylabel('Loss Coefficient, k')
legend('k vs. Voltage','Line of Best Fit','location', 'northwest')
hold off

[ypred1,sigma_ypred1,slope1,yintercept1,Slopeuncertainty1,Yinterceptuncertainty1,Q1]=linefit(velocity',k');
figure
hold on
scatter(velocity,k)
plot(velocity, ypred1)
title('Loss Coefficient vs. Airspeed')
xlabel('Airspeed [m/s]')
ylabel('Loss Coefficient, k')
legend('k vs. Velocity','Line of Best Fit', 'location', 'northwest')
hold off




% Plot calculated values against one another

figure
scatter(V,velocity)
title('Airspeed vs. Voltage')
xlabel('Voltage [V]')
ylabel('Airspeed [m/s]')


figure
scatter(V,k)
title('Loss Coefficient vs. Voltage')
xlabel('Voltage [V]')
ylabel('Loss Coefficient, k')

figure
scatter(velocity,k)
title('Loss Coefficient vs. Airspeed')
xlabel('Airspeed [m/s]')
ylabel('Loss Coefficient, k')

% For comparison purposes, a table for calculated airspeeds with their
% respective voltages and the calculated loss coefficeint at each airspeed
VoltageVelocityKvalue = [V' SRv' k' velocity' PSv'];
Values = sortrows(VoltageVelocityKvalue,1);

fprintf('               Static Ring        Loss       Static Ring Airspeed with    Pitot-Static\n')
fprintf('Voltage [V]   Airspeed [m/s]   Coefficient    Loss Coefficient [m/s]     Airspeed [m/s]\n')
fprintf('===========   ==============   ===========   =========================   ==============\n')
fprintf('%11.2f %16.4f %13.4f %27.4f %16.4f\n' , Values')
