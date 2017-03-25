function [V1, kV1, V2, kV2, V3, kV3, V4, kV4, V5, kV5, velocity1, velocity2, velocity3, velocity4, velocity5, PSv1, PSv2, PSv3, PSv4, PSv5, SRv1, SRv2, SRv3, SRv4, SRv5] = Calc_v( data1, data2 )

avgData2 = mean(data2);
Patm = avgData2(1); % [K]
Tatm = avgData2(2); % [K]
R = 287; % [ J/kg/K]

[V1, kV1, V2, kV2, V3, kV3, V4, kV4, V5, kV5, avgSR1, avgSR2, avgSR3, avgSR4, avgSR5, avgPS1, avgPS2, avgPS3, avgPS4, avgPS5] = Calc_k(data1,data2);

velocity1 = sqrt(2*(avgSR1-kV1)*((R*Tatm)/Patm));
velocity2 = sqrt(2*(avgSR2-kV2)*((R*Tatm)/Patm));
velocity3 = sqrt(2*(avgSR3-kV3)*((R*Tatm)/Patm));
velocity4 = sqrt(2*(avgSR4-kV4)*((R*Tatm)/Patm));
velocity5 = sqrt(2*(avgSR5-kV5)*((R*Tatm)/Patm));

% Velocity based on pitot-static probe measurments
PSv1 = sqrt(2*avgPS1*((R*Tatm)/Patm));
PSv2 = sqrt(2*avgPS2*((R*Tatm)/Patm));
PSv3 = sqrt(2*avgPS3*((R*Tatm)/Patm));
PSv4 = sqrt(2*avgPS4*((R*Tatm)/Patm));
PSv5 = sqrt(2*avgPS5*((R*Tatm)/Patm));

% Velocity based on Static Ring measurements
SRv1 = sqrt(2*(avgSR1)*((R*Tatm)/Patm));
SRv2 = sqrt(2*(avgSR2)*((R*Tatm)/Patm));
SRv3 = sqrt(2*(avgSR3)*((R*Tatm)/Patm));
SRv4 = sqrt(2*(avgSR4)*((R*Tatm)/Patm));
SRv5 = sqrt(2*(avgSR5)*((R*Tatm)/Patm));

end