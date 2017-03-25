function [V1, kV1, V2, kV2, V3, kV3, V4, kV4, V5, kV5, avgSR1, avgSR2, avgSR3, avgSR4, avgSR5, avgPS1, avgPS2, avgPS3, avgPS4, avgPS5] = Calc_k( data1, data2 )

V = data1(:,7);

V1 = V(1);
i=1; PS1=0; SR1=0;
while i<=20
    PS1 = PS1 + data1(i,3);% Differential Pressure from Pitot-Static Probe
    SR1 = SR1 + data2(i,3);% Differential Pressure from Static Pressure Ring
    i=i+1;
end
avgPS1 = PS1/(i-1); % Average the Differential Pressure from Pitot-Static Probe
avgSR1 = SR1/(i-1); % Average the Differential Pressure from Static Pressure Ring
kV1 = avgSR1-avgPS1; % k = (Patm-Ps)-(P0-Ps)

V2 = V(i);
PS2=0; SR2 = 0;
while i<=40
    PS2 = PS2 + data1(i,3);
    SR2 = SR2 + data2(i,3);
    i=i+1;
end
avgPS2 = PS2/(i-1);
avgSR2 = SR2/(i-1);
kV2 = avgSR2 - avgPS2;

V3 = V(i);
PS3=0; SR3 = 0;
while i<=60
    PS3 = PS3 + data1(i,3);
    SR3 = SR3 + data2(i,3);
    i=i+1;
end
avgPS3 = PS3/(i-1);
avgSR3 = SR3/(i-1);
kV3 = avgSR3 - avgPS3;

V4 = V(i);
PS4=0; SR4 = 0;
while i<=80
    PS4 = PS4 + data1(i,3);
    SR4 = SR4 + data2(i,3);
    i=i+1;
end
avgPS4 = PS4/(i-1);
avgSR4 = SR4/(i-1);
kV4 = avgSR4 - avgPS4;

V5 = V(i);
PS5=0; SR5 = 0;
while i<=100
    PS5 = PS5 + data1(i,3);
    SR5 = SR5 + data2(i,3);
    i=i+1;
end
avgPS5 = PS5/(i-1);
avgSR5 = SR5/(i-1);
kV5 = avgSR5 - avgPS5;

end
