function [Cp] = pressureDistribution(ports,pitotDynamic,xpp,angle)
numGroups = 4;

c = 3.5; %cord length in inches
Cp = zeros(9,17,numGroups);
for i = 1:9
    for j = 1:17
        for k = 1:numGroups
            Cp(i,j,k) = ports(i,j,k)/pitotDynamic(i,1,k);
        end
    end
end
cordWise = xpp/c;

%pressure distribution by angle of attack
cc=hsv(numGroups*3);
figure(1)
hold on
axis();
for i=1:numGroups
    for j=1:3
        plot(cordWise,Cp(j*3-2,:,i),'o--','color',cc(i*3+j-3,:))
        ylabel('C_{p}')
        xlabel('Chord Position')
        legend_strings{i*3+j-3} = sprintf('Angle: %d',angle(j*3-2,1,i));
    end
    
end
legend(legend_strings)
end

