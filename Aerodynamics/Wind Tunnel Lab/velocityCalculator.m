function [velocityHeight] = velocityCalculator(data)
%This calculates the velocity at the various heights in the inputted data
%file

samples = length(data)/12;

%The gas constant for air
R = 287; %[J/Kg/K]

%Recording the pressure difference into an easily usable matrix
pDifference = data(:,4);

%Recording the atmospheric tempreature into an easily usable matrix
TATM = data(:,2);

%Recording the atmospheric pressure into an easily usable matrix
PATM = data(:,1);

%Recording the height of ELD probe into an easily usable matrix
heightMat = data(:,6);

%Preallocating a matrix to store all the velocity data and the
%corresponding height data
velocityHeight = zeros(12,2);

for i = 1:12
    if (i == 1)
        pDelta = mean(pDifference(i:samples*i));
        tatmos = mean(TATM(i:samples*i));
        patmos = mean(PATM(i:samples*i));
        height = round(2*mean(heightMat(i:samples*i)))/2;
        
    else
        pDelta = mean(pDifference(samples*(i-1)+1:samples*(i)));
        tatmos = mean(TATM(samples*(i-1)+1:samples*(i)));
        patmos = mean(PATM(samples*(i-1)+1:samples*(i)));
        height = round(2*mean(heightMat(samples*(i-1)+1:samples*(i))))/2;
       
        
    end
    
    %Calculating the velocity at the point where the ELD probe is located
    v = sqrt(2*pDelta*(R*tatmos/patmos));
    
    %Storing the velocity and corresponding height in a single matrix
    velocityHeight(i,1) = v;
    velocityHeight(i,2) = height;

end

