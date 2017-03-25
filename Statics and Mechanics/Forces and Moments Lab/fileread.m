%Sara Park
%2001 Statics and Mechanics
%Lab 1
%Reading data 


clear all

%make fileID
fileID = fopen('dataFile2.txt');
%read data into matlab using textscan; skip lines starting with #
%since they are comments and are not data that is needed
C = textscan(fileID,'%s','CommentStyle','#');
fclose(fileID);


%take the first two inputs and rename to forces and moments
Forces = str2num(C{1}{1});
Moments = str2num(C{1}{2});

%Make zero matrices for all needed values
[LocationForces] = zeros(Forces,3);
[ForcesMagDir] = zeros(Forces,4);
[LocationMoments] = zeros(Moments,3);
[MomentsMagDir] = zeros(Moments,4);
[LocationSupports] = zeros(6,3);
[ForcesReactionDir] = zeros(Forces,3);
[MomentsReactionDir] = zeros(Moments,3);

%For loop to input the Location of Forces into an organized Matrix
%2 is needed because there are 2 points that start the file that are not
%needed in any matrices. It is then multiplied by i to find the point prior
%to the first point needed. Then j is added to find the point needed. A
%counter is then added in order to pass onto the next row.
counter = 0;
for i = 1:Forces
    for j = 1:3
        LocationForces(i,j)=str2num(C{1}{2*i+j+counter});
    end
    %counter is needed after the first row
    counter = counter + 1;
end

%For loop to input the Magnitude and Direction of the Forces
%The counter starts at the last point that was inputed into the prior
%matrix.
%Another counter is needed to pass through the rows since there is 4
%columns now.
counter = Forces*3;
counter2 = 0;
for i = 1:Forces
    for j = 1:4
        ForcesMagDir(i,j)=str2num(C{1}{2*i+j+counter+counter2});
    end
    counter = counter + 1;
    counter2 = counter2 + 1;
end

%For loop to input the Location of Moments
counter = Forces*3 + Forces*4;
for i = 1:Moments
    for j = 1:3
        LocationMoments(i,j)=str2num(C{1}{2*i+j+counter});
    end
    counter = counter + 1;
   
end

%For loop to input magnitude and direction of external couple moments

counter = Forces*3 + Forces*4 + Moments*3;
counter2 = 0;
for i = 1:Moments
    for j = 1:4
        MomentsMagDir(i,j)= str2num(C{1}{2*i+j+counter+counter2});
    end
    counter = counter + 1;
    counter2 = counter2 + 1;
end

%For loop to input the location of supports
counter = Forces*3 + Forces*4 + Moments*3 + Moments*4;
for i= 1:6
    for j = 1:3
        LocationsSupports(i,j) = str2num(C{1}{2*i+j+counter});
    end
    counter = counter +1;
end

%For loop to input the Direction
counter = Forces*3 +Forces*4 + Moments*3 + Moments*4 + 6*3;
counter2 = 1;
counter3=0;
while C{1}{41+counter3}== 'F';
    for j= 1:3
        ForcesReactionDir(counter2,j)= str2num(C{1}{41+counter2})
        counter2 = counter2 +1;
    end
    
    counter3 = counter3 +1;
end

while C{1}{i+2+counter} == 'M';
    for j= 2:4
        MomentsReactionDir(i,j)= str2num(C{1}{2*i+j+counter+counter2})
    end
    counter = counter +1;
    counter2 = counter2 +1;
end