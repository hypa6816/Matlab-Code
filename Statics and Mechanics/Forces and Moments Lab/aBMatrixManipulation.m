function [A,b] = aBMatrixManipulation(numberOfForces, forcesCoordinates, externalForces, externalMoments, locationSupports, uVecReactionDir, momentOrForce)
%preallocating the matrices to improve speed
A = zeros(6);
b = zeros(6,1);
%creating a for loop that adds all the external forces together
%the moments are just force coordinates cross the external forces 
%for loop always this process to be done index by index
for i = 1:numberOfForces(1,1)
    b(1) = b(1) + externalForces(i,1);
    b(2) = b(2) + externalForces(i,2);
    b(3) = b(3) + externalForces(i,3);
    Mi = cross(forcesCoordinates(i,:), externalForces(i,:));
    b(4) = b(4) + Mi(1);
    b(5) = b(5) + Mi(2);
    b(6) = b(6) + Mi(3);
end
for i = 1:numberOfForces(1,2)
    b(4) = b(4) + externalMoments(i,1);
    b(5) = b(5) + externalMoments(i,2);
    b(6) = b(6) + externalMoments(i,3);
end
b = -b;

numReactionForces = sum(strcmp(momentOrForce,'F'));
numReactionMoments = sum(strcmp(momentOrForce,'M'));


reactionForcesResultantMoments = zeros(numReactionForces,3);
for i = 1:numReactionForces
    reactionForcesResultantMoments(i,:) = cross(locationSupports(i,:),uVecReactionDir(i,:));
end

for i = 1:3
    for j = 1:numReactionForces
        A(i,j) = uVecReactionDir(j,i);
    end
end
for i = 4:6
    for j = 1:numReactionForces
        A(i,j) = reactionForcesResultantMoments(j,i-3);
    end
    for j = 1:numReactionMoments
        A(i,j+numReactionForces) = uVecReactionDir(j+numReactionForces,i-3);
    end
end