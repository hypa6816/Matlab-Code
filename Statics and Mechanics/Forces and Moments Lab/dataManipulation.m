function [numberOfForcesandMoments, forcesCoordinates, externalForces, externalMoments, locationSupports, uVecReactionDir, momentOrForce] = dataManipulation(arrayOfLines)


numberOfForcesandMoments = str2double(arrayOfLines{1,1})';
numberOfForces = numberOfForcesandMoments(1,1);
numberOfMoments = numberOfForcesandMoments(1,2);
forcesCoordinates = zeros(numberOfForcesandMoments(1,1),3);
for i = 1:numberOfForces
    forcesCoordinates(i,:) = str2double(arrayOfLines{1+i,:})';
end
position = 2*numberOfForces;
magDirExternalForces = zeros(numberOfForces,4);
for i = 1:numberOfForces
    magDirExternalForces(i,:) = str2double(arrayOfLines{(position+i),:})'
end
position = 1+2*numberOfForces+4*numberOfForces;
locationExternalCoupleMoments = zeros(numberOfMoments,3);
for i = 1:numberOfMoments
    locationExternalCoupleMoments(i,:) = str2double(arrayOfLines{(position+i),:})'
end
position= 1+2*numberOfForces+numberOfMoments;
magDirExternalCoupleMoments = zeros(numberOfForcesandMoments(1,2),4);
for i = 1:numberOfMoments
    magDirExternalCoupleMoments(i,:) = str2double(arrayOfLines{(position+i),:})'
end

locationSupports = zeros(6,3);
for i = 1:6
    locationSupports(i,:) = str2double(arrayOfLines{(1+2*numberofForces+2*numberOfMoments+i),:})';
end

holdCell = cell(6,1);
momentOrForce = cell(6,1);
dirReaction = zeros(6,3);
for i = 1:6
    holdCell{i,:} = arrayOfLines{(7+2*numberOfForces+2*numberOfMoments(1,2)+i),1}';
    momentOrForce(i,1) = cellstr(holdCell{i}{1,1});
    for j = 1:3
        dirReaction(i,j) = str2double(holdCell{i}{1,j+1});
    end
end

% Put external forces and moments into cartesian form
externalForces = zeros(numberOfForcesandMoments(1,1),3);
externalMoments = zeros(numberOfForcesandMoments(1,2),3);
uVecReactionDir = zeros(6,3);

for i = 1:numberOfForcesandMoments(1,1)
    externalForces(i,:) = magDirExternalForces(i,1)*magDirExternalForces(i,2:4)/norm(magDirExternalForces(i,2:4));
end
for i = 1:numberOfForcesandMoments(1,2)
    externalMoments(i,:) = magDirExternalCoupleMoments(i,1)*magDirExternalCoupleMoments(i,2:4)/norm(magDirExternalCoupleMoments(i,2:4));
end
for i = 1:6
    uVecReactionDir(i,:) = dirReaction(i,:)/norm(dirReaction(i,:));
end