filename = fopen('dataFile.txt');
line = fgetl(filename);
n=1;
arrayOfLines = cell(30,1);
while ischar(line)
    splitLine = strread(line, '%s');
    if strcmp(splitLine(1), '#') == 0;
        arrayOfLines{n} = splitLine;
        n = n + 1;
    end
    line = fgetl(filename);
end
fclose(filename);



numberOfForces = str2double(arrayOfLines{1,1})';
forcesCoordinates = zeros(numberOfForces(1,1),3);
for i = 1:numberOfForces(1,1)
    forcesCoordinates(i,:) = str2double(arrayOfLines{1+i,:})';
end

magDirExternalForces = zeros(numberOfForces(1,1),4);
for i = 1:numberOfForces(1,1)
    magDirExternalForces(i,:) = str2double(arrayOfLines{(2+numberOfForces(1,1)+i),:})';
end

locationExternalCoupleMoments = zeros(numberOfForces(1,2),3);
for i = 1:numberOfForces(1,2)
    locationExternalCoupleMoments(i,:) = str2double(arrayOfLines{(2+2*numberOfForces(1,1)+i),:})';
end

magDirExternalCoupleMoments = zeros(numberOfForces(1,2),4);
for i = 1:numberOfForces(1,2)
    magDirExternalCoupleMoments(i,:) = str2double(arrayOfLines{(3+2*numberOfForces(1,1)+numberOfForces(1,2)+i),:})';
end

locationSupports = zeros(6,3);
for i = 1:6
    locationSupports(i,:) = str2double(arrayOfLines{(3+2*numberOfForces(1,1)+2*numberOfForces(1,2)+i),:})';
end

holdCell = cell(6,1);
type = cell(6,1);
dirReaction = zeros(6,3);
for i = 1:6
    holdCell{i,:} = arrayOfLines{(9+2*numberOfForces(1,1)+2*numberOfForces(1,2)+i),1}';
    type(i,1) = cellstr(holdCell{i}{1,1});
    for j = 1:3
        dirReaction(i,j) = str2double(holdCell{i}{1,j+1});
    end
end

% Put external forces and moments into cartesian form
externalForces = zeros(numberOfForces(1,1),3);
externalMoments = zeros(numberOfForces(1,2),3);
uVecReactionDir = zeros(6,3);

for i = 1:numberOfForces(1,1)
    externalForces(i,:) = magDirExternalForces(i,1)*magDirExternalForces(i,2:4)/norm(magDirExternalForces(i,2:4));
end
for i = 1:numberOfForces(1,2)
    externalMoments(i,:) = magDirExternalCoupleMoments(i,1)*magDirExternalCoupleMoments(i,2:4)/norm(magDirExternalCoupleMoments(i,2:4));
end
for i = 1:6
    uVecReactionDir(i,:) = dirReaction(i,:)/norm(dirReaction(i,:));
end
    delete dirReaction;
    
    
A = zeros(6);
b = zeros(6,1);

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

reactionForces = sum(strcmp(type,'F'));
reactionMoments = sum(strcmp(type,'M'));
for i = 1:3
    for j = 1:reactionForces
        A(i,j) = uVecReactionDir(
    end
end

for i = 1:6
    for j = 1:reactionMoments
        
    end
end


























