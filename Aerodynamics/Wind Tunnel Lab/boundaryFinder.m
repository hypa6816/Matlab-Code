function [] = boundaryFinder(vheight)

%Loops through the data from each port
for i = 1:5
    
    %Calculates the free-stream velocity by averaging the two velocities at
    %the center
    freeStream = mean(vheight(23:24,1,i));
    
    %This calculates the first index at which the velocity is %95 of the
    %free stream velocity
    boundaryIndices = find(0.95*freeStream <= vheight(:,:,i));
    
    %The boundary layer thus occurs at the following index
    boundaryIndex = boundaryIndices(1);
    
    %The boundary layer height is thus
    boundaryHeight = num2str(vheight(boundaryIndex,2,i));
    
    %Port#
    port = num2str(i+6);
    
    %Prints out the results
    fprintf(['The height of the boundary layer for port ' port ' is ' boundaryHeight 'mm\n']);

end