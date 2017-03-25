%Setting the filename for input
%Setting a standard template
fileTemp = 'BoundaryLayer_S013_';

%Preallocate velocity and height matrix for EVERY HEIGHT AND PORT
vHeight = zeros(24,2,5);

%These are temporary placeholders
vh1 = zeros(24,2);
vh2 = zeros(24,2);

%Looping through the group names and storing the data
for i = 1:10
    
    if (i~=10)
        
        groupNumber = num2str(i);
        groupString = ['G0' groupNumber '.csv'];
        fileName = [fileTemp groupString];
        data = csvread(fileName,1);
        
    else
    
        fileName = [fileTemp 'G10.csv'];
        data = csvread(fileName,1);
        
    end
    
    
    
    %This uses the boundaryFinder function to collect data on the
    %velocities at various heights for the data collected
    vHeightTemp = velocityCalculator(data);
    if (mod(i,2) == 0)
        for j = 1:length(vHeightTemp)
            
            vh1(j+(j-1),:) = vHeightTemp(j,:);
            vHeight(:,:,(i/2)) = vh1+vh2;
        end
    else
        
        for j = 1:length(vHeightTemp)
            vh2(j*2,:) = vHeightTemp(j,:);
        end
        
        
    end
    
end

boundaryFinder(vHeight);





