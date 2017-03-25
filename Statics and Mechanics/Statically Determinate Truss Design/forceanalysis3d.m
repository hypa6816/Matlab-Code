function [barforces,reacforces]=forceanalysis3d(joints,connectivity,reacjoints,reacvecs,loadjoints,loadvecs)
% function [barforces,reacforces]=forceanalysis3d(joints,connectivity,reacjoints,reacvecs,loadjoints,loadvecs)
%
% compute forces in bars and reaction forces
%
% input:  joints       - coordinates of joints
%         connectivity - connectivity 
%         reacjoints   - joint id where reaction acts on
%         reacvecs     - unit vector associated with reaction force
%         loadjoints   - joint id where external load acts on
%         loadvecs     - load vector
%
% output: barforces    - force magnitude in bars
%         reacforces   - reaction forces
%
% 

% extract number of joints, bars, reactions, and loads
numjoints = size(joints,1);
numbars   = size(connectivity,1);
numreact  = size(reacjoints,1);
numloads  = size(loadjoints,1);

% number of equations
numeqns = 3 * numjoints;

% allocate arrays for linear system
Amat = zeros(numeqns);
bvec = zeros(numeqns,1);
%input known values for finding the total weight at each joint
magnetlength = .0127; % [m]
diameterjoint = .0127; % [m]
weightmagnet = .0018; % [kg]
weightjoint = .0084; % [kg]
densitybar = .00114; %[m/kg^3]
gravity = 9.81; % [m/s^2]

% allocate sum of the weight at each joint
sumweightatjoint = weightjoint.*ones(numjoints,1);



% build Amat - loop over all joints

for i=1:numjoints
    
   % equation id numbers
   idx = 3*i-2;
   idy = 3*i-1;
   idz = 3*i;
   
   % get all bars connected to joint
   [ibar,ijt]=find(connectivity==i);
   
   % loop over all bars connected to joint
   for ib=1:length(ibar)
       
       % get bar id
       barid=ibar(ib);
       
       % get coordinates for joints "i" and "j" of bar "barid"
       joint_i = joints(i,:);
       if ijt(ib) == 1
           jid = connectivity(barid,2);
       else ijt(ib)
           jid = connectivity(barid,1);
       end

       joint_j = joints(jid,:);
       % Find lengths of each member connected the the joint
       x1 = joint_i(1);
       y1 = joint_i(2);
       z1 = joint_i(3);
       
       x2 = joint_j(1);
       y2 = joint_j(2);
       z2 = joint_j(3);
       % Find the true length utlilizing the Pathogreom Theorem
       truelength = sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2);
       % To find length of the bar alone subtrant magnet length and radius
       % of the joint on both sides
       lengthbar = truelength - magnetlength -diameterjoint;
       % Find the weight of each bar
       weightbar = 1/2*lengthbar*densitybar*gravity;
       % Add the weight of the bar and magnet
       sumweightbar = weightbar + weightmagnet;
       % Add the weight of the bar and magnet to a running sum to compute
       % the sum of the weight forces on each joint
       sumweightatjoint(i) = sumweightatjoint(i) + sumweightbar;
       
       
       % compute unit vector pointing away from joint i
       vec_ij = joint_j - joint_i;
       uvec   = vec_ij/norm(vec_ij);
       
       % add unit vector into Amat
       Amat([idx idy idz],barid)=uvec;
   end
end

% build contribution of support reactions 
for i=1:numreact
    
    % get joint id at which reaction force acts
    jid=reacjoints(i);

    % equation id numbers
    idx = 3*jid-2;
    idy = 3*jid-1;
    idz = 3*jid;

    % add unit vector into Amat
    Amat([idx idy idz],numbars+i)=reacvecs(i,:);
end

% build load vector
for i=1:numloads
    
    % get joint id at which external force acts
    jid=loadjoints(i);

    % equation id numbers
    idx = 3*jid-2;
    idy = 3*jid-1;
    idz = 3*jid;

    % add unit vector into bvec (sign change)
    bvec([idx idy idz])=-loadvecs(i,:);
end

%check for invertability of Amat
if rank(Amat) ~= numeqns
    error('Amat is rank defficient: %d < %d\n',rank(Amat),numeqns);
end

%For loop to add in the weights into the z compartment for each joint of the bvec (external forces)
for i = 1:numjoints;
    bvec(3*i) = bvec(3*i)+sumweightatjoint(i);
end

% solve system
xvec=Amat\bvec;

% extract forces in bars and reaction forces
barforces=xvec(1:numbars);
reacforces=xvec(numbars+1:end);

end