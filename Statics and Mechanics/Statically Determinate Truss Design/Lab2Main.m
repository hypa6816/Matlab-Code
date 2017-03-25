% 2001 
% Lab 2
% Group 34
function Lab2Main(inputfile,outputfile)
%inputfile = 'test3d_1.inp'
%outputfile = 'test3d_1.out'
% Load all information using function read input
[joints,connectivity,reacjoints,reacvecs,loadjoints,loadvecs]=readinput3d(inputfile);
% Analyze functions using function force anaalysis.
[barforces,reacforces]=forceanalysis3d(joints,connectivity,reacjoints,reacvecs,loadjoints,loadvecs);
% Write an output file
writeoutput3d(outputfile,inputfile,barforces,reacforces,joints,connectivity,reacjoints,reacvecs,loadjoints,loadvecs)
% Find the failure Probability using the monte carlo simulation
FailureProbability=truss3dmcs(inputfile);
% Convert all joint measurements from metric standards to U.S. standards
joints = joints./0.0254;
% Plot the truss
plottruss(joints,connectivity,barforces,reacjoints,3*[0.025,0.04,0.05],[0 0 0 1])
% Find the Safety Factor using the failure probability
FOS = FoS(FailureProbability);
