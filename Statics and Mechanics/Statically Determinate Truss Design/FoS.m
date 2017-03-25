function FOS = FoS(FailureProbability)
% Calculate the factor of safety for our truss given:
% Strength of ball joints: 4.8 ± .4 N.

%% Fdsr = icdf('normal',Pdsr,µF,sigF)
Pdsr = FailureProbability; %utilizing the failure probability found in force anaylsis
uF = 4.8; % Strength of ball joints
sigF = .4; % Error in the strength of ball joints
Fdsr = icdf('Normal',Pdsr,uF,sigF);

%% Calculating n, the number of Monte Carlo Trials
n = (uF-Fdsr)/sigF;

%% Calculating the Factor of Safety
FOS = uF/(uF-n*sigF);
fprintf('The Factor of Safety is %f',FOS)
end