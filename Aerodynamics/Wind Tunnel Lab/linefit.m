function [ypred,sigma_ypred,slope,yintercept,Slopeuncertainty,Yinterceptuncertainty,Q]=linefit(x,y)
%This function makes a weighted line of best fit and the uncertainties in
%the y intercept and slope

% UID: L81F82BI76
% Date Created: 10/18/15
% Date Modified: 10/23/15

%inverting the materices to 
% Making matrices in form [A][Pls] = [y]
% First making matrix A: since y is linear, x and ones would be the only
% values
A = [x ones(length(x),1)];
% Then making the matrix of Pls which is equal to [A]^-1*y
Pls = inv(A'*A)*A'*y;
ypred = Pls(1)*x +  Pls(2);
%Pls meanings
slope = Pls(1);
yintercept = Pls(2);
sigma_ypred = sqrt((1/(length(y)-2))*sum((ypred-y).^2));
% finding the uncertainties in Y intercept and slope
% finding the Weight matrix since the uncertainty in Plt values is the diagnal of W
WDiag(1:length(x),1) = 1/(sigma_ypred^2);
W=diag(WDiag);
Q = inv(A'*W*A);
Slopeuncertainty = sqrt(Q(1,1));
Yinterceptuncertainty = sqrt(Q(2,2));







end

