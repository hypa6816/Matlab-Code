%% Facial Recognition Project
% Sara Park
% Elysia
% Jeremy
% Date created: 4/15/16
% Date edited: 4/22/16 

% Project Description: This project will explore the concept of facial
% recognition. 

% Functions included: 
% 1. Recognizing if there is a face or not. 
% 2. Comparing face to a previous face
% 3. 

%% Housekeeping:
close all;
clear clc;

%% Procedure to find Eigen faces:
% 1. Input a database of various photos (black and white)
% 2. Take the SVD approximation of the photos
% 3. Find 10 largest singular values and turn the rest to zero
% 4. Find signature vectors
% 5. Take the norm of vector 1 - vector 2

% Read in data from various pictures using standard 'photo_###'
NumberOfPhotos = 10;
Photo1 = imread('1.pgm');
[Row,Column]=size(Photo1);
Database = zeros(Row,Column,6);
for photo = 1:NumberOfPhotos
    filename = sprintf('%d.pgm',photo);
    Database(:,:,photo) = imread(filename);
    % Covert to double
    Database(:,:,photo) = double(Database(:,:,photo));
   

end

range = 7;

% Take the singular decomposition
[U_1,S_1,V_1] = svd(Database(:,:,1));
% Only take in the 10 largest singular values
S_1(range:end,range:end) = 0;
S_Matrices(:,:,1)=S_1;
% Multiplying to get lower Dimension in form: A_reduced = U*A*V_transpose
DatabaseReduced(:,:,1) = U_1*S_1*V_1';

[U_2,S_2,V_2] = svd(Database(:,:,2));
% Only take in the 10 largest singular values
S_2(range:end,range:end) = 0;
S_Matrices(:,:,2)=S_2;
% Multiplying to get lower Dimension in form: A_reduced = U*A*V_transpose
DatabaseReduced(:,:,2) = U_2*S_2*V_2';

[U_3,S_3,V_3] = svd(Database(:,:,3));
% Only take in the 10 largest singular values
S_3(range:end,range:end) = 0;
S_Matrices(:,:,3)=S_3;
% Multiplying to get lower Dimension in form: A_reduced = U*A*V_transpose
DatabaseReduced(:,:,3) = U_3*S_3*V_3';

[U_4,S_4,V_4] = svd(Database(:,:,4));
% Only take in the 10 largest singular values
S_4(range:end,range:end) = 0;
S_Matrices(:,:,4)=S_4;
% Multiplying to get lower Dimension in form: A_reduced = U*A*V_transpose
DatabaseReduced(:,:,4) = U_4*S_4*V_4';

[U_5,S_5,V_5] = svd(Database(:,:,5));
% Only take in the 10 largest singular values
S_5(range:end,range:end) = 0;
S_Matrices(:,:,5)=S_5;
% Multiplying to get lower Dimension in form: A_reduced = U*A*V_transpose
DatabaseReduced(:,:,5) = U_5*S_5*V_5';

[U_6,S_6,V_6] = svd(Database(:,:,6));
% Only take in the 10 largest singular values
S_6(range:end,range:end) = 0;
S_Matrices(:,:,6)=S_6;
% Multiplying to get lower Dimension in form: A_reduced = U*A*V_transpose
DatabaseReduced(:,:,6) = U_6*S_6*V_6';

[U_7,S_7,V_7] = svd(Database(:,:,7));
% Only take in the 10 largest singular values
S_7(range:end,range:end) = 0;
S_Matrices(:,:,7)=S_7;
% Multiplying to get lower Dimension in form: A_reduced = U*A*V_transpose
DatabaseReduced(:,:,7) = U_7*S_7*V_7';

[U_8,S_8,V_8] = svd(Database(:,:,8));
% Only take in the 10 largest singular values
S_8(range:end,range:end) = 0;
S_Matrices(:,:,8)=S_8;
% Multiplying to get lower Dimension in form: A_reduced = U*A*V_transpose
DatabaseReduced(:,:,8) = U_8*S_8*V_8';

[U_9,S_9,V_9] = svd(Database(:,:,9));
% Only take in the 10 largest singular values
S_9(range:end,range:end) = 0;
S_Matrices(:,:,9)=S_9;
% Multiplying to get lower Dimension in form: A_reduced = U*A*V_transpose
DatabaseReduced(:,:,9) = U_9*S_9*V_9';

[U_10,S_10,V_10] = svd(Database(:,:,10));
% Only take in the 10 largest singular values
S_10(range:end,range:end) = 0;
S_Matrices(:,:,10)=S_10;
% Multiplying to get lower Dimension in form: A_reduced = U*A*V_transpose
DatabaseReduced(:,:,10) = U_10*S_10*V_10';
%% Plotting the Reduced Images

figure();
for i=1:9
    subplot(3,3,i);
    imshow(uint8(DatabaseReduced(:,:,i)));
    name = sprintf('Face %d',i);
    title(name);
end
figure();
subplot(2,2,1);
imshow(uint8(Database(:,:,1)));
name = 'Face 1 (Original)';
title(name);
subplot(2,2,2);
imshow(uint8(DatabaseReduced(:,:,1)));
name1 = 'Face 1 (Reduced Quality)';
title(name1);





%% Finding signature vectors

SignatureVec1 = zeros(1,range);
SignatureVec2 = zeros(1,range);
SignatureVec3 = zeros(1,range);
SignatureVec4 = zeros(1,range);
SignatureVec5 = zeros(1,range);
SignatureVec6 = zeros(1,range);
SignatureVec7 = zeros(1,range);
SignatureVec8 = zeros(1,range);
SignatureVec9 = zeros(1,range);
SignatureVec10 = zeros(1,range);

for j = 1:range
    SignatureVec1(1,j) = S_1(j,j);
end
for j = 1:range
    SignatureVec2(1,j) = S_2(j,j);
end
for j = 1:range
    SignatureVec3(1,j) = S_3(j,j);
end
for j = 1:range
    SignatureVec4(1,j) = S_4(j,j);
end
for j = 1:range
    SignatureVec5(1,j) = S_5(j,j);
end
for j = 1:range
    SignatureVec6(1,j) = S_6(j,j);
end
for j = 1:range
    SignatureVec7(1,j) = S_7(j,j);
end
for j = 1:range
    SignatureVec8(1,j) = S_8(j,j);
end
for j = 1:range
    SignatureVec9(1,j) = S_9(j,j);
end
for j = 1:range
    SignatureVec10(1,j) = S_10(j,j);
end
%% test with newFace
SignatureVec_new = zeros(1,range);
newFace = imread('New.pgm');
newFace = double(newFace);
[U_new,S_new,V_new] = svd(newFace);
S_new(range:end,range:end)=0;
newFaceReduced = U_new*S_new*V_new';

for i = 1:range;
    SignatureVec_new(1,i) = S_new(i,i);
end



%comparing newFace to previous faces
norms = zeros(1,NumberOfPhotos);
% for i = 1: NumberOfPhotos
%     norms(:,i)= norm(SignatureVec_new - SignatureVec(i));
% end

norms(:,1) = norm(SignatureVec_new - SignatureVec1);
norms(:,2) = norm(SignatureVec_new - SignatureVec2);
norms(:,3) = norm(SignatureVec_new - SignatureVec3);
norms(:,4) = norm(SignatureVec_new - SignatureVec4);
norms(:,5) = norm(SignatureVec_new - SignatureVec5);
norms(:,6) = norm(SignatureVec_new - SignatureVec6);
norms(:,7) = norm(SignatureVec_new - SignatureVec7);
norms(:,8) = norm(SignatureVec_new - SignatureVec8);
norms(:,9) = norm(SignatureVec_new - SignatureVec9);
norms(:,10) = norm(SignatureVec_new - SignatureVec10);

[closestPicval,closestPic] = min(norms);

figure();
subplot(2,2,1);
imshow(uint8(DatabaseReduced(:,:,closestPic)));
name = sprintf('Face %.d',closestPic);
title(name);

subplot(2,2,2);
imshow(uint8(newFaceReduced));
name = sprintf('Face 10');
title(name);


        
     
    


