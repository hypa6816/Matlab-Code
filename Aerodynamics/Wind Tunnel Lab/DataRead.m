function [ data ] = DataRead( filenameIn)
%DataRead reads in the data .csv file and imports it into an array
data = csvread(filenameIn,1);
end

