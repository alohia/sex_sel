function [ X, Hk_data, n ] = getdata( caste )
%getdata This function gets the data for a caste
%   It returns the X = [income, no. of boys] and Hk

% Import data
%set(0,'DefaultFigureVisible','off');  % all subsequent figures "off"
filepath = strcat('./data/713_10class/caste', num2str(caste), '.csv');
mydata = csvread(filepath, 1);
data = struct('inc', [], 'm', [],'f', [],'t', [], 'csr', [], 'Hk', [], 'X', []);
data.inc = mydata(:,1);
data.m  = mydata(:,2);
data.f  = mydata(:,3);
data.t = mydata(:,4);
data.csr = mydata(:,5);
data.Hk = mydata(:,6);
data.X = mydata(:,7);
n = length(data.inc);
X = [data.inc data.X];
Hk_data = data.Hk;
end

