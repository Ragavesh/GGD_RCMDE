% This is sample code to test GGD-DE, GGD-MDE & GGD-RCMDE
% and DE, MDE & RCMDE entropy values 
% the following values are from the literature [3]
clc;
clear;
close all;
rng(1,'twister');
x=[9,8,1,12,5,-3,2.5,8.01,2.99,4,-1,10];
DisEn_NCDF(x,2,3,1)  % Results 1.8462 as per lit [3]

%change the beta value as per requirement
DisEn_GCDF(x,2,3,1,2) % when beta=2, GGD - normal distribution

a=4; b=-4;% generate random number between +/- 10
x=a + (b-a).*rand(1,1024);  % sample sgnal length
RCMDE(x,2,3,1,15) % RCMDE values as per lit [2]
mod_RCMDE(x,2,3,1,15,2) % beta=2, GGD-normal distribution
