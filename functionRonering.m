function R = functionRonering(M,theta,spreadStdDeg,antennaDist)
%Generate covariance matrix with the one-ring model.
%
%This Matlab function is used in the article:
%
%Emil Bjornson, Jakob Hoydis, Luca Sanguinetti, ?Massive MIMO has Unlimited
%Capacity,? IEEE Transactions on Wireless Communications, to appear.
%
%Download article: https://arxiv.org/pdf/1705.00538
%
%This is version 1.0 (Last edited: 2017-10-27)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original article listed above.
%
%
%INPUT:
%M             = Number of antennas
%theta         = Main angle
%spreadStdDeg  = Standard deviation of the angular spread around the main
%                angle (measured in degrees)
%antennaDist   = Distance between antennas (in wavelengths)
%
%OUTPUT:
%R             = M x M channel covariance matrix


%Compute angular spread in radians
spreadStd = spreadStdDeg*pi/180;

%If the antenna distance is not specified
if  nargin < 4
    
    %Half a wavelength distance
    antennaDist = 1/2;
    
end

%The covariance matrix has the Toeplitz structure, so we only need to
%compute the first row.
firstRow = zeros(M,1);

%Go through all columns in the first row
for col = 1:M
    
    %Distance from the first antenna
    distance = antennaDist*(col-1);
    
    %Set the upper and lower limit of the uniform distribution
    limits = sqrt(3)*spreadStd;
    
    %Define integrand of Eq. (XX)
    F = @(Delta)exp(1i*2*pi*distance*sin(theta+Delta))/(2*limits);
    
    %Compute the integral in Eq. (XX), over the entire interval
    firstRow(col) = integral(F,-limits,limits);
    
end

%Compute the covarince matrix by utilizing the Toeplitz structure
R = toeplitz(firstRow);
