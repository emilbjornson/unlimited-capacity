%This Matlab script can be used to generate all the figures in the article:
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


%Initialization
close all;
clear;


%Number of antennas
M = 100;

%Correlation factor in exponential correlation model
r = 0.5;

%Generate range of nominal AoAs
angles = linspace(0,2*pi,101);

%Angular deviation from nominal AoAs
Delta = 15;

%Prepare to save simulation results
eigenvaluesExponential = zeros(M,length(angles)-1);
eigenvaluesOneRing = zeros(M,length(angles)-1);


%%Go through all AoAs
for k = 1:length(angles)-1
    
    %Display simulation progress
    disp([num2str(k) ' angles out of ' num2str(length(angles)-1)]);
    
    %Generate covariance matrix with the exponential correlation model
    R = toeplitz((r*exp(1i*angles(k))).^(0:M-1));
    
    %Extract eigenvalues
    eigenvaluesExponential(:,k) = sort(eig(R),'descend');
    
    
    %Generate covariance matrix with the one-ring model
    R = functionRonering(M,angles(k)/2,Delta/sqrt(3),0.5);
    
    %Extract eigenvalues
    eigenvaluesOneRing(:,k) = sort(eig(R),'descend');
    
end


%Number of realizations of large-scale fading variations
realizations = 100;

%Standard deviation
stdvalue = 2;

%Generate eigenvalues with uncorrelated fading and large-scale fading
%variations over the array
fadingOverArray = randn(M,realizations);
eigenvaluesFadingUnsorted = 10.^(stdvalue*fadingOverArray/10);
eigenvaluesFading = sort(eigenvaluesFadingUnsorted,1,'descend');


%% Plot simulation results
figure; hold on; box on;

plot(1:M,mean(eigenvaluesOneRing,2),'b-.','LineWidth',1);

plot(1:M,mean(eigenvaluesExponential,2),'r--','LineWidth',1);

plot(1:M,mean(eigenvaluesFading,2),'k','LineWidth',1);

legend('One-ring model (\Delta=15^o)','Exponential correlation model (r=0.5)','Uncorrelated, fading variations (2 dB)','Location','SouthEast');

xlabel('Eigenvalue (decreasing order)');
ylabel('Normalized value');
set(gca,'Yscale','log');
ylim([1e-4 1e2]);
