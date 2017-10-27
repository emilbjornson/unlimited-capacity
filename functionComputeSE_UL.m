function [SE_MR,SE_MMMSE,SE_SMMSE,SE_MZF,powerlevels] = functionComputeSE_UL(Hhat,C,R,tau_c,tau_p,nbrOfRealizations,M,K,L,p)
%Compute uplink SE for different receive combining schemes.
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
%Hhat              = M x nbrOfRealizations x K x L x L matrix with the MMSE
%                    channel estimates
%C                 = M x M x K x L x L matrix with estimation error
%                    correlation matrices when using MMSE estimation
%R                 = M x M x K x L x L matrix with spatial correlation
%                    matrices
%tau_c             = Length of coherence block
%tau_p             = Length of pilot sequences
%nbrOfRealizations = Number of channel realizations
%M                 = Number of antennas per BS
%K                 = Number of UEs per cell
%L                 = Number of BSs and cells
%p                 = Uplink transmit power per UE (same for everyone)
%
%OUTPUT:
%SE_MR       = K x L matrix where element (k,l) is the uplink SE of UE k in
%              cell l achieved with MR combining
%SE_MMMSE    = Same as SE_MR but with M-MMSE combining
%SE_SMMSE    = Same as SE_MR but with S-MMSE combining
%SE_MZF      = Same as SE_MR but with M-ZF combining
%powerlevels = 4 x 3 matrix where each row is one scheme, the first column
%              is the desired signal power, second column is interference
%              power from users with the same pilot and the third column is
%              interfernece from users with different pilots. This part
%              only supports K=2.


%Store identity matrices of different sizes
eyeM = eye(M);
eyeKL = eye(K*L);

%Compute the pre-log factor (normalized with number of channel realizations)
%assuming only uplink transmission
prelogFactor = (tau_c-tau_p)/(tau_c*nbrOfRealizations);

%Compute sum of estimation error covariance matrices for all cells
C_totM = reshape(p*sum(sum(C,3),4),[M M L]);
C_totM1 = reshape(p*sum(C(:,:,1,:,:),4),[M M L]); %Only for user 1
C_totM2 = reshape(p*sum(C(:,:,2,:,:),4),[M M L]); %Only for user 2

%Compute sum of single-cell estimation error covariance matrices for each
%cell
CR_totS = zeros(M,M,L);

for j = 1:L
    CR_totS(:,:,j) = p*(sum(C(:,:,:,j,j),3)+sum(sum(R(:,:,:,[1:j-1 j+1:end],j),3),4));
end


%Prepare to store simulation results
SE_MR = zeros(K,L);
SE_SMMSE = zeros(K,L);
SE_MMMSE = zeros(K,L);
SE_MZF = zeros(K,L);
powerlevels = zeros(4,3);


%% Go through all channel realizations
for n = 1:nbrOfRealizations
    
    %Go through all cells
    for j = 1:L
        
        %Extract channel estimate realizations from all users to BS j
        Hallj = reshape(Hhat(:,n,:,:,j),[M K*L]);
        
        %Compute MR matrix
        V_MR = Hallj(:,K*(j-1)+1:K*j);
        
        %Compute M-MMSE matrix
        V_MMMSE = p*(p*(Hallj*Hallj')+C_totM(:,:,j)+eyeM)\V_MR;
        
        %Compute S-MMSE matrix
        V_SMMSE = p*(p*(V_MR*V_MR')+CR_totS(:,:,j)+eyeM)\V_MR;
        
        %Compute M-ZF matrix
        V_MZF = (Hallj/((Hallj'*Hallj+1e-14*eyeKL)))*eyeKL(:,K*(j-1)+1:K*j);
        
        
        %Go through all users in cell j
        for k = 1:K
            
            %MR combining
            v = V_MR(:,k); %Extract combining vector
            
            %Compute signal and interference+noise terms
            signal = p*abs(v'*Hhat(:,n,k,j,j))^2;
            interfnoise = p*sum(abs(v'*Hallj).^2) + v'*(C_totM(:,:,j)+eyeM)*v - signal;
            
            %Compute instantaneous achievable SE for one realization
            SE_MR(k,j) = SE_MR(k,j) + prelogFactor*real(log2(1+signal/interfnoise));
            
            
            %Extract power levels
            if nargout>4
                
                %Channel realizations of user 1 in each cell
                Hhatallj1 = reshape(Hhat(:,n,1,:,j),[M L]);
                
                %Channel realizations of user 2 in each cell
                Hhatallj2 = reshape(Hhat(:,n,2,:,j),[M L]);
                
                %Noise amplification by receive combining
                noisegain = norm(v)^2;
                
                %Store desired signal power
                powerlevels(1,1) = powerlevels(1,1) + signal/noisegain/nbrOfRealizations;
                
                if k == 1 %If user 1
                    
                    %Interference from users with the same pilot
                    powerlevels(1,2) = powerlevels(1,2) + real(p*sum(abs(v'*Hhatallj1).^2) + v'*C_totM1(:,:,j)*v - signal)/noisegain/nbrOfRealizations;
                    
                    %Interference from users with other pilot
                    powerlevels(1,3) = powerlevels(1,3) + real(p*sum(abs(v'*Hhatallj2).^2) + v'*C_totM2(:,:,j)*v)/noisegain/nbrOfRealizations;
                    
                elseif k == 2 %If user 1
                    
                    %Interference from users with the same pilot
                    powerlevels(1,2) = powerlevels(1,2) + real(p*sum(abs(v'*Hhatallj2).^2) + v'*C_totM2(:,:,j)*v - signal)/noisegain/nbrOfRealizations;
                    
                    %Interference from users with other pilot
                    powerlevels(1,3) = powerlevels(1,3) + real(p*sum(abs(v'*Hhatallj1).^2) + v'*C_totM1(:,:,j)*v)/noisegain/nbrOfRealizations;
                    
                end
                
            end
            
            
            
            
            %M-MMSE combining
            v = V_MMMSE(:,k); %Extract combining vector
            
            %Compute signal and interference+noise terms
            signal = p*abs(v'*Hhat(:,n,k,j,j))^2;
            interfnoise = p*sum(abs(v'*Hallj).^2) + v'*(C_totM(:,:,j)+eyeM)*v - signal;
            
            %Compute instantaneous achievable SE for one realization
            SE_MMMSE(k,j) = SE_MMMSE(k,j) + prelogFactor*real(log2(1+signal/interfnoise));
            
            
            %Extract power levels
            if nargout>4
                
                %Noise amplification by receive combining
                noisegain = norm(v)^2;
                
                %Store desired signal power
                powerlevels(3,1) = powerlevels(3,1) + signal/noisegain/nbrOfRealizations;
                
                if k == 1 %If user 1
                    
                    %Interference from users with the same pilot
                    powerlevels(3,2) = powerlevels(3,2) + real(p*sum(abs(v'*Hhatallj1).^2) + v'*C_totM1(:,:,j)*v  - signal)/noisegain/nbrOfRealizations;
                    
                    %Interference from users with other pilot
                    powerlevels(3,3) = powerlevels(3,3) + real(p*sum(abs(v'*Hhatallj2).^2) + v'*C_totM2(:,:,j)*v)/noisegain/nbrOfRealizations;
                    
                elseif k == 2 %If user 2
                    
                    %Interference from users with the same pilot
                    powerlevels(3,2) = powerlevels(3,2) + real(p*sum(abs(v'*Hhatallj2).^2) + v'*C_totM2(:,:,j)*v  - signal)/noisegain/nbrOfRealizations;
                    
                    %Interference from users with other pilot
                    powerlevels(3,3) = powerlevels(3,3) + real(p*sum(abs(v'*Hhatallj1).^2) + v'*C_totM1(:,:,j)*v)/noisegain/nbrOfRealizations;
                    
                end
                
            end
            
            
            
            
            %S-MMSE combining
            v = V_SMMSE(:,k); %Extract combining vector
            
            %Compute signal and interference+noise terms
            signal = p*abs(v'*Hhat(:,n,k,j,j))^2;
            interfnoise = p*sum(abs(v'*Hallj).^2) + v'*(C_totM(:,:,j)+eyeM)*v - signal;
            
            %Compute instantaneous achievable SE for one realization
            SE_SMMSE(k,j) = SE_SMMSE(k,j) + prelogFactor*real(log2(1+signal/interfnoise));
            
            
            %Extract power levels
            if nargout>4
                
                %Noise amplification by receive combining
                noisegain = norm(v)^2;
                
                %Store desired signal power
                powerlevels(2,1) = powerlevels(2,1) + signal/noisegain/nbrOfRealizations;
                
                if k == 1 %If user 1
                    
                    %Interference from users with the same pilot
                    powerlevels(2,2) = powerlevels(2,2) + real(p*sum(abs(v'*Hhatallj1).^2) + v'*C_totM1(:,:,j)*v  - signal)/noisegain/nbrOfRealizations;
                    
                    %Interference from users with other pilot
                    powerlevels(2,3) = powerlevels(2,3) + real(p*sum(abs(v'*Hhatallj2).^2) + v'*C_totM2(:,:,j)*v)/noisegain/nbrOfRealizations;
                    
                elseif k == 2 %If user 2
                    
                    %Interference from users with the same pilot
                    powerlevels(2,2) = powerlevels(2,2) + real(p*sum(abs(v'*Hhatallj2).^2) + v'*C_totM2(:,:,j)*v  - signal)/noisegain/nbrOfRealizations;
                    
                    %Interference from users with other pilot
                    powerlevels(2,3) = powerlevels(2,3) + real(p*sum(abs(v'*Hhatallj1).^2) + v'*C_totM1(:,:,j)*v)/noisegain/nbrOfRealizations;
                    
                end
                
            end
            
            
            
            
            %M-ZF combining
            v = V_MZF(:,k); %Extract combining vector
            
            %Compute signal and interference+noise terms
            signal = p*abs(v'*Hhat(:,n,k,j,j))^2;
            interfnoise = p*sum(abs(v'*Hallj).^2) + v'*(C_totM(:,:,j)+eyeM)*v - signal;
            
            %Compute instantaneous achievable SE for one realization
            SE_MZF(k,j) = SE_MZF(k,j) + prelogFactor*real(log2(1+signal/interfnoise));
            
            
            %Extract power levels
            if nargout>4
                
                %Noise amplification by receive combining
                noisegain = norm(v)^2;
                
                %Store desired signal power
                powerlevels(4,1) = powerlevels(4,1) + signal/noisegain/nbrOfRealizations;
                
                if k == 1 %If user 1
                    
                    %Interference from users with the same pilot
                    powerlevels(4,2) = powerlevels(4,2) + real(p*sum(abs(v'*Hhatallj1).^2) + v'*C_totM1(:,:,j)*v  - signal)/noisegain/nbrOfRealizations;
                    
                    %Interference from users with other pilot
                    powerlevels(4,3) = powerlevels(4,3) + real(p*sum(abs(v'*Hhatallj2).^2) + v'*C_totM2(:,:,j)*v)/noisegain/nbrOfRealizations;
                    
                elseif k == 2 %If user 2
                    
                    %Interference from users with the same pilot
                    powerlevels(4,2) = powerlevels(4,2) + real(p*sum(abs(v'*Hhatallj2).^2) + v'*C_totM2(:,:,j)*v  - signal)/noisegain/nbrOfRealizations;
                    
                    %Interference from users with other pilot
                    powerlevels(4,3) = powerlevels(4,3) + real(p*sum(abs(v'*Hhatallj1).^2) + v'*C_totM1(:,:,j)*v)/noisegain/nbrOfRealizations;
                    
                end
                
            end
            
        end
        
    end
    
end
