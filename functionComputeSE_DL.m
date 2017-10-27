function [SE_MR,SE_MMMSE,SE_SMMSE,SE_MZF] = functionComputeSE_DL(H,Hhat,C,R,tau_c,tau_p,nbrOfRealizations,M,K,L,p,rho)
%Compute downlink SE for different transmit precoding schemes.
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
%H                 = M x nbrOfRealizations x K x L x L matrix with the
%                    channel realizations
%Hhat              = M x nbrOfRealizations x K x L x L matrix with the 
%                    channel estimates
%C                 = M x M x K x L x L matrix with estimation error
%                    correlation matrices
%R                 = M x M x K x L x L matrix with spatial correlation
%                    matrices
%tau_c             = Length of coherence block
%tau_p             = Length of pilot sequences
%nbrOfRealizations = Number of channel realizations
%M                 = Number of antennas per BS
%K                 = Number of UEs per cell
%L                 = Number of BSs and cells
%p                 = Uplink transmit power per UE (same for everyone)
%rho               = Downlink transmit power per UE (same for everyone)
%
%OUTPUT:
%SE_MR       = K x L matrix where element (k,l) is the downlink SE of UE k
%              in cell l achieved with MR precoding
%SE_MMMSE    = Same as SE_MR but with M-MMSE precoding
%SE_SMMSE    = Same as SE_MR but with S-MMSE precoding
%SE_MZF      = Same as SE_MR but with M-ZF precoding


%Store identity matrices of different sizes
eyeM = eye(M);
eyeKL = eye(K*L);

%Compute the prelog factor (normalized with number of realizations)
%assuming only uplink transmission
prelogFactor = (tau_c-tau_p)/(tau_c);

%Compute sum of estimation error covariance matrices for all cells
C_totM = reshape(p*sum(sum(C,3),4),[M M L]);

%Compute sum of single-cell estimation error covariance matrices for each
%cell
CR_totS = zeros(M,M,L);

for j = 1:L
    CR_totS(:,:,j) = p*(sum(C(:,:,:,j,j),3)+sum(sum(R(:,:,:,[1:j-1 j+1:end],j),3),4));
end


%Prepare to store simulation results for signal gains
signal_MR = zeros(K,L);
signal_SMMSE = zeros(K,L);
signal_MMMSE = zeros(K,L);
signal_MZF = zeros(K,L);


%Prepare to store simulation results for sum interference powers
interf_MR = zeros(K,L);
interf_SMMSE = zeros(K,L);
interf_MMMSE = zeros(K,L);
interf_MZF = zeros(K,L);


%% Go through all channel realizations
for n = 1:nbrOfRealizations
    
    %Go through all cells
    for j = 1:L
        
        %Extract channel realizations from all UEs to BS j
        Hallj = reshape(H(:,n,:,:,j),[M K*L]);
        
        %Extract channel estimate realizations from all UEs to BS j
        Hhatallj = reshape(Hhat(:,n,:,:,j),[M K*L]);
        
        %Compute MR matrix
        V_MR = Hhatallj(:,K*(j-1)+1:K*j);
        
        %Compute M-MMSE matrix
        if nargout > 1
            V_MMMSE = p*(p*(Hhatallj*Hhatallj')+C_totM(:,:,j)+eyeM)\V_MR;
        end
        
        %Compute S-MMSE matrix
        if nargout > 2
            V_SMMSE = p*(p*(V_MR*V_MR')+CR_totS(:,:,j)+eyeM)\V_MR;
        end
        
        %Compute M-ZF matrix
        if nargout > 3
            V_MZF = (Hhatallj/((Hhatallj'*Hhatallj+1e-14*eyeKL)))*eyeKL(:,K*(j-1)+1:K*j);
        end
        
        
        %Go through all UEs in cell j
        for k = 1:K
            
            if norm(V_MR(:,k))>0
                
                %MR precoding
                w = V_MR(:,k)/norm(V_MR(:,k)); %Extract precoding vector
                
                %Compute realizations of expectations in signal and
                %interference terms
                signal_MR(k,j) = signal_MR(k,j) + (w'*H(:,n,k,j,j))/nbrOfRealizations;
                interf_MR = interf_MR + rho*reshape(abs(w'*Hallj).^2,[K L])/nbrOfRealizations;
                
                
                %M-MMSE precoding
                if nargout > 1
                    
                    w = V_MMMSE(:,k)/norm(V_MMMSE(:,k)); %Extract precoding vector
                    
                    %Compute realizations of expectations in signal and
                    %interference terms
                    signal_MMMSE(k,j) = signal_MMMSE(k,j) + (w'*H(:,n,k,j,j))/nbrOfRealizations;
                    interf_MMMSE = interf_MMMSE + rho*reshape(abs(w'*Hallj).^2,[K L])/nbrOfRealizations;
                    
                end
                
                
                %S-MMSE precoding
                if nargout > 2
                    
                    w = V_SMMSE(:,k)/norm(V_SMMSE(:,k)); %Extract precoding vector
                    
                    %Compute realizations of expectations in signal and
                    %interference terms
                    signal_SMMSE(k,j) = signal_SMMSE(k,j) + (w'*H(:,n,k,j,j))/nbrOfRealizations;
                    interf_SMMSE = interf_SMMSE + rho*reshape(abs(w'*Hallj).^2,[K L])/nbrOfRealizations;
                    
                end
                
                
                %M-ZF precoding
                if nargout > 3
                    
                    w = V_MZF(:,k)/norm(V_MZF(:,k)); %Extract precoding vector
                    
                    %Compute realizations of expectations in signal and
                    %interference terms
                    signal_MZF(k,j) = signal_MZF(k,j) + (w'*H(:,n,k,j,j))/nbrOfRealizations;
                    interf_MZF = interf_MZF + rho*reshape(abs(w'*Hallj).^2,[K L])/nbrOfRealizations;
                    
                end
                
            end
            
        end
        
    end
    
end


%Compute SE with MR
SE_MR = prelogFactor*real(log2(1+(rho*abs(signal_MR).^2) ./ (interf_MR - rho*abs(signal_MR).^2 + 1)));

%Compute SE with M-MMSE
if nargout > 1
    SE_MMMSE = prelogFactor*real(log2(1+(rho*abs(signal_MMMSE).^2) ./ (interf_MMMSE - rho*abs(signal_MMMSE).^2 +1)));
end

%Compute SE with S-MMSE
if nargout > 2
    SE_SMMSE = prelogFactor*real(log2(1+(rho*abs(signal_SMMSE).^2) ./ (interf_SMMSE - rho*abs(signal_SMMSE).^2 + 1)));
end

%Compute SE with M-ZF
if nargout > 3
    SE_MZF = prelogFactor*real(log2(1+(rho*abs(signal_MZF).^2) ./ (interf_MZF - rho*abs(signal_MZF).^2 + 1)));
end
