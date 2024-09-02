tic %Start Code-Execution
clc %Clear Terminal
clear %Clear Workspace
close all %Close Side Operations
snrdB=-5:1:30; %Iteration of SNR(dB) : Channel-Variation
snr=10.^(snrdB/10); %SNR from SNR(dB)
ne=zeros(1,numel(snrdB));
N=input('Enter sample no. (no. of communications) : '); %Sample-Number
for j=1:numel(snr) %Iteration (upto Individual SNR)
    sigma=sqrt(1/(2*snr(j))); %Standard Deviation of Gaussian-Noise
    E=0; %Error-count Initiated
    for i=1:N %Iteration for Individual Samples
        %Transmitter-modelling
        d=round(rand(1));%Random bit-generation 
        x=(2*d)-1;%Unipolar(Bits)->Bipolar(Signal)
        %Channel-modelling
        n=sigma*randn;%Noise-Signal
        %Receiver-modelling
        y=x+n;%Noise-corrupted Signal
        %Decision-Logic (at End)
        if y>=0 
            d_rcv=1;
        else
            d_rcv=0;
        end
        %Error-Detection
        if(d_rcv~=d) %Bit-Mismatch (After Communicaton)
            E=E+1; %Error-count Upgrading
        end
    end
    ne(j)=E; %j-th Error, with respect to the j-th SNR (Array-Generation)
end
ber_sim=(ne/N); %Simulation-BER (from Error-count)
ber_theo=qfunc(sqrt(2*snr)); %Theoretical-BER (from SNR)
semilogy(snrdB,ber_sim,snrdB,ber_theo,'--'); %Semi-log Plot (Log. in Y-axis)
axis([-5 30 0.0001 1]); %Axis-Mapping (with Thickness=1)
xlabel('SNR in dB'); %X-Label
ylabel('BER'); %Y-Label
legend('Simulation','Theoretical'); %Naming of the Different Plots
toc %End Code-execution : Calculate Total Time Elapsed

