clear all; clc;
%close all;
% Create Rayleigh fading channel object.


chan = rayleighchan(1/1000000,50,1.0e-004 * [0 0.0400 0.0800 0.1200],[0 -3 -6 -9]);
% ts= 1/1000000 is the sample time of the input signal
% Maximum Doppler Shifts=50
%delayVector = 1.0e-004 * [0 0.0400 0.0800 0.1200];
%gainVector = [0 -3 -6 -9]; % Average path gains (dB)


% Generate data and apply fading channel.
M = 2; % DBPSK modulation order
hMod = comm.BPSKModulator; % Create a DPSK modulator
%modulates using the binary phase shift keying method
hDemod = comm.BPSKDemodulator; % Create a DPSK demodulator
% demodulates a signal that was modulated using the binary phase shift keying method

tx = randi([0 M-1],500,1); % Generate a random bit stream
% Make a random matrix of 500x1 with value as 0 and 1

dpskSig = step(hMod, tx); % DPSK modulate the signal
% demodulates input tx, with the BPSK demodulator hMod, and returns dpskSig
% if tx is 1 then dpskSig value is -1.0 + 0.0i
% if tx is 0 then dpskSig value is 1.0 + 0.0i

constellation(hDemod); 
%Displays the Constellation diagram of BPSK Demodulator

y = zeros(size(dpskSig));
for i=1:length(dpskSig)
    if dpskSig(i) == 1.0000+0.0000i
        y(i) = 1;
    else
        y(i) = -1;
    end
end

Y = eye(M); %Eye Diagram 
yxx = zeros(500,2); %Initializing with 0 to 500x2 matrix
for i=1:length(dpskSig)
    if y(i) == 1
        yxx(i,:) = Y(1,:);
    else
        yxx(i,:) = Y(2,:);
    end
end

fadedSig = filter(chan,dpskSig); % Apply the channel effects
%processes the baseband signal vector dpsk with the rayleigh channel object chan. 
%result is the signal vector fadedSig. 

% Compute error rate for different values of SNR
SNR = 0:2:20; % Range of SNR values, in dB.
numSNR = length(SNR); %length of the SNR vector created

berVec = zeros(3, numSNR);
% Create an AWGNChannel and ErrorRate calculator System object
hChan = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)');
hErrorCalc = comm.ErrorRate; %Communication Error rate calculation
ber2=zeros(numSNR,1);

for n = 1:numSNR %Looping over the SNR
    
hChan.SNR = SNR(n); %Initializing the SNR value to that of AWGN Channel
rxSig = step(hChan,fadedSig); % Add Gaussian noise by using hChan (AWGN Channel)
rx = step(hDemod, rxSig); % Demodulate using BPSK Demodulator
X = [real(rxSig), imag(rxSig)]; %Seperating the real and imaginary part
reset(hErrorCalc) %Reset the Error Calculation
% Compute error rate.
berVec(:,n) = step(hErrorCalc,tx,rx); %Computation of error rate
output2=mynn(X,yxx); %Passing it to find the signal which was sent
[aa,ii]=max(output2);
yt=(ii==1);
rxf = step(hDemod, 2*double(yt')-1);
ber2(n)=sum(abs(tx-rxf))/500;
end

BER = berVec(1,:); %Bit error
figure;
 %Compute theoretical performance results, for comparison.
BERtheory = berfading(SNR,'dpsk',M,1); 
% berfading calculates Bit error rate (BER) for Rayleigh fading channels
% SNR is the ratio of bit energy to noise power spectral density, in dB
% For DPSK

% Plot BER results
semilogy(SNR,BERtheory,'b-');
%Plotting the SNR values for the theoretical values

xlabel('SNR (dB)'); ylabel('BER');
title('Binary Signal over Rayleigh Fading Channel');
 %hold;
% semilogy(SNR,ber2,'g^-');
legend('BER');

tar = [real(dpskSig), imag(dpskSig)];
input_target = [X, tar];
input_target2 = [X, tar(:,1)];

