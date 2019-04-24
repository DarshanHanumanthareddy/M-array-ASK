clc;
clear all;
close all;
format long;
SNR = 0:1:22; % Signal to Noise Ratio[dB]
N = 1e6; % total Number of symbols
M = 16; % Number of Symbols
k = log2(M); % Number of bits per Symbol
alphabet = [-(M-1) : 2 : (M-1)];
Eavg = (1/M)*(sum(alphabet.^2)); % average power of transmitted signal
Eb_N0_dB = SNR; % signal to noise ratio
Eb_N0 = 10.^(Eb_N0_dB/10);
Es_N0 = Eb_N0*k; % symbols energy to noise ratio
St = randsrc(1,N,alphabet);
St_norm = St/sqrt(Eavg);
% received signal
sigma2 = 1./(Es_N0);
N0 = sigma2;
for i = 1:length(SNR)
n = sqrt(sigma2(i)/2)*(randn(1,length(St))+j*randn(1,length(St)));
Sr_norm(i,:) = St_norm + n;
end
% Decision regions comparison
Sr = Sr_norm*sqrt(Eavg);
DM = [-(M-2):2:(M-2)]; % Decision Margins
for i = 1:length(SNR)
So_I(find(real(Sr(i,:)) < DM(1))) = alphabet(1);
if (length(DM) > 1)
for k = 2:length(DM)
So_I(find((real(Sr(i,:)) > DM(k-1))&(real(Sr(i,:)) < DM(k)))) = alphabet(k);
end
end
So_I(find(real(Sr(i,:)) > DM(length(DM)))) = alphabet(length(alphabet));
So(i,:) = So_I;
Pe16(i) = symerr(St,So(i,:))/N;
Pe16_analytic(i) = (2*(M-1)/M)*qfunc(sqrt((6*log2(M)/(M^2-1))*(Eb_N0(i)))); %analytical result
end
'M=16 Simulation Completed'
% M=8 simulation
M = 8; % Number of Symbols
k = log2(M); % Number of bits per Symbol
alphabet = [-(M-1) : 2 : (M-1)];
Eavg = (1/M)*(sum(alphabet.^2)); % average power of transmitted signal
Eb_N0_dB = SNR; % signal to noise ratio
Eb_N0 = 10.^(Eb_N0_dB/10);
Es_N0 = Eb_N0*k; % symbols energy to noise ratio
% M-ASK transmitted signal
St = randsrc(1,N,alphabet);
St_norm = St/sqrt(Eavg);
%received signal
sigma2 = 1./(Es_N0);
N0 = sigma2;
for i = 1:length(SNR)
n = sqrt(sigma2(i)/2)*(randn(1,length(St))+j*randn(1,length(St)));
Sr_norm(i,:) = St_norm + n;
end
%Decision regions comparison
Sr = Sr_norm*sqrt(Eavg);
DM = [-(M-2):2:(M-2)]; % Decision Margins
for i = 1:length(SNR)
So_I(find(real(Sr(i,:)) < DM(1))) = alphabet(1);
if (length(DM) > 1)
for k = 2:length(DM)
So_I(find((real(Sr(i,:)) > DM(k-1))&(real(Sr(i,:)) < DM(k)))) = alphabet(k);
end
end
So_I(find(real(Sr(i,:)) > DM(length(DM)))) = alphabet(length(alphabet));
So(i,:) = So_I;
Pe8(i) = symerr(St,So(i,:))/N;
Pe8_analytic(i) = (2*(M-1)/M)*qfunc(sqrt((6*log2(M)/(M^2-1))*(Eb_N0(i)))); %analytical result
end
'M=8 Simulation Completed'
%M=4 simultion
M = 4; % Number of Symbols
k = log2(M); % Number of bits per Symbol
alphabet = [-(M-1) : 2 : (M-1)];
Eavg = (1/M)*(sum(alphabet.^2)); % average power of transmitted signal
Eb_N0_dB = SNR; % signal to noise ratio
Eb_N0 = 10.^(Eb_N0_dB/10);
Es_N0 = Eb_N0*k; % symbols energy to noise ratio
% M-ASK transmitted signal
St = randsrc(1,N,alphabet);
St_norm = St/sqrt(Eavg);
% received signal
sigma2 = 1./(Es_N0);
N0 = sigma2;
for i = 1:length(SNR)
n = sqrt(sigma2(i)/2)*(randn(1,length(St))+j*randn(1,length(St)));
Sr_norm(i,:) = St_norm + n;
end
%Decision regions comparison
Sr = Sr_norm*sqrt(Eavg);
DM = [-(M-2):2:(M-2)]; % Decision Margins
for i = 1:length(SNR)
So_I(find(real(Sr(i,:)) < DM(1))) = alphabet(1);
if (length(DM) > 1)
for k = 2:length(DM)
So_I(find((real(Sr(i,:)) > DM(k-1))&(real(Sr(i,:)) < DM(k)))) = alphabet(k);
end
end
So_I(find(real(Sr(i,:)) > DM(length(DM)))) = alphabet(length(alphabet));
So(i,:) = So_I;
Pe4(i) = symerr(St,So(i,:))/N;
Pe4_analytic(i) = (2*(M-1)/M)*qfunc(sqrt((6*log2(M)/(M^2-1))*(Eb_N0(i)))); %analytical result
end
'M=4 Simulation Completed'
%M=2 simulation
M = 2; % Number of Symbols
k = log2(M); % Number of bits per Symbol
alphabet = [-(M-1) : 2 : (M-1)];
Eavg = (1/M)*(sum(alphabet.^2)); % average power of transmitted signal
Eb_N0_dB = SNR; % signal to noise ratio
Eb_N0 = 10.^(Eb_N0_dB/10);
Es_N0 = Eb_N0*k; % symbols energy to noise ratio
% M-ASK transmitted signal
St = randsrc(1,N,alphabet);
St_norm = St/sqrt(Eavg);
% received signal
sigma2 = 1./(Es_N0);
N0 = sigma2;
for i = 1:length(SNR)
n = sqrt(sigma2(i)/2)*(randn(1,length(St))+j*randn(1,length(St)));
Sr_norm(i,:) = St_norm + n;
end
% Optimum Receiver Structure /or Decision regions comparison
% Decision Structure
Sr = Sr_norm*sqrt(Eavg); % deNormalization of received signal
DM = [-(M-2):2:(M-2)]; % Decision Margins
for i = 1:length(SNR)
So_I(find(real(Sr(i,:)) < DM(1))) = alphabet(1);
if (length(DM) > 1)
for k = 2:length(DM)
So_I(find((real(Sr(i,:)) > DM(k-1))&(real(Sr(i,:)) < DM(k)))) = alphabet(k);
end
end
So_I(find(real(Sr(i,:)) > DM(length(DM)))) = alphabet(length(alphabet));
So(i,:) = So_I;
Pe2(i) = symerr(St,So(i,:))/N;
Pe2_analytic(i) = (2*(M-1)/M)*qfunc(sqrt((6*log2(M)/(M^2-1))*(Eb_N0(i)))); %analytical result
end
'M=2 Simulation Completed'
semilogy(SNR, Pe16, '-*', SNR, Pe16_analytic,'-o',SNR, Pe8, '-*', SNR, Pe8_analytic,'-o',SNR, Pe4, '-*', SNR, Pe4_analytic,'-o',SNR, Pe2, '-*', SNR, Pe2_analytic,'-o');
ylim([1e-6 1e-1]);
title('Pe: M-ASK, Analytical and Symulation result')
legend('M=16(simulation)', 'M=16(Analytic)','M=8(simulation)', 'M=8(Analytic)','M=4(simulation)', 'M=4(Analytic)','M=2(simulation)', 'M=2(Analytic)');
xlabel('SNR [dB]');
ylabel('Symbol error rate');