%% BERs plots for BPSK and QPSK both theortical and through simulations

clear
clc

N0=10;
Eb=[];
x=0.0001;
Eb=[Eb x];
for i=1:1:5 %generating a vector of Eb values ranging from 0.0001 to 
    x=x*10;
    for j=1:1:9
        Eb=[Eb j*x];
    end
end

%% BER Theortical

dmin_BPSK=2*sqrt(Eb); % minimum distance between any two points in the constellationof BPSK  
dmin_QPSK=2*sqrt(Eb); % minimum distance between any two points in the constellationof QPSK

y=dmin_BPSK/(2*sqrt(N0));
BER_BPSK=0.5*erfc(y); % calculate the BER using the error function 
y=dmin_QPSK/(2*sqrt(N0));
BER_QPSK=0.5*erfc(y);

subplot(2,1,1)
semilogy(20*log10(Eb/N0),BER_BPSK)
hold on 
semilogy(20*log10(Eb/N0),BER_QPSK)
title('BERs of BPSK and QPSK Theortical')
xlabel('Eb/N0(dB)')
ylabel('Bit Error Rate')
legend('BPSK','QPSK')
grid on 
xlim([-100 50])
%% BER through Simulation 

BER_BPSK_sim=[];
BER_QPSK_sim=[];

for i =1:1:length(Eb) % run the BPSK modulation for different valus of Eb/N0
    b=BPSK_mode(N0,Eb(i));
    BER_BPSK_sim = [BER_BPSK_sim b];
end

subplot(2,1,2)
semilogy(20*log10(Eb/N0),BER_BPSK_sim)
hold on 

for i =1:1:length(Eb)% run the QPSK modulation for different valus of Eb/N0
    b=QPSK_mode(N0,Eb(i));
    BER_QPSK_sim=[BER_QPSK_sim b];
end

semilogy(20*log10(Eb/N0),BER_QPSK_sim)

title('BERs of BPSK and QPSK Simulation')
xlabel('Eb/N0(dB)')
ylabel('Bit Error Rate')
legend('BPSK','QPSK')
grid on 
xlim([-100 50])

