function BER=BPSK_mode(N0,Eb) %% inputs are energy per bit and the psd of the noise 
%% setting parameters 

Tb=5; %% bit duration in secs  
%Eb=5
N_bit = 5000; %%number of samples per bit 
t_bit=linspace(0,Tb,N_bit); %% time base for each bit  
msg_l = 5000; %% number of bits sent 
t_signal = linspace(0,msg_l*Tb,msg_l*N_bit);%% total duration of the messag 
fc=2/Tb; %% frequency of the carrier  
%N0=5

%% message source ------> generates a rondom stream of 0s and 1s

message= randi([0 1],1,msg_l);

%% signal transimission encoder --------> polar non return tozero: (1)->1,(0)->-1

encodedMessage=[];

for i=1:1:msg_l
    if message(i)==1
        signal_seg=sqrt(Eb)*ones(1,N_bit);
    elseif message(i) == 0
        signal_seg=-sqrt(Eb)*ones(1,N_bit);
    end
    encodedMessage=[encodedMessage signal_seg];
end 

% plot(encodedMessage)

%% BPSK modulator---> multiply with the carrier 

carrier=sqrt(2/Tb)*cos(2*pi*fc*t_signal); %% normalize basis function 

modulatedSignal=encodedMessage.*carrier;

%plot(t_signal,modulatedSignal)

%% constellation of transsmitted BPSK

basis_func=sqrt(2/Tb)*cos(2*pi*fc*t_bit);
si_vector=[];

for i=1:N_bit:length(modulatedSignal)
    vec=modulatedSignal(i:i+N_bit-1);
    vec=vec.*basis_func;
    intg=trapz(t_bit,vec); %% seperation is tb 
    si_vector=[si_vector intg];
end  


%scatterplot(si_vector)
%% calculate the minimum distance between constellation points 

distance=[];

for i=1:1:length(si_vector)
    d=sqrt(si_vector(i)^2);% calculate the distnace of each point to the origin
    distance =[distance d];
end

[dmin,index]=min(distance);
distance(index)=Inf;
dmin=dmin+min(distance);

%% adding whie additive Gussian noise

%recievedSignal=awgn(modulatedSignal,0.1,'measured'); %% 10 snr in dB is used 
recievedSignal = modulatedSignal+normrnd(0,N0/2,1,length(modulatedSignal));

%% demodulation using correlator and constellation of recieved BPSK

basis_func=sqrt(2/Tb)*cos(2*pi*fc*t_bit);
xi_vector=[];

for i=1:N_bit:length(recievedSignal)
    vec=recievedSignal(i:i+N_bit-1);
    vec=vec.*basis_func;
    intg=trapz(t_bit,vec); %% seperation is tb 
    xi_vector=[xi_vector intg];
end  


%scatterplot(xi_vector)



%% signal transimission decoder -----> using ML rule 

rec_signal_Decoded = [];
d1=[];
d2=[];

for i=1:1:length(xi_vector)    
x1 = sqrt((xi_vector(i)-sqrt(Eb))^2);
x2 = sqrt((xi_vector(i)+sqrt(Eb))^2);
d1=[d1 x1];
d2=[d2 x2];
end
d1=abs(d1);
d2=abs(d2);

for i=1:1:length(xi_vector)
    if d1(i)<d2(i)
        rec_signal_Decoded = [rec_signal_Decoded 1];
    else
        rec_signal_Decoded = [rec_signal_Decoded 0];
    end
    
end


y = dmin/(2*sqrt(N0));
BER = 0.5*erfc(y);

end 