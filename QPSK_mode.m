function BER=QPSK_mode(N0,Eb)
%% setting parameters 

Tb = 5; %% bit duration in secs 
T=2*Tb;  %%dibit duration 
%Eb = 5 %% enerdy per bit
E = 2*Eb; %%energy per symbol

N_bit = 1000; %%number of samples per bit 
N_dibit = 2 * N_bit;

t_bit=linspace(0,Tb,N_bit); %% time base for each bit  
t_dibit = linspace(0,T,N_dibit);

msg_l = 100; %% number of bits sent which has to be even 
t_signal = linspace(0,msg_l*Tb,msg_l*N_bit);%% total duration of the messag 

fc=2/Tb; %% frequency of the carrier 
%N0=5
%% message source ------> generates a rondom stream of 0s and 1s

message= randi([0 1],1,msg_l);
odd_bits = [];
even_bits = [];

for i=1:1:msg_l
    if mod(i,2)== 1
        odd_bits=[odd_bits message(i)];
    else
        even_bits=[even_bits message(i)];
    end
end

%% signal transimission encoder --------> polar non return tozero: (1)->1,(0)->-1

encodedodd=[];
encodedeven=[];

for i=1:1:length(odd_bits)
    if odd_bits(i)==1
        signal_seg=sqrt(E)*ones(1,2*N_bit);
    elseif odd_bits(i) == 0
        signal_seg=-sqrt(E)*ones(1,2*N_bit);
    end
    encodedodd=[encodedodd signal_seg];
end 

for i=1:1:length(even_bits)
    if even_bits(i)==1
        signal_seg=sqrt(E)*ones(1,2*N_bit);
    elseif even_bits(i) == 0
        signal_seg=-sqrt(E)*ones(1,2*N_bit);
    end
    encodedeven=[encodedeven signal_seg];
end 

%% QPSK modulator---> multiply with the carrier 

carrier_i=sqrt(2/T)*cos(2*pi*fc*t_signal); %% normalize basis function 
carrier_Q=sqrt(2/T)*sin(2*pi*fc*t_signal); %% normalize basis function 

modulatedSignal=encodedodd.*carrier_i + encodedeven.*carrier_Q;

%plot(t_signal,modulatedSignal)

%% constellation of transsmitted QPSK

basis_func_i=sqrt(2/T)*cos(2*pi*fc*t_dibit);
basis_func_Q=sqrt(2/T)*sin(2*pi*fc*t_dibit);
si1_vector=[];
si2_vector=[];

for i=1:N_dibit:length(modulatedSignal)
    vec=modulatedSignal(i:i+N_dibit-1);
    vec=vec.*basis_func_i;
    intg=trapz(t_dibit,vec); %% seperation is tb 
    si1_vector=[si1_vector intg];
end  

for i=1:N_dibit:length(modulatedSignal)
    vec=modulatedSignal(i:i+N_dibit-1);
    vec=vec.*basis_func_Q;
    intg=trapz(t_dibit,vec); %% seperation is tb 
    si2_vector=[si2_vector intg];
end 

si_vector=[si1_vector ; si2_vector];
si_vector = si_vector/2; %% cause we have integrted over 2Tb not one 
%scatterplot(transpose(si_vector))


%% calculate the minimum distance between constellation points 

distance=[];
[R,L] = size(si_vector);

for i=1:1:L
    d=sqrt(si_vector(1,i)^2 + si_vector(2,i)^2);% calculate the distnace of each point to the origin
    distance =[distance d];
end

[dmin,index]=min(distance);
distance(index)=Inf;
dmin=dmin+min(distance);

%% adding whie additive Gussian noise

%recievedSignal=awgn(modulatedSignal,0.1,'measured'); %% 10 snr in dB is used 
recievedSignal = modulatedSignal+unifrnd(0,N0/2,1,length(modulatedSignal));

%% constellation of recieved QPSK

basis_func_i=sqrt(2/T)*cos(2*pi*fc*t_dibit);
basis_func_Q=sqrt(2/T)*sin(2*pi*fc*t_dibit);
xi1_vector=[];
xi2_vector=[];

for i=1:N_dibit:length(recievedSignal)
    vec=recievedSignal(i:i+N_dibit-1);
    vec=vec.*basis_func_i;
    intg=trapz(t_dibit,vec); %% seperation is tb 
    xi1_vector=[xi1_vector intg];
end  

for i=1:N_dibit:length(recievedSignal)
    vec=recievedSignal(i:i+N_dibit-1);
    vec=vec.*basis_func_Q;
    intg=trapz(t_dibit,vec); %% seperation is tb 
    xi2_vector=[xi2_vector intg];
end 

xi_vector=[xi1_vector ; xi2_vector];
xi_vector=xi_vector/2; %% cause we have interated over 2 Tb not one 
%scatterplot(transpose(xi_vector))

%% signal transimission decoder -----> using ML rule

decode_odd = [];

d1=[];
d2=[];

for i=1:1:length(xi1_vector)    
x1 = sqrt((xi1_vector(i)-sqrt(Eb))^2);
x2 = sqrt((xi1_vector(i)+sqrt(Eb))^2);
d1=[d1 x1];
d2=[d2 x2];
end
d1=abs(d1);
d2=abs(d2);

for i=1:1:length(xi2_vector)
    if d1(i)<d2(i)
        decode_odd = [decode_odd 1 0];
    else
        decode_odd = [decode_odd 0 0];
    end
    
end

decode_even = [];

d1=[];
d2=[];

for i=1:1:length(xi2_vector)    
x1 = sqrt((xi2_vector(i)-sqrt(Eb))^2);
x2 = sqrt((xi2_vector(i)+sqrt(Eb))^2);
d1=[d1 x1];
d2=[d2 x2];
end
d1=abs(d1);
d2=abs(d2);

for i=1:1:length(xi2_vector)
    if d1(i)<d2(i)
        decode_even = [decode_even  0 1];
    else
        decode_even = [decode_even  0 0];
    end
    
end


rec_signal_Decoded = decode_even + decode_odd;


y = dmin/(2*sqrt(N0));
BER = 0.5*erfc(y);

end 