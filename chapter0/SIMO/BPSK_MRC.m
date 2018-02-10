%% ������ʼ��
clc;clear all; close all;
nd = 2^16; %���÷��澫��
snr_in_dB = 0:16 ;
ber = zeros(1,length(snr_in_dB));
theory_rayleigh = zeros(1,length(snr_in_dB));
ber_EGC = zeros(1,length(snr_in_dB));
h = waitbar(1,'��Ⱥ����ڼ���...');

%% ��ʼ����
for snr_num=1:length(snr_in_dB)
    SNR=exp(snr_in_dB(snr_num)*log(10)/10);
    noe = 0;noeegc = 0; % ������
    nod = 0;nodegc = 0; % ���������
    %BPSK�ź�����
    data1=rand(1,nd)>0.5;
    data2=2.*data1-1;
    
    %ģ�� 1Tx 2Rx
    sigma=calNoisePower(2,data2,snr_in_dB(snr_num));
    %��1��
    n = randn(1,nd) + 1j*randn(1,nd) ;
    h1 = 1/sqrt(2)*(randn(1,nd) + 1j*randn(1,nd)); % ���������ŵ�
    data41=data2.*h1+sigma.*n; %���澭��AWGN������˥���ŵ�
    h11=conj(h1); %�����ŵ�����ָ���ĸ���������õ���MRC��
    h11egc=conj(h1)./abs(h1);%�����ŵ�����ָ���ĸ����EGCӦ�õ��abs(h1)
    data411 = data41.*h11; %����MRC��Ϻ�ļ�ֵ
    data411egc = data41.*h11egc; %����EGC��Ϻ�ļ�ֵ
    %��2��
    n = randn(1,nd) + 1j*randn(1,nd) ;
    h2 = 1/sqrt(2)*(randn(1,nd) + 1j*randn(1,nd)); 
    data42=data2.*h2+sigma.*n;
    h22=conj(h2);
    h22egc=conj(h2)./abs(h2);
    data422 =data42.*h22;
    data422egc = data42.*h22egc; 
    data4=data411+data422; %����������ص��ŵ��µ��źŽ���MRC�ϲ�
    data4egc=data411egc+data422egc;

    % BPSK ��� for MRC
    demodata1=data4 > 0;
    noe2=sum(abs(data1-demodata1));
    nod2=length(data1);
    noe=noe+noe2;
    nod=nod+nod2;
    % BPSK ��� for EGC
    demodata1egc=data4egc > 0;
    noe2egc=sum(abs(data1-demodata1egc));
    nod2egc=length(data1);
    noeegc=noeegc+noe2egc;
    nodegc=nodegc+nod2egc;
        
    % �����ʼ�¼
    ber(snr_num) = noe/nod;
    ber_EGC(snr_num) = noeegc/nodegc;
    theory_rayleigh(snr_num) = berfading(snr_in_dB(snr_num),'psk',2,2);
end 

%% �����ͼ
close(h);
figure;
semilogy(snr_in_dB,ber, 'r*','LineWidth',2 );hold on;
semilogy(snr_in_dB,ber_EGC, 'b+','LineWidth',2 );hold on;
semilogy(snr_in_dB,0.5*erfc(sqrt(2*10.^(snr_in_dB/10))/sqrt(2)), 'b+-' );hold on;
semilogy(snr_in_dB,theory_rayleigh, 'k-');
axis normal tight;ylabel( 'BER' );xlabel( 'Eb/N0 (dB)' );
title('BPSK��L=2��SIMOϵͳ�µķ���');
legend1=legend( 'MRC�ϲ�����ֵ','EGC�ϲ�����ֵ','AWGN����ֵ' ,'˫����������ֵ');
set(legend1,...
    'Position',[0.147976807959614 0.134315949800738 0.248214289460863 0.165476194109236]);