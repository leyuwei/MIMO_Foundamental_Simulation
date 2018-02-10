%% ������ʼ��
clc;clear all; close all;
nd = 2^16; %���÷��澫��
realnd = nd/4;
ebnodb = 0:16;
maxpower=1;
ber = zeros(1,length(ebnodb));
theory_rayleigh = zeros(1,length(ebnodb));
ber_EGC = zeros(1,length(ebnodb));
h = waitbar(1,'��Ⱥ����ڼ���...');

%% ��ʼ����
for snr_num=1:length(ebnodb)
    SNR=exp(ebnodb(snr_num)*log(10)/10);
    noe = 0;noeegc = 0; % ������
    nod = 0;nodegc = 0; % ���������
    %16QAM�ź�����
    data1=rand(1,nd)>0.5;
    data2=qam16_modulation(data1);
    data2=data2./maxpower; %������һ��
    
    %ģ�� 1Tx 2Rx
    sigma=calNoisePower(16,data2,ebnodb(snr_num));
    %��1��
    n = randn(1,realnd) + 1j*randn(1,realnd) ;
    h1 = 1/sqrt(2)*(randn(1,realnd) + 1j*randn(1,realnd)); % ���������ŵ�
    data41=data2.*h1+sigma.*n; %���澭��AWGN������˥���ŵ�
    h11=conj(h1)./realnd; %�����ŵ�����ָ���ĸ���������õ���MRC��
    h11egc=conj(h1)./(abs(h1));%�����ŵ�����ָ���ĸ����EGCӦ�õ��abs(h1)
    data411 = data41.*h11; %����MRC��Ϻ�ļ�ֵ
    data411egc = data41.*h11egc; %����EGC��Ϻ�ļ�ֵ
    %��2��
    n = randn(1,realnd) + 1j*randn(1,realnd) ;
    h2 = 1/sqrt(2)*(randn(1,realnd) + 1j*randn(1,realnd)); 
    data42=data2.*h2+sigma.*n;
    h22=conj(h2)./realnd;
    h22egc=conj(h2)./(abs(h2));
    data422 =data42.*h22;
    data422egc = data42.*h22egc; 
    data4=data411+data422; %����������ص��ŵ��µ��źŽ���MRC�ϲ�
    data4egc=data411egc+data422egc;

    % 16QAM ��� for MRC
    [demodata1,uu]=qam16_demodulation(data4,maxpower/2);
    noe2=sum(abs(data1-demodata1));
    nod2=length(data1);
    noe=noe+noe2;
    nod=nod+nod2;
    % 16QAM ��� for EGC
    [demodata1egc,uu]=qam16_demodulation(data4egc,maxpower/2);
    noe2egc=sum(abs(data1-demodata1egc));
    nod2egc=length(data1);
    noeegc=noeegc+noe2egc;
    nodegc=nodegc+nod2egc;
        
    % �����ʼ�¼
    ber(snr_num) = noe/nod;
    ber_EGC(snr_num) = noeegc/nodegc;
    theory_rayleigh(snr_num) = berfading(ebnodb(snr_num),'qam',16,2);
end 

%% �����ͼ
close(h);
figure;
semilogy(ebnodb,ber, 'r*','LineWidth',2 );hold on;
semilogy(ebnodb,ber_EGC, 'b+','LineWidth',2 );hold on;
semilogy(ebnodb,0.5*erfc(sqrt(2*10.^(ebnodb/10))/sqrt(2)), 'b+-' );hold on;
semilogy(ebnodb,theory_rayleigh, 'k-');
axis normal tight;ylabel( 'BER' );xlabel( 'Eb/N0 (dB)' );
title('16QAM��L=2��SIMOϵͳ�µķ���');
legend1=legend( 'MRC�ϲ�����ֵ','EGC�ϲ�����ֵ','AWGN����ֵ' ,'˫����������ֵ');
set(legend1,...
    'Position',[0.147976807959614 0.134315949800738 0.248214289460863 0.165476194109236]);