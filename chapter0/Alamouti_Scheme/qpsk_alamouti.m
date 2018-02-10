% qpsk BER�������
% AWGN��Rayleigh˥���ŵ�MISOϵͳ���棨Alamouti�ԱȰ汾��

%% ����׼��
clear all; close all; clc; 
N = 2^16;
m = N/2;
nTx = 2;
maxlen = sqrt(2);
ebn0db = -2:16;
Xs = [1+1j -1+1j -1-1j 1-1j]; % �����㣬��ML�о�ʹ��
BER = zeros(1,length(ebn0db));
theoryBer_nRx22 = zeros(1,length(ebn0db));
% ����qpsk����
ip = rand(1,N)>0.5;
ST = qpsk_modulation(ip);
ST=ST./maxlen;
hh = waitbar(0,'���Ե� ���ڼ���...');

%% ��ʼ����Alamouti
for n=1:length(ebn0db)
    waitbar(n/length(ebn0db),hh,['Alamouti���㣺' num2str(100*n/length(ebn0db)) '%']);
    ST_estimated=zeros(1,length(ST));
    for k=1:2:length(ST)
        % ��ԭʼ�źŽ��д����任����������·
        s0=ST(k);
        s1=ST(k+1);
        
        % ��ʼ������
        sigma=calNoisePower(4,ST,ebn0db(n));
        Noise0=sigma.*(randn+1j*randn)./sqrt(2);
        Noise1=sigma.*(randn+1j*randn)./sqrt(2);
        
        % ��ʼ������iid������˥���ŵ�ϵ��
        h0=(randn+1j*randn)*sqrt(1/2);
        h1=(randn+1j*randn)*sqrt(1/2);
        hEq=(abs(h0))^2+(abs(h1))^2;

        % ����Alamoutiԭʼ���� ʽ(11) �����źű�ʾ����
        r0=h0*s0+h1*s1+Noise0;  
        r1=-h0*conj(s1)+h1*conj(s0)+Noise1;      

        % ����Alamoutiԭʼ���� ʽ(12) �ϲ�ģʽ��ʾ����
        S0=conj(h0).*r0+h1.*conj(r1);
        S1=conj(h1).*r0-h0.*conj(r1);
        
        % ����ML׼����о�
        S=[S0,S1]./hEq;
        % ����ML׼���Ҫ�󣬼���Squared Euclidean Distance��ƽ��ŷ�Ͼ���
        % ����˵�����б��������һ����������
        d11=squaredis(S(1),[1 1],maxlen);
        d12=squaredis(S(1),[-1 1],maxlen);
        d13=squaredis(S(1),[-1 -1],maxlen);
        d14=squaredis(S(1),[1 -1],maxlen);
        d21=squaredis(S(2),[1 1],maxlen);
        d22=squaredis(S(2),[-1 1],maxlen);
        d23=squaredis(S(2),[-1 -1],maxlen);
        d24=squaredis(S(2),[1 -1],maxlen);
        % ����ML׼����о�Ҫ���ҳ���С����
        [ds1_min,position1]=min([d11,d12,d13,d14]);
        [ds2_min,position2]=min([d21,d22,d23,d24]);
        % ����
        s1_estemated=Xs(position1);         
        s2_estemated=Xs(position2);
        ST_estimated(k)=s1_estemated;
        ST_estimated(k+1)=s2_estemated; 
    end
    % ���
    [R,~]=qpsk_demodulation(ST_estimated);
    % ����BER
    no_errors=size(find(ip-R),2);
    BER(1,n)=no_errors/length(ip);
    theoryBer_nRx22(n) = berfading(ebn0db(n),'psk',4,2);
end

%% ��ʼ����Beamforming
nErr2 = zeros(1,length(ebn0db));
for ii = 1:length(ebn0db)
   ip = rand(1,N)>0.5;
   s = qpsk_modulation(ip); % ʹ�ñ���������qpsk����ͼ
   s = s./maxlen; % ������һ����ֻ�Ķ���������
   n = 1/sqrt(2)*(randn(1,m) + 1j*randn(1,m)); % ��˹������
   h = 1/sqrt(2)*(randn(nTx,m) + 1j*randn(nTx,m)); % �����ŵ�
   hEff = h.*exp(-1j*angle(h)); % ���Ͷ˲�������
   sigma=calNoisePower(4,s,ebn0db(ii)); % ��������
   sr = (1/sqrt(nTx))*kron(ones(nTx,1),s); %�����������ߵķ������� ��Ϊһ��2*N�ľ��� Ҳ����˵��һ����ά����
   y2 = sum(hEff.*sr,1) + sigma*n; 
   y2Hat = maxlen.*y2./sum(hEff,1); 
   [ip2Hat,uu2] = qpsk_demodulation(y2Hat);
   nErr2(ii) = size(find(ip- ip2Hat),2);
end
simBer2 = nErr2/N; % ����BER ����Beamforming��

%% ��ͼ
close(hh);figure;
semilogy(ebn0db,theoryBer_nRx22,'k*-','LineWidth',1); hold on;
semilogy(ebn0db,simBer2,'bx','LineWidth',2); hold on;
semilogy(ebn0db,BER,'r*','LineWidth',2);
legend('L=2 ������','2��1�� Beamforming����','2��1�� Alamouti����');
title('QPSK��AWGN������˥���ŵ��µ�MISOϵͳ���棨Alamouti��');
axis normal tight; grid on;
xlabel('Eb/N0 (dB)'); ylabel('BER');