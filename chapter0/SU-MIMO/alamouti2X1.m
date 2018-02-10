clc;clear all;close all;

%% ��ʼ����Alamouti
ebn0db=-2:1:16;
precision=2^16;
sim=zeros(1,length(ebn0db));
theo=zeros(1,length(ebn0db));
for ii=1:length(ebn0db)
    sim(ii)=alamouti_sim(precision,16,ebn0db(ii));
    theo(ii)=berfading(ebn0db(ii),'qam',16,2);
end

%% ��ʼ����Beamforming
maxlen=sqrt(10);
N=precision;
m=N/log2(16);
nTx=2;
nErr2 = zeros(1,length(ebn0db));
for ii = 1:length(ebn0db)
   ip = rand(1,N)>0.5;
   s = qam16_modulation(ip); % ʹ�ñ���������qam16����ͼ
   s = s./maxlen; % ������һ����ֻ�Ķ���������
   n = 1/sqrt(2)*(randn(1,m) + 1j*randn(1,m)); % ��˹������
   h = 1/sqrt(2)*(randn(nTx,m) + 1j*randn(nTx,m)); % �����ŵ�
   hEff = h.*exp(-1j*angle(h)); % ���Ͷ˲�������
   sigma=calnoisesigma(16,s,ebn0db(ii),1); % ��������
   sr = kron(ones(nTx,1),s)./sqrt(2); %�����������ߵķ������� ��Ϊһ��2*N�ľ��� Ҳ����˵��һ����ά����
   y2 = sum(hEff.*sr,1) + sigma*n; 
   y2Hat = maxlen.*y2./sum(hEff,1); 
   [ip2Hat,uu2] = qam16_demodulation(y2Hat,maxlen);
   nErr2(ii) = size(find(ip- ip2Hat),2);
end
simBer2 = nErr2/N; % ����BER ����Beamforming��

%% ��ͼ
semilogy(ebn0db,theo,'k+-','LineWidth',1);hold on;
semilogy(ebn0db,sim,'r*','LineWidth',2);hold on;
semilogy(ebn0db,simBer2,'b+','LineWidth',2);
axis normal tight; grid on;
xlabel('Eb/N0 (dB)'); ylabel('BER');
title('16QAM��Alamouti����');
legend('˫������ֵ','Alamouti�����','Beamforming�����');