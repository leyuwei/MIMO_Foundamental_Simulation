% qam16 BER�������
% AWGN��Rayleigh˥���ŵ�MISOϵͳ���棨Alamouti�ԱȰ汾��

%% ����׼��
clear all; close all; clc; 
N = 2^16;
m = N/4;
nTx = 2;
maxlen = sqrt(10);
ebn0db = -2:16;
Xs = [1+1j 3+1j 3+3j 1+3j -1+3j -1+1j -3+3j -3+1j -3-1j -3-3j -1-3j -1-1j 1-1j 3-1j 1-3j 3-3j]; % �����㣬��ML�о�ʹ��
BER = zeros(1,length(ebn0db));
theoryBer_nRx22 = zeros(1,length(ebn0db));
% ����qam16����
ip = rand(1,N)>0.5;
ST = qam16_modulation(ip);
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
        sigma=calNoisePower(16,ST,ebn0db(n));
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
        Xss=[real(Xs);imag(Xs)];
        d11=squaredis(S(1),[Xss(1,1) Xss(2,1)],maxlen);
        d12=squaredis(S(1),[Xss(1,2) Xss(2,2)],maxlen);
        d13=squaredis(S(1),[Xss(1,3) Xss(2,3)],maxlen);
        d14=squaredis(S(1),[Xss(1,4) Xss(2,4)],maxlen);
        d15=squaredis(S(1),[Xss(1,5) Xss(2,5)],maxlen);
        d16=squaredis(S(1),[Xss(1,6) Xss(2,6)],maxlen);
        d17=squaredis(S(1),[Xss(1,7) Xss(2,7)],maxlen);
        d18=squaredis(S(1),[Xss(1,8) Xss(2,8)],maxlen);
        d19=squaredis(S(1),[Xss(1,9) Xss(2,9)],maxlen);
        d110=squaredis(S(1),[Xss(1,10) Xss(2,10)],maxlen);
        d111=squaredis(S(1),[Xss(1,11) Xss(2,11)],maxlen);
        d112=squaredis(S(1),[Xss(1,12) Xss(2,12)],maxlen);
        d113=squaredis(S(1),[Xss(1,13) Xss(2,13)],maxlen);
        d114=squaredis(S(1),[Xss(1,14) Xss(2,14)],maxlen);
        d115=squaredis(S(1),[Xss(1,15) Xss(2,15)],maxlen);
        d116=squaredis(S(1),[Xss(1,16) Xss(2,16)],maxlen);
        d21=squaredis(S(2),[Xss(1,1) Xss(2,1)],maxlen);
        d22=squaredis(S(2),[Xss(1,2) Xss(2,2)],maxlen);
        d23=squaredis(S(2),[Xss(1,3) Xss(2,3)],maxlen);
        d24=squaredis(S(2),[Xss(1,4) Xss(2,4)],maxlen);
        d25=squaredis(S(2),[Xss(1,5) Xss(2,5)],maxlen);
        d26=squaredis(S(2),[Xss(1,6) Xss(2,6)],maxlen);
        d27=squaredis(S(2),[Xss(1,7) Xss(2,7)],maxlen);
        d28=squaredis(S(2),[Xss(1,8) Xss(2,8)],maxlen);
        d29=squaredis(S(2),[Xss(1,9) Xss(2,9)],maxlen);
        d210=squaredis(S(2),[Xss(1,10) Xss(2,10)],maxlen);
        d211=squaredis(S(2),[Xss(1,11) Xss(2,11)],maxlen);
        d212=squaredis(S(2),[Xss(1,12) Xss(2,12)],maxlen);
        d213=squaredis(S(2),[Xss(1,13) Xss(2,13)],maxlen);
        d214=squaredis(S(2),[Xss(1,14) Xss(2,14)],maxlen);
        d215=squaredis(S(2),[Xss(1,15) Xss(2,15)],maxlen);
        d216=squaredis(S(2),[Xss(1,16) Xss(2,16)],maxlen);
        % ����ML׼����о�Ҫ���ҳ���С����
        [ds1_min,position1]=min([d11,d12,d13,d14,d15,d16,d17,d18,d19,d110,d111,d112,d113,d114,d115,d116]);
        [ds2_min,position2]=min([d21,d22,d23,d24,d25,d26,d27,d28,d29,d210,d211,d212,d213,d214,d215,d216]);
        % ����
        s1_estemated=Xs(position1);         
        s2_estemated=Xs(position2);
        ST_estimated(k)=s1_estemated;
        ST_estimated(k+1)=s2_estemated; 
    end
    % ���
    [R,~]=qam16_demodulation(ST_estimated,maxlen);
    % ����BER
    no_errors=size(find(ip-R),2);
    BER(1,n)=no_errors/length(ip);
    theoryBer_nRx22(n) = berfading(ebn0db(n),'qam',16,2);
end

%% ��ʼ����Beamforming
nErr2 = zeros(1,length(ebn0db));
for ii = 1:length(ebn0db)
   ip = rand(1,N)>0.5;
   s = qam16_modulation(ip); % ʹ�ñ���������qam16����ͼ
   s = s./maxlen; % ������һ����ֻ�Ķ���������
   n = 1/sqrt(2)*(randn(1,m) + 1j*randn(1,m)); % ��˹������
   h = 1/sqrt(2)*(randn(nTx,m) + 1j*randn(nTx,m)); % �����ŵ�
   hEff = h.*exp(-1j*angle(h)); % ���Ͷ˲�������
   sigma=calNoisePower(16,s,ebn0db(ii)); % ��������
   sr = (1/sqrt(nTx))*kron(ones(nTx,1),s); %�����������ߵķ������� ��Ϊһ��2*N�ľ��� Ҳ����˵��һ����ά����
   y2 = sum(hEff.*sr,1) + sigma*n; 
   y2Hat = maxlen.*y2./sum(hEff,1); 
   [ip2Hat,uu2] = qam16_demodulation(y2Hat,maxlen);
   nErr2(ii) = size(find(ip- ip2Hat),2);
end
simBer2 = nErr2/N; % ����BER ����Beamforming��

%% ��ͼ
close(hh);figure;
semilogy(ebn0db,theoryBer_nRx22,'k*-','LineWidth',1); hold on;
semilogy(ebn0db,simBer2,'bx','LineWidth',2); hold on;
semilogy(ebn0db,BER,'r*','LineWidth',2);
legend('L=2 ������','2��1�� Beamforming����','2��1�� Alamouti����');
title('16QAM��AWGN������˥���ŵ��µ�MISOϵͳ���棨Alamouti��');
axis normal tight; grid on;
xlabel('Eb/N0 (dB)'); ylabel('BER');