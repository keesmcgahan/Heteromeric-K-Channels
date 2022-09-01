clear all
clc

V=[-100:.01:100];
%%%%%%%%%%%%%
%terms that are probv# are the data probability curve pulled frm the paper
%terms that are prob4v# are the power four probability curve fir to data

%%%%% %open as function of a
syms x


%%%%%%%best fit
%   prob4vrm32x=1./(1+exp(-(x+59.18)/15.61));
%   prob4vrm1x=1./(1+exp(-(x+44.08)/24.75));

%%%%Furthest fitting
%prob4vrm1x=1./(1+exp(-(x+37.38)/20.75));
prob4vrm1x=1./(1+exp(-(x+50.73)/28.75));
%prob4vrm32x=1./(1+exp(-(x+65.83)/19.61));
prob4vrm32x=1./(1+exp(-(x+54.14)/12.61));

S25=zeros(1,41);
S50=zeros(1,41);
S75=zeros(1,41);
S10=zeros(1,41);
S35=zeros(1,41);
S85=zeros(1,41);

a=zeros(1,41);
for i=1:41
    a(i)=(i-1)/10;
    eq25=prob4vrm1x.^(4-a(i)).*prob4vrm32x.^a(i)==.25;
    eq50= prob4vrm1x.^(4-a(i)).*prob4vrm32x.^a(i)==.5;
    eq75=prob4vrm1x.^(4-a(i)).*prob4vrm32x.^a(i)==.75;
    eq10=prob4vrm1x.^(4-a(i)).*prob4vrm32x.^a(i)==.1;
    eq35=prob4vrm1x.^(4-a(i)).*prob4vrm32x.^a(i)==.35;
    eq85=prob4vrm1x.^(4-a(i)).*prob4vrm32x.^a(i)==.85;
    S25(i)=single(vpasolve(eq25));
    S50(i)=single(vpasolve(eq50));
    S75(i)=single(vpasolve(eq75));
    S10(i)=single(vpasolve(eq10));
    S35(i)=single(vpasolve(eq35));
    S85(i)=single(vpasolve(eq85));
end    
a4=[4,3,1,0];
prob50data=[-32,-20,-2,-1];
prob25data=[-45.18,-37.57,-21.87,-20.67];
prob75data=[-18.8,-2.42,16.67,19.87];
prob10data=[-58.3,-55.15,-39.4,-42.74];
prob35data=[-39.4,-29.9,-12.5,-12.76];
prob85data=[-11.18,7.75,27.4,31.9];

figure(1)
hold on
plot(a,S10)
plot(a,S25)
plot(a,S35)
plot(a,S50)
plot(a,S75)
plot(a,S85)
plot(a4,prob10data,'bp')
plot(a4,prob25data,'rd')
plot(a4,prob35data,'yx')
plot(a4,prob50data,'mo')
plot(a4,prob75data,'g^')
plot(a4,prob85data,'cs')
hold off
title('%open as function of a','fontsize',20)
legend('10% open','25% open','35% open','50% open','75% open','85% open')
xlabel('a')
ylabel('Voltage where Probability is achieved')


%%%%%%%%%%%%%%%%%%%%LOOKING AT IT WITH THE DATA FROM PAPER 1.1/1.2%%%%%
probvrm32=1./(1+exp(-(V+32)/12));
probvrm20=1./(1+exp(-(V+20)/16));
probvrm18=1./(1+exp(-(V+18)/12));
probvrm2=1./(1+exp(-(V+2)/17));
probvrm1=1./(1+exp(-(V+1)/19));

%%%%%%%%%%%%Optimally Fit set for kv1.1 data to Al-Sabi data
 % prob4vrm32=1./(1+exp(-(V+59.18)/15.61));
 % prob4vrm1=1./(1+exp(-(V+44.08)/24.75));




%%%%%%%%%%%FURTHEST FITTING K VALUES THAT STILL GIVE SAME MAG ERR
%prob4vrm1=1./(1+exp(-(V+37.38)/20.75));
prob4vrm1=1./(1+exp(-(V+50.73)/28.75));
%prob4vrm32=1./(1+exp(-(V+65.83)/19.61));
prob4vrm32=1./(1+exp(-(V+54.14)/12.61));

het321r40=prob4vrm32.^4.*prob4vrm1.^0;
het321r31=prob4vrm32.^3.*prob4vrm1.^1;
het321r22=prob4vrm32.^2.*prob4vrm1.^2;
het321r13=prob4vrm32.^1.*prob4vrm1.^3;
het321r04=prob4vrm32.^0.*prob4vrm1.^4;
ran321cursum5050=.0625*het321r40+.25*het321r31+.375*het321r22+.25*het321r13+.0625*het321r04;
ran321cursum2575=.0039*het321r40+.0469*het321r31+.2109*het321r22+.4219*het321r13+.3164*het321r04;
ran321cursum7525=.0039*het321r04+.0469*het321r13+.2109*het321r22+.4219*het321r31+.3164*het321r40;
hom7525rat=.75*het321r40+.25*het321r04;
hom2575rat=.25*het321r40+.75*het321r04;



figure(2)
hold on
plot(V,probvrm32,'g')
plot(V,prob4vrm32.^4,'g--')
plot(V,probvrm1,'m')
plot(V,prob4vrm1.^4,'m--')
 plot(V,probvrm20,'r')
 plot(V,prob4vrm32.^3.*prob4vrm1.^1,'r--')
 plot(V,probvrm2,'c')
 plot(V,prob4vrm32.^1.*prob4vrm1.^3,'c--')
%plot(V,prob4vrm32.^2.*prob4vrm1.^2)
legend('4:0 data','4:0 model fit','0:4 data','0:4 model fit','3:1 data','3:1 model pred.','1:3','1:3 model pred.')
xlabel('Voltage (mV)')
ylabel('Open Probability')
title('%open as function of a','fontsize',20)
axis([-90 90 0 1])
hold off


% %%%%%%%%%%%%1999 coexpression DADAMO work%%%%%
probvrm28=1./(1+exp(-(V+28.5)/7.8));
prob4vrm28=1./(1+exp(-(V+46.2)/10.1));

probvrm13=1./(1+exp(-(V+13.4)/6.5));
prob4vrm13=1./(1+exp(-(V+28.1)/8.4));

probcoex2813=1./(1+exp(-(V+16.3)/7.9));
probconcat2813=1./(1+exp(-(V+18.3)/8.3));
prob4coex2813=1./(1+exp(-(V+34.9)/10.3));
prob4concat2813=1./(1+exp(-(V+36.8)/10.8));

het2813r40=prob4vrm28.^4.*prob4vrm13.^0;
het2813r31=prob4vrm28.^3.*prob4vrm13.^1;
het2813r22=prob4vrm28.^2.*prob4vrm13.^2;
het2813r13=prob4vrm28.^1.*prob4vrm13.^3;
het2813r04=prob4vrm28.^0.*prob4vrm13.^4;
ran2813cursum=.0625*het2813r40+.25*het2813r31+.375*het2813r22+.25*het2813r13+.0625*het2813r04;
het28135050=.5*het2813r40+.5*het2813r04;


figure(3)
hold on
plot(V,het2813r40)
plot(V,het2813r31,'m')
plot(V,het2813r22)
plot(V,het2813r13)
plot(V,het2813r04)
plot(V,ran2813cursum)
plot(V,prob4coex2813.^4)
plot(V,het28135050,'k')
hold off
legend('4:0','3:1','2:2','1:3','0:4','rand','coexpression data','50:50')
title('dadamo 1.1/1.2','fontsize',20)
xlabel('Voltage (mV)')
ylabel('Open Probability')
axis([-60 20 0 1])




%%%%%%%%%%%%%%%%%%%%%%%IMBRICI 1.1/1.2 coex%%%%%%%%%%%
probvrm26=1./(1+exp(-(V+26.4)/6));
prob4vrm26=1./(1+exp(-(V+40)/7.8));

probvrm8=1./(1+exp(-(V+8.3)/10));
prob4vrm8=1./(1+exp(-(V+30.9)/13));

probcoex268=1./(1+exp(-(V+18)/12));

het268r40=prob4vrm26.^4.*prob4vrm8.^0;
het268r31=prob4vrm26.^3.*prob4vrm8.^1;
het268r22=prob4vrm26.^2.*prob4vrm8.^2;
het268r13=prob4vrm26.^1.*prob4vrm8.^3;
het268r04=prob4vrm26.^0.*prob4vrm8.^4;
ran268cursum=.0625*het268r40+.25*het268r31+.375*het268r22+.25*het268r13+.0625*het268r04;
het2685050=.5*het268r40+.5*het268r04;

figure(4)
hold on
plot(V,het268r40)
plot(V,het268r31,'m')
plot(V,het268r22)
plot(V,het268r13)
plot(V,het268r04)
plot(V,ran268cursum)
plot(V,probcoex268)
plot(V,het2685050,'k')
hold off
legend('4:0','3:1','2:2','1:3','0:4','rand','coexpression data','50:50')
xlabel('Voltage (mV)')
ylabel('Open Probability')
title('imbirici 1.1/1.2')
axis([-60 20 0 1])

%%%%%%%%%%7.2 and Mutant STORY MIcel 2013%%%%%%%%%%%%%%%%%%%%%%%%%%
%probvrm30=1./(1+exp(-(V+30)/10));
%probvrm18=1./(1+exp(-(V+18)/10));
%probvrm5=1./(1+exp(-(V+5)/10));
%probvr5=1./(1+exp(-(V-5)/10));
probvr35=1./(1+exp(-(V-35.4)/23));
probvr45=1./(1+exp(-(V-45)/16));
probvr17=1./(1+exp(-(V-17.7)/19.3));
probvrm22=1./(1+exp(-(V+22.9)/12.3));
probvrm6=1./(1+exp(-(V+6.8)/16.7));

%prob4vrm30=(1./(1+exp(-(V+54)/14.559)));
prob4vr45=(1./(1+exp(-(V-6.9)/22.9)));
prob4vrm22=(1./(1+exp(-(V+51.9)/17.44)));
prob4vrm6=(1./(1+exp(-(V+45.1)/22.9)));
prob4vr35=(1./(1+exp(-(V+16.1)/30.9)));
prob4vr17=(1./(1+exp(-(V+26)/26.27)));

%%%%%%CUR SUM equations 7.2/mutant
Het7q40=prob4vrm22.^4;
Het7q31=prob4vrm22.^3.*prob4vr45.^1;
Het7q22=prob4vrm22.^2.*prob4vr45.^2;
Het7q13=prob4vrm22.^1.*prob4vr45.^3;
Het7q04=prob4vr45.^4;
randcursum2q=.0625*Het7q40+.25*Het7q31+.375*Het7q22+.25*Het7q13+.0625*Het7q04;

Het7w40=prob4vrm22.^4;
Het7w31=prob4vrm22.^3.*prob4vr35.^1;
Het7w22=prob4vrm22.^2.*prob4vr35.^2;
Het7w13=prob4vrm22.^1.*prob4vr35.^3;
Het7w04=prob4vr35.^4;
randcursum2w=.0625*Het7w40+.25*Het7w31+.375*Het7w22+.25*Het7w13+.0625*Het7w04;



figure(5)
hold on
plot(V,prob4vrm22.^4)
plot(V,prob4vrm22.^3.*prob4vr45.^1,'m')
plot(V,prob4vrm22.^2.*prob4vr45.^2)
plot(V,prob4vrm22.^1.*prob4vr45.^3)
plot(V,prob4vr45.^4)
plot(V,randcursum2q)
plot(V,prob4vr17.^4)
plot(V,.5*prob4vrm22.^4+.5*prob4vr45.^4,'k')
title('7.2 and 7.2q 2013','fontsize',20)
legend('4:0','3:1','2:2','1:3','0:4','rand','coexpression data','50:50')
xlabel('Voltage (mV)')
ylabel('Open Probability')
hold off
axis([-60 100 0 1])

figure(6)
hold on
plot(V,prob4vrm22.^4)
plot(V,prob4vrm22.^3.*prob4vr35.^1,'m')
plot(V,prob4vrm22.^2.*prob4vr35.^2)
plot(V,prob4vrm22.^1.*prob4vr35.^3)
plot(V,prob4vr35.^4)
plot(V,randcursum2w)
plot(V,prob4vrm6.^4)
plot(V,.5*prob4vrm22.^4+.5*prob4vr35.^4,'k')
legend('4:0','3:1','2:2','1:3','0:4','rand','coexpression data','50:50')
title('7.2 and 7.2w 2013','fontsize',20)
xlabel('Voltage (mV)')
ylabel('Open Probability')
hold off
axis([-60 100 0 1])


%%%%%%%%%%%%%%%%%%%%%%%%%7.4/7.5 STORY%%%%%%%%%%%%%%
probvrm46=1./(1+exp(-(V+46.2)/12.3));
probvrm28=1./(1+exp(-(V+28.6)/10.8));
probvrm34=1./(1+exp(-(V+34.0)/11.3));

prob4vrm46=1./(1+exp(-(V+75.2)/17.4));
prob4vrm28=1./(1+exp(-(V+54.5)/15.5));
prob4vrm34=1./(1+exp(-(V+60.9)/16.2));


Het745r40=prob4vrm46.^4;
Het745r31=prob4vrm46.^3.*prob4vrm28.^1;
Het745r22=prob4vrm46.^2.*prob4vrm28.^2;
Het745r13=prob4vrm46.^1.*prob4vrm28.^3;
Het745r04=prob4vrm28.^4;
randcursum745=.0625*Het745r40+.25*Het745r31+.375*Het745r22+.25*Het745r13+.0625*Het745r04;
hom5050745=.5*Het745r40+.5*Het745r04;


figure(7)
hold on
plot(V,prob4vrm46.^4)
plot(V,prob4vrm46.^3.*prob4vrm28.^1,'m')
plot(V,prob4vrm46.^2.*prob4vrm28.^2)
plot(V,prob4vrm46.^1.*prob4vrm28.^3)
plot(V,prob4vrm28.^4)
plot(V,randcursum745)
plot(V,prob4vrm34.^4)
plot(V,hom5050745,'k')
legend('4:0','3:1','2:2','1:3','0:4','rand','coexpression data','50:50')
xlabel('Voltage (mV)')
ylabel('Open Probability')
title('7.4 and 7.5','fontsize',20)
axis([-80 10 0 1])
hold off


%%%%%%%%%%7.2 and Mutant STORY MIcel 2015%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%Homomeric expression that needs to be refit and refits
%kv7.2
probvrm26=1./(1+exp(-(V+26.7)/12.6));
prob4vrm26=1./(1+exp(-(V+55.2)/16.4));
%kv7.2from 2013
probvrm22=1./(1+exp(-(V+22.9)/12.3));
prob4vrm22=1./(1+exp(-(V+50.7)/16.));
%kv7.2Q
probvrm46=1./(1+exp(-(V+46.9)/18.2));
prob4vrm46=1./(1+exp(-(V+87.8)/23.5));
%kv7.2H
probvrm56=1./(1+exp(-(V+56.1)/25.7));
prob4vrm56=1./(1+exp(-(V+111)/31.85));
%kv7.3
probvrm38=1./(1+exp(-(V+38.2)/7.2));
prob4vrm38=1./(1+exp(-(V+54.4)/9.3));

%%%%%%%%%%%%Coexpression Data to compare to
%kv72/kv72q
probvrm37=1./(1+exp(-(V+37.7)/15.3));
%kv72/kv72h
probvrm41=1./(1+exp(-(V+41.2)/15));
%kv72/kv73
probvrm33=1./(1+exp(-(V+33)/12.1));
%kv72/kv73 from 
probvrm30=1./(1+exp(-(V+30.7)/11.7));



%%%%%%%7.2 and 7.2q heteromer predictions Miceli 2016
%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%
Het7q40=prob4vrm26.^4;
Het7q31=prob4vrm26.^3.*prob4vrm46.^1;
Het7q22=prob4vrm26.^2.*prob4vrm46.^2;
Het7q13=prob4vrm26.^1.*prob4vrm46.^3;
Het7q04=prob4vrm46.^4;
randcursumq=.0625*Het7q40+.25*Het7q31+.375*Het7q22+.25*Het7q13+.0625*Het7q04;
hom50507q=.5*Het7q40+.5*Het7q04;

figure(8)
hold on
plot(V,Het7q40)
plot(V,Het7q31,'m')
plot(V,Het7q22)
plot(V,Het7q13)
plot(V,Het7q04)
plot(V,randcursumq)
plot(V,probvrm37)
plot(V,hom50507q,'k')
hold off
xlabel('Voltage (mV)')
ylabel('Open Probability')
title('7.2 and 7.2q 2016','fontsize',20)
legend('4:0','3:1','2:2','1:3','0:4','rand','coexpression data','50:50')
axis([-100 20 0 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%7.2 and 7.2H heteromer predictions Miceli 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%
Het7h40=prob4vrm26.^4;
Het7h31=prob4vrm26.^3.*prob4vrm56.^1;
Het7h22=prob4vrm26.^2.*prob4vrm56.^2;
Het7h13=prob4vrm26.^1.*prob4vrm56.^3;
Het7h04=prob4vrm56.^4;
randcursumh=.0625*Het7h40+.25*Het7h31+.375*Het7h22+.25*Het7h13+.0625*Het7h04;
hom50507h=.5*Het7h40+.5*Het7h04;

figure(9)
hold on
plot(V,Het7h40)
plot(V,Het7h31,'m')
plot(V,Het7h22)
plot(V,Het7h13)
plot(V,Het7h04)
plot(V,randcursumh)
plot(V,probvrm41)
plot(V,hom50507h,'k')
hold off
xlabel('Voltage (mV)')
ylabel('Open Probability')
title('7.2 and 7.2H 2016','fontsize',20)
legend('4:0','3:1','2:2','1:3','0:4','rand','coexpression data','50:50')
axis([-100 20 0 1])

%%%%%%%7.2 and 7.3 heteromer predictions Miceli 2016
Het7273r40=prob4vrm26.^4;
Het7273r31=prob4vrm26.^3.*prob4vrm38.^1;
Het7273r22=prob4vrm26.^2.*prob4vrm38.^2;
Het7273r13=prob4vrm26.^1.*prob4vrm38.^3;
Het7273r04=prob4vrm38.^4;
randcursum7273=.0625*Het7273r40+.25*Het7273r31+.375*Het7273r22+.25*Het7273r13+.0625*Het7273r04;
hom50507273=.5*Het7273r40+.5*Het7273r04;

figure(10)
hold on
plot(V,Het7273r40)
plot(V,Het7273r31,'m')
plot(V,Het7273r22)
plot(V,Het7273r13)
plot(V,Het7273r04)
plot(V,randcursum7273)
plot(V,probvrm33)
plot(V,hom50507273,'k')
hold off
xlabel('Voltage (mV)')
ylabel('Open Probability')
title('7.2 and 7.3 2016','fontsize',20)
legend('4:0','3:1','2:2','1:3','0:4','rand','coexpression data','50:50')
axis([-60 20 0 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%1.1/1.2 mutant 2017 hasan experiment %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%kv1WT
probvrm25=1./(1+exp(-(V+25)/12.6));
prob4vrm25=1./(1+exp(-(V+53.5)/16.4));
%%kv1f
probvr11=1./(1+exp(-(V-11)/8));
prob4vr11=1./(1+exp(-(V+7.1)/10.4));
%%%%%coexpression data
probvrm13=1./(1+exp(-(V+13)/9));



Het1f40=prob4vrm25.^4;
Het1f31=prob4vrm25.^3.*prob4vr11.^1;
Het1f22=prob4vrm25.^2.*prob4vr11.^2;
Het1f13=prob4vrm25.^1.*prob4vr11.^3;
Het1f04=prob4vr11.^4;
randcursum1f=.0625*Het1f40+.25*Het1f31+.375*Het1f22+.25*Het1f13+.0625*Het1f04;
hom50501F=.5*Het1f40+.5*Het1f04;


figure(11)
hold on
plot(V,Het1f40)
plot(V,Het1f31,'m')
plot(V,Het1f22)
plot(V,Het1f13)
plot(V,Het1f04)
plot(V,randcursum1f)
plot(V,probvrm13)
plot(V,hom50501F,'k')
hold off
xlabel('Voltage (mV)')
ylabel('Open Probability')
title('1.1 and 1.1F303','fontsize',20)
legend('4:0','3:1','2:2','1:3','0:4','rand','coexpression data','50:50')
axis([-60 50 0 1])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%1.1 mutant 2021 Miceli experiment epilepsy %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Kv1.1
probvrm26v1=1./(1+exp(-(V+26)/8.7));
prob4vrm26v1=1./(1+exp(-(V+45.7)/11.3));

%Kv11T
probvrm46vT=1./(1+exp(-(V+46.2)/8.6));
prob4vrm46vT=1./(1+exp(-(V+65.7)/11.2));

%Kv11S
probvr37vS=1./(1+exp(-(V-37.8)/8.4));
prob4vr37vS=1./(1+exp(-(V-18.8)/10.9));

%%%%%Coexpression Data%%%%%%%%
%%%Kv1/Kv1T
probvrm43T1=1./(1+exp(-(V+43.3)/9.3));

%%%Kv1/Kv1S
probvr13S1=1./(1+exp(-(V-13.1)/10.6));



figure(7)
hold on 
plot(V,probvrm26v1)
plot(V,prob4vrm26v1.^4)
plot(V,probvrm46vT)
plot(V,prob4vrm46vT.^4)
plot(V,probvr37vS)
plot(V,prob4vr37vS.^4)
hold off

%%%%%%%%%%%%%%%
%Kv1 and Kv1T mutant
%%%%%%%%%%%%%%%%
Het1T40=prob4vrm26v1.^4;
Het1T31=prob4vrm26v1.^3.*prob4vrm46vT.^1;
Het1T22=prob4vrm26v1.^2.*prob4vrm46vT.^2;
Het1T13=prob4vrm26v1.^1.*prob4vrm46vT.^3;
Het1T04=prob4vrm46vT.^4;
randcursum1T=.0625*Het1T40+.25*Het1T31+.375*Het1T22+.25*Het1T13+.0625*Het1T04;
hom50501T=.5*Het1T40+.5*Het1T04;


figure(12)
hold on
plot(V,Het1T40)
plot(V,Het1T31,'m')
plot(V,Het1T22)
plot(V,Het1T13)
plot(V,Het1T04)
plot(V,randcursum1T)
plot(V,probvrm43T1)
plot(V,hom50501T,'k')
hold off
xlabel('Voltage (mV)')
ylabel('Open Probability')
title('1.1 1.1T','fontsize',20)
legend('4:0','3:1','2:2','1:3','0:4','rand','coexpression data','50:50')
axis([-80 20 0 1])

%%%%%%%%%%%%%%%
%Kv1 and Kv1S mutant
%%%%%%%%%%%%%%%%
Het1S40=prob4vrm26v1.^4;
Het1S31=prob4vrm26v1.^3.*prob4vr37vS.^1;
Het1S22=prob4vrm26v1.^2.*prob4vr37vS.^2;
Het1S13=prob4vrm26v1.^1.*prob4vr37vS.^3;
Het1S04=prob4vr37vS.^4;
randcursum1S=.0625*Het1S40+.25*Het1S31+.375*Het1S22+.25*Het1S13+.0625*Het1S04;
hom50501S=.5*Het1S40+.5*Het1S04;

figure(13)
hold on
plot(V,Het1S40)
plot(V,Het1S31,'m')
plot(V,Het1S22)
plot(V,Het1S13)
plot(V,Het1S04)
plot(V,randcursum1S)
plot(V,probvr13S1)
plot(V,hom50501S,'k')
hold off
xlabel('Voltage (mV)')
ylabel('Open Probability')
title('1.1 1.1S','fontsize',20)
legend('4:0','3:1','2:2','1:3','0:4','rand','coexpression data','50:50')
axis([-60 60 0 1])


