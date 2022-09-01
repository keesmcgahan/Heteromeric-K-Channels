clear all
clc

%%%Must download data files and API from ChannelPedia.com (https://channelpedia.epfl.ch/)
%%%%Downlaod the 

%%%%List of Kchannel data downloaded and looked at
%Kv1.1=[8018,9391,9395,9405,9411,9416,9412,9413,9414]
%Kv1.2=[11350,11361,5401]
%kv7=[25845,10290,4930,4933]
%kv3=[2899,5516,1892]
 
%%%%%expers gives the list of channels to be run through the function

%%%%%The length of the error matrices and following for loop is based on
%%%%%length of expers.
expers=[9391,9395,9411,9412,9413,9414];
errormat=zeros(18,5);
rerrmat=zeros(6,5);
werrmat=zeros(6,5);
serrmat=zeros(6,5);

for i=1:6
   [time, data]=nwbGetProtocolTraces(expers(i), 'Activation',1);
   [timeout,dataout,fitcurve,errors,weightederrors,shorterrors]=power4(time, data);
   figure(i)
   hold on
   plot(timeout,dataout)
   plot(fitcurve{1},'g')
   plot(fitcurve{2},'m')
   plot(fitcurve{3},'r')
   plot(fitcurve{4},'k')
   plot(fitcurve{5},'y')
   hold off
xlabel('Time (s)')
ylabel('P_o(t)')
legend('Vclamp Data','p=1','p=2','p=3','p=4','p=5')
   errormat(3*(i-1)+1,:)=errors;
   errormat(3*(i-1)+2,:)=weightederrors;
   errormat(3*(i-1)+3,:)=shorterrors;
   rerrmat(i,:)=errors;
   werrmat(i,:)=weightederrors;
   serrmat(i,:)=shorterrors;
end    


function [timeout,dataout,fitcurve,errors,weightederrors,shorterrors]=power4(time,data)

adjto0=zeros(30,18);
for i=1:30
    adjto0(i,:)=data(i,:);
end 

avgadj=mean(adjto0);
data1=data-avgadj;

Voltage=zeros(1,18);
for i=1:18
    Voltage(i)=i*10-100;
end

DrivingForce=Voltage+90-1;




voltprotonly=zeros(3000,18);
time1=[0:.0001:.30];
for j=1:3001
voltprotonly(j,:)=data1(j+997,:)./DrivingForce;
 
end  

endogcur=zeros(10,18);
for i=1:10
    endogcur(i,:)=voltprotonly(i,:);
end

endogcuravg=mean(endogcur);

voltprotmend=voltprotonly-endogcuravg;
maxcur=max(max(voltprotmend(:,8:18)));

probonly=voltprotmend/maxcur;

%%%%Change the number in intrun to determine the voltage step
intrun=probonly(:,7);

%%%%Change t1 based on voltage step. THis is the time value at which the channels begin to open.
%%%%Values used for paper were 10,15,20,30,40,50 for intrun set to
%%%%12,11,10,9,8,7 respectively
t1=50;

endcur=mean(intrun(end-500:end));

fo = fitoptions('Method','NonlinearLeastSquares','Lower',[endcur^(1/4),0],'Upper',[endcur^(1/4),15],...
               'StartPoint',[endcur^(1/4),.1]);
go = fitoptions('Method','NonlinearLeastSquares','Lower',[endcur^(1/3),0],'Upper',[endcur^(1/3),15],...
               'StartPoint',[endcur^(1/3),.1]);           
ho = fitoptions('Method','NonlinearLeastSquares','Lower',[endcur^(1/2),0],'Upper',[endcur^(1/2),15],...
               'StartPoint',[endcur^(1/2),.1]);
jo = fitoptions('Method','NonlinearLeastSquares','Lower',[endcur^(1/5),0],'Upper',[endcur^(1/5),15],...
               'StartPoint',[endcur^(1/5),.1]);   
ko = fitoptions('Method','NonlinearLeastSquares','Lower',[endcur,0],'Upper',[endcur,15],...
               'StartPoint',[endcur,.1]);            
           
ft = fittype('(a*(1-exp(-x/b)))^4','options',fo);
[curve4,fof1]=fit(time1',intrun,ft);
gt = fittype('(a*(1-exp(-x/b)))^3','options',go);
[curve3,gof1]=fit(time1',intrun,gt);
ht = fittype('(a*(1-exp(-x/b)))^2','options',ho);
[curve2,hof1]=fit(time1',intrun,ht);
jt = fittype('(a*(1-exp(-x/b)))^5','options',jo);
[curve5,jof1]=fit(time1',intrun,jt);
kt = fittype('(a*(1-exp(-x/b)))','options',ko);
[curve1,kof1]=fit(time1',intrun,kt);

timeout=time1;
dataout=intrun;
fitcurve={curve1,curve2,curve3,curve4,curve5};

errors=[kof1.rmse^2,hof1.rmse^2,gof1.rmse^2,fof1.rmse^2,jof1.rmse^2];


timeshort=time1(1:t1);
timelong=time1(t1:length(time1)-1);
curve2(timeshort);
curve2(timelong);
fitshort=intrun(1:t1);
fitlong=intrun(t1:length(time1)-1);

weighterrp5=100*immse(fitshort,curve5(timeshort))+immse(fitlong,curve5(timelong));
weighterrp4=100*immse(fitshort,curve4(timeshort))+immse(fitlong,curve4(timelong));
weighterrp2=100*immse(fitshort,curve2(timeshort))+immse(fitlong,curve2(timelong));
weighterrp3=100*immse(fitshort,curve3(timeshort))+immse(fitlong,curve3(timelong));
weighterrp1=100*immse(fitshort,curve1(timeshort))+immse(fitlong,curve1(timelong));

shorterrp5=immse(fitshort,curve5(timeshort));
shorterrp4=immse(fitshort,curve4(timeshort));
shorterrp3=immse(fitshort,curve3(timeshort));
shorterrp2=immse(fitshort,curve2(timeshort));
shorterrp1=immse(fitshort,curve1(timeshort));

weightederrors=[weighterrp1,weighterrp2,weighterrp3,weighterrp4,weighterrp5];
shorterrors=[shorterrp1,shorterrp2,shorterrp3,shorterrp4,shorterrp5];
end