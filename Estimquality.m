
clear
% **
% make sure you get the sign of q_e right. Also note that f_d_qe, f_w_qe have the PARTY of Challenger.
% **


E_V_july8=csvread('qualitydata.csv',1,0);

%Estimates from the first stage
% load('./Est1.mat');
 load('./Est2.mat');
% load('./Est3.mat');
% load('./Est411.mat');
% load('./Est412.mat');
% load('./Est413.mat');
% load('./Est421.mat');
% load('./Est422.mat');
% load('./Est423.mat');
% load('./Est4311.mat');
% load('./Est4312.mat');
% load('./Est4313.mat');
% load('./Est4321.mat');
% load('./Est4322.mat');
% load('./Est4323.mat');
% load('./Est51.mat');
% load('./Est52.mat');
% load('./Est53.mat');

%Estimates from second stage
load('./q_C_E_VCT.mat');
%load ('./presseq.mat')
%load ('./DEL_Pen.mat')
%load ('./Continue.mat')
load ('./coef3rdstep.mat')
load ('./est2ndstage.mat')
load ('./est3rdstage.mat')
load ('./estopen.mat')
%Others
%load ('./retire.txt');
%load ('./stateevol.mat')



dele3=find(sum(isnan(E_V_july8),2)>0); 
E_V_july8(dele3,:)=[]; %Drop NaN
E_V_july8(E_V_july8(:,2)>2002,:)=[]; %Drop year 2004 and on

%deleC2=find(E_V_july8(:,12)<5000|E_V_july8(:,12)>1200000);
deleC2=find(E_V_july8(:,9)<5000);
E_V_july8(deleC2,:)=[];%Drop if nonpositive saving
deleC4=find((E_V_july8(:,15)<5000)&(E_V_july8(:,12)==1|E_V_july8(:,4)==1));
%deleC4=find((E_V_july8(:,16)<5000|E_V_july8(:,16)>1200000)&E_V_july8(:,8)==1);
E_V_july8(deleC4,:)=[];%Drop if nonpositive challenger spending


%Construct dataset
%E_VContestFUL=find(E_V_july8(:,66)==0);           %% E_VContestFUL is the elements in which Quality of challenger A not observed.
%E_VNContestFUL=find(E_V_july8(:,8)==1);          %% E_VNContestFUL is the elements in which contest==1 :297 distinct districts

%Win_E_V=(E_V_july8(:,13)>=0.5);  % Dummy for whether incumbent won the election in (t).
%Win_E_VCT=Win_E_V;
%Win_E_VCT(E_VContestFUL,:)=[];
%IND5=find(Win_E_VCT==0);  % index for elections with incumbent loss (among those with entry).
                          % Number of distinct districts within contest=0 is 27.
                          % Number of distinct districts within contest=1 is 295.




PartyE_V=3*ones(length(E_V_july8),1)-2*E_V_july8(:,1);  %%PartyE_V=1  if candidate i is Democrat and PartyE_V=-1 if candidate i is Republican.
PrespartyE_V=3*ones(length(E_V_july8),1)-2*E_V_july8(:,3);
SameE_V=PartyE_V==PrespartyE_V; %Same party if the two match
DifE_V=PartyE_V~=PrespartyE_V;
SameE_V=SameE_V-DifE_V;
PresdumE_V=E_V_july8(:,22); %Presdum=1 if president's incumbency is on second period
MidtermE_V=E_V_july8(:,23);
YEARE_V=E_V_july8(:,2);           % YEARE_V is the year of the election.
XSEV_=[E_V_july8(:,5),E_V_july8(:,10)];    % E_V(:,5)=unemployment (%), E_V(:,10)=partisanship
TenureE_V=E_V_july8(:,11);  %log(realtotal)
LOGD_E_V=log(max(ones(length(E_V_july8),1),E_V_july8(:,8)));  %log(Spend)
LOGW_E_V=log(max(ones(length(E_V_july8),1),E_V_july8(:,6)));  %log(begining cash)
LOGW_NXT_E_V=log(max(ones(length(E_V_july8),1),E_V_july8(:,9)));  %log(realendcash)
LOGTotal_E_V=log(max(ones(length(E_V_july8),1),E_V_july8(:,7)));  %log(realtotal)
LOGD_E_VC=log(max(ones(length(E_V_july8),1),E_V_july8(:,15)));  %log(Spend)
LOGW_NXT_E_VC=log(max(ones(length(E_V_july8),1),E_V_july8(:,17)));  %log(realendcash)
LOGTotal_E_VC=log(max(ones(length(E_V_july8),1),E_V_july8(:,14)));  %log(realtotal)
qi=E_V_july8(:,20);
qc=E_V_july8(:,21);
idi=E_V_july8(:,18);
idc=E_V_july8(:,19);
contested=E_V_july8(:,12);
open=E_V_july8(:,4);
%%%%%%%%%%%%%%%%%%%%%%%% 
% Defining variables: done
%%%%%%%%%%%%%%%%%%%%%%%%

%%

%Construct datset
datasetV=[LOGW_NXT_E_V,LOGTotal_E_V,LOGD_E_V,PartyE_V,SameE_V,PresdumE_V,MidtermE_V,XSEV_,TenureE_V,LOGW_NXT_E_VC,LOGTotal_E_VC,LOGD_E_VC];
q=[qi,qc];
id=[idi,idc];
%dataset for contested election
datasetc=datasetV(contested==1,:);
qcon=q(contested==1,:);
idcon=id(contested==1,:);
%dataset for openseat
dataseto=datasetV(open==1,:);
qo=q(open==1,:);
ido=id(open==1,:);

%%
%Estimate quality for contested
qidistc=zeros(length(datasetc),1);
qcdistc=qidistc;

q0=[0.05,-0.05];
tic
for i=1:length(datasetc);
options=optimset('MaxIter',200000,'MaxFunEvals',1000000,'Display','none');
[minq,SRR,flag]=fminsearch(@(XQEV) Minimizequality(Est2,mintheta2ndstage,est3rdstage,coef3rdstep,estopen,q_C_E_VCT,datasetc(i,:),XQEV,1),q0,options);
%[minq,SRR,flag]=fminunc(@(XQEV) Minimizequality(Est2,mintheta2ndstage,est3rdstage,coef3rdstep,estopen,q_C_E_VCT,datasetc(i,:),XQEV,1),q0,options);
if flag==1
qidistc(i)=minq(1);
qcdistc(i)=minq(2);
end
end
toc
save qcondist2 qidistc qcdistc 

%%
%Estimate quality for open
qidisto=zeros(length(dataseto),1);
qcdisto=qidisto;
q0=[0,0];
tic
for i=1:length(dataseto);
options=optimset('MaxIter',200000,'MaxFunEvals',1000000,'Display','none');
[minq,SRR,flag]=fminsearch(@(XQEV) Minimizequality(Est2,mintheta2ndstage,est3rdstage,coef3rdstep,estopen,q_C_E_VCT,dataseto(i,:),XQEV,3),q0,options);
%[minq,SRR,flag]=fminunc(@(XQEV) Minimizequality(Est2,mintheta2ndstage,est3rdstage,coef3rdstep,estopen,q_C_E_VCT,datasetc(i,:),XQEV,1),q0,options);
if flag==1
qidisto(i)=minq(1);
qcdisto(i)=minq(2);
end
end
toc
save qopendist2 qidisto qcdisto 
%%
%OUT specification check
datasetV=datasetc;
LOGW_NXT_E_V=datasetV(1,1);
LOGTotal_E_V=datasetV(1,2);
LOGD_E_V=datasetV(1,3);
PartyE_V=datasetV(1,4);
SameE_V=datasetV(1,5);
PresdumE_V=datasetV(1,6);
MidtermE_V=datasetV(1,7);
XSEV_=datasetV(1,8:9);
TenureE_V=datasetV(1,10);
LOGW_NXT_E_VC=datasetV(1,11);
LOGTotal_E_VC=datasetV(1,12);
LOGD_E_VC=datasetV(1,13);

vdelta=0.90;
% thetaS=thetain(1:2,1);
% B_I=thetain(1,1);
% B_C=thetain(2,1);
% B_T=thetain(5,1);
% thetaS2=thetain(3:4,1);


alpha=1/2;
beta=2;
ben1=abs(mintheta2ndstage(2));
ben2=ben1;
sig=abs(mintheta2ndstage(3));

costi=abs(mintheta2ndstage(1));
a=mintheta2ndstage(4);
b=mintheta2ndstage(5);
c=mintheta2ndstage(6);

costc=abs(est3rdstage(1));
ac=est3rdstage(2);
bc=est3rdstage(3);
cc=est3rdstage(4);


B_I=Est2(1);
B_C=Est2(2);
B_O=abs(estopen(1));
B_state=Est2(3:4);
B_T=Est2(5);
Derivi=@(qi)[zeros(length(LOGW_NXT_E_V),1),ones(length(LOGW_NXT_E_V),1),2*LOGW_NXT_E_V/10,3*LOGW_NXT_E_V.^2/100,4*LOGW_NXT_E_V.^3/1000,...
     zeros(length(LOGW_NXT_E_V),1),qi,2*LOGW_NXT_E_V/10.*qi,3*LOGW_NXT_E_V.^2/100.*qi,4*LOGW_NXT_E_V.^3/1000.*qi,zeros(length(LOGW_NXT_E_V),5)]*coef3rdstep;
Derivc=@(qc)[zeros(length(LOGW_NXT_E_VC),1),ones(length(LOGW_NXT_E_VC),1),2*LOGW_NXT_E_VC/10,3*LOGW_NXT_E_VC.^2/100,4*LOGW_NXT_E_VC.^3/1000,...
     zeros(length(LOGW_NXT_E_VC),1),qc,2*LOGW_NXT_E_VC/10.*qc,3*LOGW_NXT_E_VC.^2/100.*qc,4*LOGW_NXT_E_VC.^3/1000.*qc,zeros(length(LOGW_NXT_E_VC),5)]*coef3rdstep;

OUTI=@(qi)((beta*costi*(alpha/beta)*ben2*(ones(length(qi),1)+a*qi+b*qi.^2+c*qi.^3).*LOGTotal_E_V.^(beta-1)).*(1./exp(LOGTotal_E_V(:,1))))./((vdelta*Derivi(qi)).*(1./exp(LOGW_NXT_E_V)));
OUTC=@(qc)((beta*costc*(alpha/beta)*ben2*(ones(length(qc),1)+ac*(qc+abs(min(q_C_E_VCT))+0.001)+bc*(qc+abs(min(q_C_E_VCT))+0.001).^2+cc*(qc+abs(min(q_C_E_VCT))+0.001).^3).*LOGTotal_E_VC.^(beta-1)).*(1./exp(LOGTotal_E_VC(:,1))))./((vdelta*Derivc(qc)).*(1./exp(LOGW_NXT_E_VC)));


space=linspace(-0.2+abs(min(q_C_E_VCT)),0.2+abs(min(q_C_E_VCT)),10000);
outmatc=zeros(length(space),1);
for i= 1:length(space)
    outmatc(i)=OUTC(space(i));
end

space2=linspace(0.05,0.065,10000);
outmati=zeros(length(space2),1);
for i= 1:length(space2)
    outmati(i)=OUTI(space2(i));
end
figure(1)
plot(space-abs(min(q_C_E_VCT)),outmatc)
figure(2)
plot(space2,outmati)
%%

%Quality data restoration
load qcondist2
load qopendist2

%Indicator of sample truncation
%for contested
outc=zeros(length(datasetc),1);
OUTIc=outc;
OUTCc=outc;
BX1ic=outc;
BX1cc=outc;
SRR11ic=outc;
SRR11cc=outc;
testc=outc;
for i=1:length(datasetc)
[~,outc(i),OUTIc(i),OUTCc(i),BX1ic(i),BX1cc(i),SRR11ic(i),SRR11cc(i),testc(i)]=Minimizequality(Est2,mintheta2ndstage,est3rdstage,coef3rdstep,estopen,q_C_E_VCT,datasetc(i,:),[qidistc(i),qcdistc(i)],1);
end
sum(outc)

%for open
outo=zeros(length(dataseto),1);
OUTIo=outo;
OUTCo=outo;
BX1io=outo;
BX1co=outo;
SRR11io=outo;
SRR11co=outo;
testo=outo;
for i=1:length(dataseto)
[~,outo(i),OUTIo(i),OUTCo(i),BX1io(i),BX1co(i),SRR11io(i),SRR11co(i),testo(i)]=Minimizequality(Est2,mintheta2ndstage,est3rdstage,coef3rdstep,estopen,q_C_E_VCT,dataseto(i,:),[qidisto(i),qcdisto(i)],3);
end
sum(outo)

%Drop truncated samples
qdistc=[qidistc,qcdistc];
idcon2=idcon;
idcon2(min(qdistc,[],2)<-0.5|outc==1,:)=[];
qcon2=qcon;
qcon2(min(qdistc,[],2)<-0.5|outc==1,:)=[];
qdistc(min(qdistc,[],2)<-0.5|outc==1,:)=[];

qdisto=[qidisto,qcdisto];
ido2=ido;
ido2(min(qdisto,[],2)<-0.5|outo==1,:)=[];
qo2=qo;
qo2(min(qdisto,[],2)<-0.5|outo==1,:)=[];
qdisto(min(qdisto,[],2)<-0.5|outo==1,:)=[];

%Define incumbent quality and challenger quality
qualinc=[idcon2(:,1),qdistc(:,1),qcon2(:,1),ones(length(idcon2),1),zeros(length(idcon2),1)];
qualchal=[idcon2(:,2),qdistc(:,2),qcon2(:,2),ones(length(idcon2),1),zeros(length(idcon2),1);...
    ido2(:,1),qdisto(:,1),qo2(:,1),zeros(length(ido2),1),ones(length(ido2),1);...
    ido2(:,2),qdisto(:,2),qo2(:,2),zeros(length(ido2),1),ones(length(ido2),1)];

%%
figure(1)
scatter(OUTIc(abs(OUTIc)<2),qidistc(abs(OUTIc)<2))%(abs(qidistc)<0.065)
figure(2)
scatter(OUTCc(abs(OUTCc)<2),qcdistc(abs(OUTCc)<2))%(abs(qcdistc)<0.065)
figure(3)
scatter(OUTCo(abs(OUTCo)<2),qcdisto(abs(OUTCo)<2))%(abs(qidistc)<0.065)

%%
%Result check
%Define variables
qdistinc=qualinc(:,2);
qdistchal=qualchal(:,2);
qidif=qdistinc-qualinc(:,3);
qcdif=qdistchal-qualchal(:,3);
qidif(qualinc(:,3)==0,:)=[];
qcdif(qualchal(:,3)==0,:)=[];

contesti=qualinc(:,4);
contestc=qualchal(:,4);

%%
%Raw distribution
bins=linspace(-0.5,0.5,100);
figure(1)
histogram(qdistinc(contesti==1,:),bins)
hold on
%figure(2)
histogram(qdistchal,bins)
%figure(3)
%hist(qdistchal(contestc==0,:))

%Compare them with original estimation
load XQEV2
figure(2)
histogram(XQEV2,bins)
hold on
histogram(q_C_E_VCT,bins)


%%
Derivfunc=@(LOGW_NXT_E_VC,qc)[zeros(length(LOGW_NXT_E_VC),1),ones(length(LOGW_NXT_E_VC),1),2*LOGW_NXT_E_VC/10,3*LOGW_NXT_E_VC.^2/100,4*LOGW_NXT_E_VC.^3/1000,...
     zeros(length(LOGW_NXT_E_VC),1),qc,2*LOGW_NXT_E_VC/10.*qc,3*LOGW_NXT_E_VC.^2/100.*qc,4*LOGW_NXT_E_VC.^3/1000.*qc,zeros(length(LOGW_NXT_E_VC),5)]*coef3rdstep;
 space=linspace(0.05,0.065,1000).';
 l=length(space);
 meanw=12*ones(l,1);
 outmat=Derivfunc(meanw,space);
 
 scatter(space,outmat);
%%
%Difference from original estimates
figure(1)
histogram(qidif)

figure(2)
histogram(qcdif)

%%
%Across time variation
aux=[qualinc;qualchal];
qdist=zeros(max(max(idi),max(idc)),1);
for i=1:max(max(idi),max(idc));
subq=aux(aux(:,1)==i,:);
if length(subq)>1
qdist(i)=sqrt(var(subq(:,2)));
if qdist(i)>0.3
    i
end
end
end
qdist=qdist(qdist>0);
histogram(qdist)

%%
%save data
qual=zeros(length(E_V_july8),3);

qual(contested==1,:)=[qidistc,qcdistc,outc];
qual(open==1,:)=[qidisto,qcdisto,outo];


data2=[E_V_july8,qual];
dlmwrite('contestedresult.csv',data2,'delimiter', ',', 'precision', 16)





% 
% 
% 
% %How many elections hit the boundary?
% load qopendist
% out=zeros(length(datasetc),1);
% OUTI=out;
% OUTC=out;
% for i=1:length(datasetc)
% [~,out(i),OUTI(i),OUTC(i)]=Minimizequality(Est2,mintheta2ndstage,est3rdstage,coef3rdstep,estopen,q_C_E_VCT,datasetc(i,:),[qidistc(i),qcdistc(i)],1);
% end
% sum(out)
% %%
% load qcondist
% data2=[datasetc,qcon,idcon,qidistc,qcdistc,out];
% dlmwrite('contestedresult.csv',data2,'delimiter', ',', 'precision', 16)
% %%
% load qcondist
% %Drop if estimated quality <-0.5
% qdist=[qidistc,qcdistc];
% idcon2=idcon;
% idcon2(min(qdist,[],2)<-0.5|out==1,:)=[];
% qcon2=qcon;
% qcon2(min(qdist,[],2)<-0.5|out==1,:)=[];
% qdist(min(qdist,[],2)<-0.5|out==1,:)=[];
% qidistc=qdist(:,1);
% qcdistc=qdist(:,2);
% 
% 
% 
% 
% 
% %Raw distribution
% figure(1)
% hist(qidistc(qidistc>-0.5))
% 
% figure(2)
% hist(qcdistc(qcdistc>-0.5))
% %%
% %Difference from original estimates
% qcheck=qcon2-[qidistc,qcdistc];
% qcheck(qcon2(:,1)==0,:)=[];
% figure(1)
% hist(qcheck(:,1))
% 
% figure(2)
% hist(qcheck(:,2))
% %%
% %Across time variation
% aux=[idcon2(:,1),qidistc;idcon2(:,2),qcdistc];
% qdist=zeros(max(max(idi),max(idc)),1);
% for i=1:max(max(idi),max(idc));
% subq=aux(aux(:,1)==i,:);
% if length(subq)>1
% qdist(i)=sqrt(var(subq(:,2)));
% if qdist(i)>0.3
%     i
% end
% end
% end
% %%
% qdist=qdist(qdist>0);
% hist(qdist)