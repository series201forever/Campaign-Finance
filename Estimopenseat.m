
clear
% **
% make sure you get the sign of q_e right. Also note that f_d_qe, f_w_qe have the PARTY of Challenger.
% **




E_V_july8=csvread('openseatdata.csv',1,0);

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
%Others
%load ('./retire.txt');
%load ('./stateevol.mat')



dele3=find(sum(isnan(E_V_july8),2)>0); 
E_V_july8(dele3,:)=[]; %Drop NaN
E_V_july8(E_V_july8(:,2)>2002,:)=[]; %Drop year 2004 and on
deleC=find((E_V_july8(:,16)==0).*E_V_july8(:,8)==1);
E_V_july8(deleC,:)=[]; %Drop if oppornent disburse=0 and contested

%deleC2=find(E_V_july8(:,12)<5000|E_V_july8(:,12)>1200000);
%deleC2=find(E_V_july8(:,12)<5000);
%E_V_july8(deleC2,:)=[];%Drop if nonpositive saving
%deleC4=find((E_V_july8(:,16)<5000)&E_V_july8(:,8)==1);
%deleC4=find((E_V_july8(:,16)<5000|E_V_july8(:,16)>1200000)&E_V_july8(:,8)==1);
%E_V_july8(deleC4,:)=[];%Drop if nonpositive challenger spending


%Construct dataset
E_VContestFUL=find(E_V_july8(:,66)==0);           %% E_VContestFUL is the elements in which Quality of challenger A not observed.
%E_VNContestFUL=find(E_V_july8(:,8)==1);          %% E_VNContestFUL is the elements in which contest==1 :297 distinct districts

%Win_E_V=(E_V_july8(:,13)>=0.5);  % Dummy for whether incumbent won the election in (t).
%Win_E_VCT=Win_E_V;
%Win_E_VCT(E_VContestFUL,:)=[];
%IND5=find(Win_E_VCT==0);  % index for elections with incumbent loss (among those with entry).
                          % Number of distinct districts within contest=0 is 27.
                          % Number of distinct districts within contest=1 is 295.




PartyE_V=3*ones(length(E_V_july8),1)-2*E_V_july8(:,3);  %%PartyE_V=1  if candidate i is Democrat and PartyE_V=-1 if candidate i is Republican.
PartyE_VCT=PartyE_V;
PartyE_VCT(E_VContestFUL,:)=[];     %PartyE_VCT is the party indicator when there is entry next period.
% PartyE_VCTwnxt=PartyE_VCT;
% PartyE_VCTwnxt(IND5,:)=[];
% PartyE_VNCT=PartyE_V;
% PartyE_VNCT(E_VNContestFUL,:)=[];    %PartyE_VCT is the party indicator when there is NO entry next period.

PrespartyE_V=3*ones(length(E_V_july8),1)-2*E_V_july8(:,39);
SameE_V=PartyE_V==PrespartyE_V; %Same party if the two match
DifE_V=PartyE_V~=PrespartyE_V;
SameE_V=SameE_V-DifE_V;
SameE_VCT=SameE_V;
SameE_VCT(E_VContestFUL,:)=[];     %SameE_VCT is the same indicator when there is entry next period.
%SameE_VCTwnxt=SameE_VCT;
%SameE_VCTwnxt(IND5,:)=[];
%SameE_VNCT=PartyE_V;
%SameE_VNCT(E_VNContestFUL,:)=[];    %SameE_VNCT is the same indicator when there is NO entry next period.

PresdumE_V=E_V_july8(:,64); %Presdum=1 if president's incumbency is on second period
PresdumE_VCT=PresdumE_V;
PresdumE_VCT(E_VContestFUL,:)=[];     
%PresdumE_VCTwnxt=PresdumE_VCT;
%PresdumE_VCTwnxt(IND5,:)=[];
%PresdumE_VNCT=PresdumE_V;
%PresdumE_VNCT(E_VNContestFUL,:)=[]; 

MidtermE_V=E_V_july8(:,65);
 MidtermE_VCT=MidtermE_V;
 MidtermE_VCT(E_VContestFUL,:)=[];     
% MidtermE_VCTwnxt=MidtermE_VCT;
% MidtermE_VCTwnxt(IND5,:)=[];
% MidtermE_VNCT=MidtermE_V;
% MidtermE_VNCT(E_VNContestFUL,:)=[]; 

YEARE_V=E_V_july8(:,2);           % YEARE_V is the year of the election.
 YEARE_VCT=YEARE_V;
 YEARE_VCT(E_VContestFUL,:)=[];
% YEARE_VCTwnxt=YEARE_VCT;
% YEARE_VCTwnxt(IND5,:)=[];
% YEARE_VNCT=YEARE_V;
% YEARE_VNCT(E_VNContestFUL,:)=[];
% 
% XSNCEV_=[E_V_july8(:,30),E_V_july8(:,27)];       % E_V(:,30)=unemployment_NC (%), E_V(:,27)=pctwhite_NC
% XSNCEVCT_=XSNCEV_;
% XSNCEVCT_(E_VContestFUL,:)=[];
% XSNCEVCTwnxt_=XSNCEVCT_;
% XSNCEVCTwnxt_(IND5,:)=[];
% XSNCEVNCT_=XSNCEV_;
% XSNCEVNCT_(E_VNContestFUL,:)=[];
% 
% XSEV_=[E_V_july8(:,24),E_V_july8(:,21)];    % E_V(:,24)=unemployment (%), E_V(:,21)=pctwhite (%)
% XS_EVCT_=XSEV_;
% XS_EVCT_(E_VContestFUL,:)=[];
% XS_EVNCT_=XSEV_;
% XS_EVNCT_(E_VNContestFUL,:)=[];

XSEV_=[E_V_july8(:,24),E_V_july8(:,63)];    % E_V(:,24)=unemployment (%), E_V(:,63)=partisanship
 XS_EVCT_=XSEV_;
 XS_EVCT_(E_VContestFUL,:)=[];
% XS_EVCTwnxt_=XS_EVCT_;
% XS_EVCTwnxt_(IND5,:)=[];
% XS_EVNCT_=XSEV_;
% XS_EVNCT_(E_VNContestFUL,:)=[];

%TenureE_V=E_V_july8(:,15);      %Tenure=# of terms from first period.
% TenureE_VCT=TenureE_V;
% TenureE_VCT(E_VContestFUL,:)=[];
% TenureE_VCTwnxt=TenureE_VCT;
% TenureE_VCTwnxt(IND5,:)=[];
% TenureE_VNCT=TenureE_V;
% TenureE_VNCT(E_VNContestFUL,:)=[];








% LOGLOGD_E_V=log(max(ones(length(E_V_july8),1),NCE_V(:,1)));
% LOGLOGD_E_VCT=LOGLOGD_E_V;
% LOGLOGD_E_VCT(E_VContestFUL,:)=[]; 
% LOGLOGD_E_VNCT=LOGLOGD_E_V;
% LOGLOGD_E_VNCT(E_VNContestFUL,:)=[];
% LOGLOGD_E_VCTwnxt=LOGLOGD_E_VCT;
% LOGLOGD_E_VCTwnxt(IND5,:)=[];
% LOGLOGTot_E_V=log(max(ones(length(E_V_july8),1),LOGTot_NCE_V));
% LOGLOGTot_E_VCT=LOGLOGTot_E_V;
% LOGLOGTot_E_VCT(E_VContestFUL,:)=[];
% LOGLOGTot_E_VNCT=LOGLOGTot_E_V;
% LOGLOGTot_E_VNCT(E_VNContestFUL,:)=[];
% LOGLOGTot_E_VCTwnxt=LOGLOGTot_E_VCT;
% LOGLOGTot_E_VCTwnxt(IND5,:)=[];

% LOGTot_NCE_V=log(max(ones(length(E_V_july8),1),E_V_july8(:,31)));
% LOGD_NCE_V=log(max(ones(length(E_V_july8),1),E_V_july8(:,32)));
% LOGW_NCE_V=log(max(ones(length(E_V_july8),1),E_V_july8(:,33)));
% LOGWNXT_NCE_V=log(max(ones(length(E_V_july8),1),E_V_july8(:,34)));
% Tenure_NCE_V=E_V_july8(:,35); % Tenure_NC


%NCE_V=[LOGD_NCE_V,LOGW_NCE_V,LOGWNXT_NCE_V];
% NCE_VCT=NCE_V;
% NCE_VCT(E_VContestFUL,:)=[];        %NCE_V limited to samples in which there is entry.
% NCE_VCTwnxt=NCE_VCT;
% NCE_VCTwnxt(IND5,:)=[];
% NCE_VNCT=NCE_V;
% NCE_VNCT(E_VNContestFUL,:)=[];       %NCE_V limited to samples in which there is no entry.
% NCE_VTenCT=Tenure_NCE_V;
% NCE_VTenCT(E_VContestFUL,:)=[];
% NCE_VTenCTwnxt=NCE_VTenCT;
% NCE_VTenCTwnxt(IND5,:)=[];
% NCE_VTenNCT=Tenure_NCE_V;
% NCE_VTenNCT(E_VNContestFUL,:)=[];
% LOGTot_NCE_VCT=LOGTot_NCE_V;
% LOGTot_NCE_VCT(E_VContestFUL,:)=[];
% LOGTot_NCE_VNCT=LOGTot_NCE_V;
% LOGTot_NCE_VNCT(E_VNContestFUL,:)=[];
% 

% RTotDE_V=(max(0,NCE_V(:,1)).^(-1/2))./LOGTot_NCE_V;
% RTotDE_V=RTotDE_V.*(exp(LOGTot_NCE_V)./exp(NCE_V(:,1)));
% RTotDE_VCT=(max(0,NCE_VCT(:,1)).^(-1/2))./LOGTot_NCE_VCT;
% RTotDE_VCT=RTotDE_VCT.*(exp(LOGTot_NCE_VCT)./exp(NCE_VCT(:,1)));
% RTotDE_VCTwnxt=RTotDE_VCT;
% RTotDE_VCTwnxt(IND5,:)=[];
% RTotDE_VNCT=(max(0,NCE_VNCT(:,1)).^(-1/2))./LOGTot_NCE_VNCT;
% RTotDE_VNCT=RTotDE_VNCT.*(exp(LOGTot_NCE_VNCT)./exp(NCE_VNCT(:,1)));


LOGD_E_V=log(max(ones(length(E_V_july8),1),E_V_july8(:,10)));  %log(Spend)
 LOGD_E_VCT=LOGD_E_V;
 LOGD_E_VCT(E_VContestFUL,:)=[];                             %log(Spend) limited to contested.
% LOGD_E_VNCT=LOGD_E_V;
% LOGD_E_VNCT(E_VNContestFUL,:)=[];                            %log(Spend) limited to uncontested.
% 


LOGW_E_V=log(max(ones(length(E_V_july8),1),E_V_july8(:,11)));  %log(begining cash)
 LOGW_E_VCT=LOGW_E_V;
 LOGW_E_VCT(E_VContestFUL,:)=[];
% LOGW_E_VCTwnxt=LOGW_E_VCT;
% LOGW_E_VCTwnxt(IND5,:)=[];
% LOGW_E_VNCT=LOGW_E_V;
% LOGW_E_VNCT(E_VNContestFUL,:)=[];


LOGW_NXT_E_V=log(max(ones(length(E_V_july8),1),E_V_july8(:,12)));  %log(realendcash)
 LOGW_NXT_E_VCT=LOGW_NXT_E_V;
 LOGW_NXT_E_VCT(E_VContestFUL,:)=[];
% LOGW_NXT_E_VCTwnxt=LOGW_NXT_E_VCT;
% LOGW_NXT_E_VCTwnxt(IND5,:)=[];
% LOGW_NXT_E_VNCT=LOGW_NXT_E_V;
% LOGW_NXT_E_VNCT(E_VNContestFUL,:)=[];
% IND6CT=find(LOGW_NXT_E_VCT==0);                      % IND6CT is the index where realendcash_i ==0
% IND6NCT=find(LOGW_NXT_E_VNCT==0);                    % IND6NCT is the index where realendcash_i ==0

LOGW_NXT_E_C=log(max(ones(length(E_V_july8),1),E_V_july8(:,18)));  %log(realendcash_c)
LOGW_NXT_E_C_VCT=LOGW_NXT_E_C;
LOGW_NXT_E_C_VCT(E_VContestFUL,:)=[];
% IND7=find(LOGW_NXT_E_VC==0|LOGW_NXT_E_VCT==0);   % IND7 is the index where realedcash_c==0 or realendcash_I==0 (among E_VContestFUL)
LOGD_E_C_VC=log(max(ones(length(E_V_july8),1),E_V_july8(:,16)));  %log(realdisburse_c)
LOGD_E_C_VCT=LOGD_E_C_VC;
LOGD_E_C_VCT(E_VContestFUL,:)=[];

LOGTotal_E_V=log(max(ones(length(E_V_july8),1),E_V_july8(:,9)));  %log(realtotal)
 LOGTotal_E_VCT=LOGTotal_E_V;
 LOGTotal_E_VCT(E_VContestFUL,:)=[];
% LOGTotal_E_VNCT=LOGTotal_E_V;
% LOGTotal_E_VNCT(E_VNContestFUL,:)=[];
LOGTot_E_C=log(max(ones(length(E_V_july8),1),E_V_july8(:,4)));  % log(realtotal_c)
LOGTot_E_C_VCT=LOGTot_E_C;
LOGTot_E_C_VCT(E_VContestFUL,:)=[];

VSEV=E_V_july8(:,13);
 VSEVCT=VSEV;
 VSEVCT(E_VContestFUL,:)=[];
% 
% 




%%B-Spline for q_I(RTotDE_V), where RTotDE_V=(NCE_VCT(:,1)./LOGTot_NCE_V)  etc.%%
% take knots to be between 0.88 to 1.025 with 8 knots (8 basis functions). (Almost all lie

%  fineness=9;
% % mesh=quantile(RTotD_NC,linspace(0,1,fineness).');
% mesh=quantile(RTotDE_V,linspace(0,1,fineness).');
% X_KnotEV1=(RTotDE_V<mesh(2,1)).*(1-(RTotDE_V-mesh(1,1))/(mesh(2,1)-mesh(1,1)));
% for i=0:(numel(mesh)-3)
%     PLUS=(RTotDE_V>=mesh(i+1,1)).*(RTotDE_V<mesh(i+2,1)).*((RTotDE_V-mesh(i+1,1))/(mesh(i+2,1)-mesh(i+1,1)))...
%         +(RTotDE_V>=mesh(i+2,1)).*(RTotDE_V<mesh(i+3,1)).*(1-(RTotDE_V-mesh(i+2,1))/(mesh(i+3,1)-mesh(i+2,1)));
%     X_KnotEV1=[X_KnotEV1,PLUS];
% end
% PLUS=(RTotDE_V>=mesh((numel(mesh)-1),1)).*(RTotDE_V-mesh((numel(mesh)-1),1))/(mesh(numel(mesh),1)-mesh((numel(mesh)-1),1));
% X_KnotEV1=[X_KnotEV1,PLUS];

% %First column:redundant
% X_KnotEV1=X_KnotEV1(:,2:fineness);
% X_KnotE_VCT=X_KnotEV1;
% X_KnotE_VCT(E_VContestFUL,:)=[];
% X_KnotE_VNCT=X_KnotEV1;
% X_KnotE_VNCT(E_VNContestFUL,:)=[];


XQEV=E_V_july8(:,66); % q_I
 XQEVCT=XQEV;
 XQEVCT(E_VContestFUL,:)=[];
% XQEVNCT2=XQEV;
% XQEVNCT2(E_VNContestFUL,:)=[];
% XQEVCTwnxt2=XQEVCT2;
% XQEVCTwnxt2(IND5,:)=[];

%%%%%%%%%%%%%%%%%%%%%%%% 
% Defining variables: done
%%%%%%%%%%%%%%%%%%%%%%%%

%%

%Continuation payoff: inputed from incumbents' problem
Continuei=max(0,[ones(length(LOGW_NXT_E_V),1),LOGW_NXT_E_V,LOGW_NXT_E_V.^2/10,LOGW_NXT_E_V.^3/100,LOGW_NXT_E_V.^4/1000,...
    XQEV,LOGW_NXT_E_V.*XQEV,LOGW_NXT_E_V.^2/10.*XQEV,LOGW_NXT_E_V.^3/100.*XQEV,LOGW_NXT_E_V.^4/1000.*XQEV,...
    zeros(length(LOGW_NXT_E_V),1),XSEV_(:,1).*SameE_V,XSEV_(:,2).*PartyE_V,...
    PresdumE_V,MidtermE_V]*coef3rdstep);
Continuei_VCT=Continuei;
Continuei_VCT(E_VContestFUL,:)=[];


%Calculate derivative of value function
 Deriv=[zeros(length(Continuei),1),ones(length(Continuei),1),2*LOGW_NXT_E_V/10,3*LOGW_NXT_E_V.^2/100,4*LOGW_NXT_E_V.^3/1000,...
     zeros(length(Continuei),1),XQEV,2*LOGW_NXT_E_V/10.*XQEV,3*LOGW_NXT_E_V.^2/100.*XQEV,4*LOGW_NXT_E_V.^3/1000.*XQEV,zeros(length(Continuei),5)]*coef3rdstep;
 Derivhist1=Deriv;
 Deriv=max(Deriv,0.00001);
Deriv_VCT=Deriv;
Deriv_VCT(E_VContestFUL,:)=[];

%%

%Specification check
figure(1)
scatter(XQEVCT(LOGW_NXT_E_VCT>0&Continuei_VCT>0),Continuei_VCT(LOGW_NXT_E_VCT>0&Continuei_VCT>0))
figure(2)
scatter(XQEVCT(LOGW_NXT_E_VCT>0&Continuei_VCT<=0),Continuei_VCT(LOGW_NXT_E_VCT>0&Continuei_VCT<=0))
figure(3)
scatter(XQEVCT(LOGW_NXT_E_VCT==0),Continuei_VCT(LOGW_NXT_E_VCT==0))
%scatter(XQEV(LOGW_NXT_E_VC>0&Continuei>0),Deriv(LOGW_NXT_E_VC>0&Continuei>0))
%scatter(XQEV(Deriv>0),Deriv(Deriv>0))
length(find(LOGW_NXT_E_VCT>0))
length(find(Continuei_VCT>0))
%length(find(Deriv>0))
%% 
 D=@(test) coef3rdstep(6)+test*coef3rdstep(7)+test.^2/10*coef3rdstep(8)+test.^3/100*coef3rdstep(9)+test.^4/1000*coef3rdstep(10);
 space1=9:0.5:15;
 mat1=D(space1);
 scatter(space1,mat1);
%%
B=@(test,test2)coef3rdstep(2)+2*test2/10*coef3rdstep(3)+3*test2.^2/100*coef3rdstep(4)+4*test2.^3/1000*coef3rdstep(5)+test*coef3rdstep(7)+2*test2/10.*test*coef3rdstep(8)...
   +3*test2.^2/100.*test*coef3rdstep(9)+4*test2.^3/1000.*test*coef3rdstep(10);
space1=0.01:0.01:0.2;
space2=5:0.5:15;
mat1=B(space1,11);
mat2=B(0.07,space2);
scatter(space2,mat2)

%% 

%Estimation of cost parameters
datasetVCT=[XQEVCT,LOGW_NXT_E_VCT,LOGTotal_E_VCT,LOGD_E_VCT,Continuei_VCT,Deriv_VCT];
%datasetVCT((LOGW_NXT_E_VCT==0)|(XQEVCT>0.15),:)=[];


%%

%Estimate beta_I
para0=0.04;
options=optimset('MaxIter',200000,'MaxFunEvals',1000000,'Display','iter');
[mintheta,SRR]=fminsearch(@(theta4) Minimizeopenseat(mintheta2ndstage,est3rdstage,theta4,q_C_E_VCT,datasetVCT),para0,options);

estopen=mintheta;
save estopen estopen

%%
%Compute q_c
thetain=mintheta2ndstage;
theta4=estopen;
theta3=est3rdstage;



vdelta=0.90;
alpha=1/2;
beta=2;
ben1=abs(thetain(2));
ben2=ben1;
sig=abs(thetain(3));


costc=abs(theta3(1));
a=theta3(2);
b=theta3(3);
c=theta3(4);
B_I=abs(theta4(1));


OUT=((beta*costc*(alpha/beta)*ben2*(ones(length(XQEVCT),1)+a*(XQEVCT+abs(min(q_C_E_VCT))+0.001)+b*(XQEVCT+abs(min(q_C_E_VCT))+0.001).^2+c*(XQEVCT+abs(min(q_C_E_VCT))+0.001).^3).*LOGTotal_E_VCT.^(beta-1)).*(1./exp(LOGTotal_E_VCT(:,1))))./((vdelta*Deriv_VCT).*(1./exp(LOGW_NXT_E_VCT)));
ako=min(max(OUT,0.000001),0.999999);
BX1=norminv(ako);



regressor=[LOGD_E_VCT,LOGD_E_C_VCT,XS_EVCT_(:,1).*SameE_VCT,XS_EVCT_(:,2).*PartyE_VCT];
coef=[B_I;(-1*B_I);Est2([3,4])];
qi=XQEVCT;
%Sample selection
%   regressor((OUT>1|OUT<0),:)=[];
%   BX1((OUT>1|OUT<0),:)=[];
%   qi((OUT>1|OUT<0),:)=[];%(LOGW_NXT_E_VCT==0)|(XQEVCT>0.15)

qc=regressor*coef+qi-sig*BX1;


save qcopen qc





