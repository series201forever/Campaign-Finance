
clear
% **
% make sure you get the sign of q_e right. Also note that f_d_qe, f_w_qe have the PARTY of Challenger.
% **




E_V_july8=csvread('E_Vjuly8partisan.csv',1,0);

%Estimates from the first stage
load('./Est1.mat');
load('./Est2.mat');
load('./Est3.mat');
load('./Est411.mat');
load('./Est412.mat');
load('./Est413.mat');
load('./Est421.mat');
load('./Est422.mat');
load('./Est423.mat');
load('./Est4311.mat');
load('./Est4312.mat');
load('./Est4313.mat');
load('./Est4321.mat');
load('./Est4322.mat');
load('./Est4323.mat');
load('./Est51.mat');
load('./Est52.mat');
load('./Est53.mat');

%Estimates from second stage
load('./q_C_E_VCT.mat');
%load ('./presseq.mat')
%load ('./DEL_Pen.mat')
load ('./Continue.mat')
%load ('./coef2ndstep.mat')
load ('./est2ndstage.mat')

%Others
%load ('./retire.txt');
%load ('./stateevol.mat')



dele3=find(sum(isnan(E_V_july8),2)>0); 
E_V_july8(dele3,:)=[]; %Drop NaN
E_V_july8(E_V_july8(:,2)>2002,:)=[]; %Drop year 2004 and on
deleC=find((E_V_july8(:,16)==0).*E_V_july8(:,8)==1);
E_V_july8(deleC,:)=[]; %Drop if oppornent disburse=0 and contested
E_V_july8(174:177,:)=[]; %Drop an Rtotd outlier
deleC2=find(E_V_july8(:,12)<5000|E_V_july8(:,12)>1200000);
%deleC2=find(E_V_july8(:,12)<5000);
E_V_july8(deleC2,:)=[];%Drop if nonpositive saving
deleC3=find(E_V_july8(:,11)<5000);
E_V_july8(deleC3,:)=[];%Drop if nonpositive begcash
%deleC4=find((E_V_july8(:,16)<5000)&E_V_july8(:,8)==1);
deleC4=find((E_V_july8(:,16)<5000|E_V_july8(:,16)>1200000)&E_V_july8(:,8)==1);
E_V_july8(deleC4,:)=[];%Drop if nonpositive challenger spending


%Construct dataset
E_VContestFUL=find(E_V_july8(:,8)==0);           %% E_VContestFUL is the elements in which contest==0 :277 distinct districts
E_VNContestFUL=find(E_V_july8(:,8)==1);          %% E_VNContestFUL is the elements in which contest==1 :297 distinct districts

Win_E_V=(E_V_july8(:,13)>=0.5);  % Dummy for whether incumbent won the election in (t).
Win_E_VCT=Win_E_V;
Win_E_VCT(E_VContestFUL,:)=[];
IND5=find(Win_E_VCT==0);  % index for elections with incumbent loss (among those with entry).
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
SameE_VCTwnxt=SameE_VCT;
SameE_VCTwnxt(IND5,:)=[];
SameE_VNCT=PartyE_V;
SameE_VNCT(E_VNContestFUL,:)=[];    %SameE_VNCT is the same indicator when there is NO entry next period.

PresdumE_V=E_V_july8(:,64); %Presdum=1 if president's incumbency is on second period
PresdumE_VCT=PresdumE_V;
PresdumE_VCT(E_VContestFUL,:)=[];     
PresdumE_VCTwnxt=PresdumE_VCT;
PresdumE_VCTwnxt(IND5,:)=[];
PresdumE_VNCT=PresdumE_V;
PresdumE_VNCT(E_VNContestFUL,:)=[]; 

MidtermE_V=E_V_july8(:,65);
MidtermE_VCT=MidtermE_V;
MidtermE_VCT(E_VContestFUL,:)=[];     
MidtermE_VCTwnxt=MidtermE_VCT;
MidtermE_VCTwnxt(IND5,:)=[];
MidtermE_VNCT=MidtermE_V;
MidtermE_VNCT(E_VNContestFUL,:)=[]; 

YEARE_V=E_V_july8(:,2);           % YEARE_V is the year of the election.
% YEARE_VCT=YEARE_V;
% YEARE_VCT(E_VContestFUL,:)=[];
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
XS_EVCTwnxt_=XS_EVCT_;
XS_EVCTwnxt_(IND5,:)=[];
XS_EVNCT_=XSEV_;
XS_EVNCT_(E_VNContestFUL,:)=[];

TenureE_V=E_V_july8(:,15);      %Tenure=# of terms from first period.
TenureE_VCT=TenureE_V;
TenureE_VCT(E_VContestFUL,:)=[];
TenureE_VCTwnxt=TenureE_VCT;
TenureE_VCTwnxt(IND5,:)=[];
TenureE_VNCT=TenureE_V;
TenureE_VNCT(E_VNContestFUL,:)=[];








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

LOGTot_NCE_V=log(max(ones(length(E_V_july8),1),E_V_july8(:,31)));
LOGD_NCE_V=log(max(ones(length(E_V_july8),1),E_V_july8(:,32)));
LOGW_NCE_V=log(max(ones(length(E_V_july8),1),E_V_july8(:,33)));
LOGWNXT_NCE_V=log(max(ones(length(E_V_july8),1),E_V_july8(:,34)));
% Tenure_NCE_V=E_V_july8(:,35); % Tenure_NC


NCE_V=[LOGD_NCE_V,LOGW_NCE_V,LOGWNXT_NCE_V];
NCE_VCT=NCE_V;
NCE_VCT(E_VContestFUL,:)=[];        %NCE_V limited to samples in which there is entry.
% NCE_VCTwnxt=NCE_VCT;
% NCE_VCTwnxt(IND5,:)=[];
NCE_VNCT=NCE_V;
NCE_VNCT(E_VNContestFUL,:)=[];       %NCE_V limited to samples in which there is no entry.
% NCE_VTenCT=Tenure_NCE_V;
% NCE_VTenCT(E_VContestFUL,:)=[];
% NCE_VTenCTwnxt=NCE_VTenCT;
% NCE_VTenCTwnxt(IND5,:)=[];
% NCE_VTenNCT=Tenure_NCE_V;
% NCE_VTenNCT(E_VNContestFUL,:)=[];
LOGTot_NCE_VCT=LOGTot_NCE_V;
LOGTot_NCE_VCT(E_VContestFUL,:)=[];
LOGTot_NCE_VNCT=LOGTot_NCE_V;
LOGTot_NCE_VNCT(E_VNContestFUL,:)=[];


RTotDE_V=(max(0,NCE_V(:,1)).^(-1/2))./LOGTot_NCE_V;
RTotDE_V=RTotDE_V.*(exp(LOGTot_NCE_V)./exp(NCE_V(:,1)));
RTotDE_VCT=(max(0,NCE_VCT(:,1)).^(-1/2))./LOGTot_NCE_VCT;
RTotDE_VCT=RTotDE_VCT.*(exp(LOGTot_NCE_VCT)./exp(NCE_VCT(:,1)));
RTotDE_VCTwnxt=RTotDE_VCT;
RTotDE_VCTwnxt(IND5,:)=[];
RTotDE_VNCT=(max(0,NCE_VNCT(:,1)).^(-1/2))./LOGTot_NCE_VNCT;
RTotDE_VNCT=RTotDE_VNCT.*(exp(LOGTot_NCE_VNCT)./exp(NCE_VNCT(:,1)));


LOGD_E_V=log(max(ones(length(E_V_july8),1),E_V_july8(:,10)));  %log(Spend)
LOGD_E_VCT=LOGD_E_V;
LOGD_E_VCT(E_VContestFUL,:)=[];                             %log(Spend) limited to contested.
LOGD_E_VNCT=LOGD_E_V;
LOGD_E_VNCT(E_VNContestFUL,:)=[];                            %log(Spend) limited to uncontested.



LOGW_E_V=log(max(ones(length(E_V_july8),1),E_V_july8(:,11)));  %log(begining cash)
LOGW_E_VCT=LOGW_E_V;
LOGW_E_VCT(E_VContestFUL,:)=[];
LOGW_E_VCTwnxt=LOGW_E_VCT;
LOGW_E_VCTwnxt(IND5,:)=[];
LOGW_E_VNCT=LOGW_E_V;
LOGW_E_VNCT(E_VNContestFUL,:)=[];


LOGW_NXT_E_V=log(max(ones(length(E_V_july8),1),E_V_july8(:,12)));  %log(realendcash)
LOGW_NXT_E_VCT=LOGW_NXT_E_V;
LOGW_NXT_E_VCT(E_VContestFUL,:)=[];
LOGW_NXT_E_VCTwnxt=LOGW_NXT_E_VCT;
LOGW_NXT_E_VCTwnxt(IND5,:)=[];
LOGW_NXT_E_VNCT=LOGW_NXT_E_V;
LOGW_NXT_E_VNCT(E_VNContestFUL,:)=[];
IND6CT=find(LOGW_NXT_E_VCT==0);                      % IND6CT is the index where realendcash_i ==0
IND6NCT=find(LOGW_NXT_E_VNCT==0);                    % IND6NCT is the index where realendcash_i ==0

LOGW_NXT_E_VC=log(max(ones(length(E_V_july8),1),E_V_july8(:,18)));  %log(realendcash_c)
LOGW_NXT_E_VC(E_VContestFUL,:)=[];
IND7=find(LOGW_NXT_E_VC==0|LOGW_NXT_E_VCT==0);   % IND7 is the index where realedcash_c==0 or realendcash_I==0 (among E_VContestFUL)
LOGD_E_VC=log(max(ones(length(E_V_july8),1),E_V_july8(:,16)));  %log(realdisburse_c)
LOGD_E_VC(E_VContestFUL,:)=[];
LOGTotal_E_V=log(max(ones(length(E_V_july8),1),E_V_july8(:,9)));  %log(realtotal)
LOGTotal_E_VCT=LOGTotal_E_V;
LOGTotal_E_VCT(E_VContestFUL,:)=[];
LOGTotal_E_VNCT=LOGTotal_E_V;
LOGTotal_E_VNCT(E_VNContestFUL,:)=[];
LOGTot_E_VC=log(max(ones(length(E_V_july8),1),E_V_july8(:,4)));  % log(realtotal_c)
LOGTot_E_VC(E_VContestFUL,:)=[];

VSEV=E_V_july8(:,13);
VSEVCT=VSEV;
VSEVCT(E_VContestFUL,:)=[];






%%B-Spline for q_I(RTotDE_V), where RTotDE_V=(NCE_VCT(:,1)./LOGTot_NCE_V)  etc.%%
% take knots to be between 0.88 to 1.025 with 8 knots (8 basis functions). (Almost all lie

 fineness=9;
% mesh=quantile(RTotD_NC,linspace(0,1,fineness).');
mesh=quantile(RTotDE_V,linspace(0,1,fineness).');
X_KnotEV1=(RTotDE_V<mesh(2,1)).*(1-(RTotDE_V-mesh(1,1))/(mesh(2,1)-mesh(1,1)));
for i=0:(numel(mesh)-3)
    PLUS=(RTotDE_V>=mesh(i+1,1)).*(RTotDE_V<mesh(i+2,1)).*((RTotDE_V-mesh(i+1,1))/(mesh(i+2,1)-mesh(i+1,1)))...
        +(RTotDE_V>=mesh(i+2,1)).*(RTotDE_V<mesh(i+3,1)).*(1-(RTotDE_V-mesh(i+2,1))/(mesh(i+3,1)-mesh(i+2,1)));
    X_KnotEV1=[X_KnotEV1,PLUS];
end
PLUS=(RTotDE_V>=mesh((numel(mesh)-1),1)).*(RTotDE_V-mesh((numel(mesh)-1),1))/(mesh(numel(mesh),1)-mesh((numel(mesh)-1),1));
X_KnotEV1=[X_KnotEV1,PLUS];

%First column:redundant
X_KnotEV1=X_KnotEV1(:,2:fineness);
X_KnotE_VCT=X_KnotEV1;
X_KnotE_VCT(E_VContestFUL,:)=[];
X_KnotE_VNCT=X_KnotEV1;
X_KnotE_VNCT(E_VNContestFUL,:)=[];



%Define incumbent quality
   thetaQ2=Est2(9:numel(Est2));
    if numel(thetaQ2)==3
        XQEV2=[RTotDE_V,RTotDE_V.^2,RTotDE_V.^3]*thetaQ2; % q_I
        XQEVCT2=XQEV2;
        XQEVCT2(E_VContestFUL,:)=[];
        XQEVNCT2=XQEV2;
        XQEVNCT2(E_VNContestFUL,:)=[];
        XQEVCTwnxt2=XQEVCT2;
        XQEVCTwnxt2(IND5,:)=[];
    else
        XQEV2=X_KnotEV1*thetaQ2;  % q_I
        XQEVCT2=XQEV2;
        XQEVCT2(E_VContestFUL,:)=[];
        XQEVNCT2=XQEV2;
        XQEVNCT2(E_VNContestFUL,:)=[];
        XQEVCTwnxt2=XQEVCT2;
        XQEVCTwnxt2(IND5,:)=[];
    end

%%%%%%%%%%%%%%%%%%%%%%%% 
% Defining variables: done
%%%%%%%%%%%%%%%%%%%%%%%%

%%

%Continuation payoff: inputed from incumbents' problem
%Value function same as incumbents
Regressand=[ones(length(NCE_V),1),LOGW_NXT_E_V,LOGW_NXT_E_V.^2/10,LOGW_NXT_E_V.^3/100,LOGW_NXT_E_V.^4/1000,...
    XQEV2,LOGW_NXT_E_V.*XQEV2,LOGW_NXT_E_V.^2/10.*XQEV2,LOGW_NXT_E_V.^3/100.*XQEV2,LOGW_NXT_E_V.^4/1000.*XQEV2,...
    TenureE_V,XSEV_(:,1).*SameE_V,XSEV_(:,2).*PartyE_V,PresdumE_V,MidtermE_V];

%With additional intercepts: not used
% Regressand=[ones(length(NCE_V),1),LOGW_NXT_E_V,LOGW_NXT_E_V.^2/10,LOGW_NXT_E_V.^3/100,LOGW_NXT_E_V.^4/1000,...
%     XQEV2,LOGW_NXT_E_V.*XQEV2,LOGW_NXT_E_V.^2/10.*XQEV2,LOGW_NXT_E_V.^3/100.*XQEV2,LOGW_NXT_E_V.^4/1000.*XQEV2,...
%     TenureE_V,TenureE_V.^2,XSEV_(:,1).*SameE_V,XSEV_(:,1).^2.*SameE_V,XSEV_(:,2).*PartyE_V,XSEV_(:,2).^2.*PartyE_V,...
%     TenureE_V.*XSEV_(:,1).*SameE_V,TenureE_V.*XSEV_(:,2).*PartyE_V,PresdumE_V,MidtermE_V];

coef3rdstep=Regressand\Continue;
test=Regressand*coef3rdstep;

Continuec=max(0,[ones(length(LOGW_NXT_E_VC),1),LOGW_NXT_E_VC,LOGW_NXT_E_VC.^2/10,LOGW_NXT_E_VC.^3/100,LOGW_NXT_E_VC.^4/1000,...
    q_C_E_VCT,LOGW_NXT_E_VC.*q_C_E_VCT,LOGW_NXT_E_VC.^2/10.*q_C_E_VCT,LOGW_NXT_E_VC.^3/100.*q_C_E_VCT,LOGW_NXT_E_VC.^4/1000.*q_C_E_VCT,...
    zeros(length(LOGW_NXT_E_VC),1),XS_EVCT_(:,1).*((-1)*SameE_VCT),XS_EVCT_(:,2).*((-1)*PartyE_VCT),...
    PresdumE_VCT,MidtermE_VCT]*coef3rdstep);

%With additional intercepts: not used
%Continuec=max(0,[ones(length(LOGW_NXT_E_VC),1),LOGW_NXT_E_VC,LOGW_NXT_E_VC.^2/10,LOGW_NXT_E_VC.^3/100,LOGW_NXT_E_VC.^4/1000,...
%    q_C_E_VCT,LOGW_NXT_E_VC.*q_C_E_VCT,LOGW_NXT_E_VC.^2/10.*q_C_E_VCT,LOGW_NXT_E_VC.^3/100.*q_C_E_VCT,LOGW_NXT_E_VC.^4/1000.*q_C_E_VCT,...
%    zeros(length(LOGW_NXT_E_VC),1),zeros(length(LOGW_NXT_E_VC),1),XS_EVCT_(:,1).*((-1)*SameE_VCT),XS_EVCT_(:,1).^2.*((-1)*SameE_VCT),XS_EVCT_(:,2).*((-1)*PartyE_VCT),XS_EVCT_(:,2).^2.*((-1)*PartyE_VCT),...
%    zeros(length(LOGW_NXT_E_VC),2),PresdumE_VCT,MidtermE_VCT]*coef3rdstep);


%%

%Specification check
figure(1)
scatter(q_C_E_VCT(LOGW_NXT_E_VC>0&Continuec>0),Continuec(LOGW_NXT_E_VC>0&Continuec>0))
figure(2)
scatter(q_C_E_VCT(LOGW_NXT_E_VC>0&Continuec<=0),Continuec(LOGW_NXT_E_VC>0&Continuec<=0))
figure(3)
scatter(q_C_E_VCT(LOGW_NXT_E_VC==0),Continuec(LOGW_NXT_E_VC==0))
%scatter(q_C_E_VCT(LOGW_NXT_E_VC>0&Continuec>0),Deriv(LOGW_NXT_E_VC>0&Continuec>0))
%scatter(q_C_E_VCT(Deriv>0),Deriv(Deriv>0))
length(find(LOGW_NXT_E_VC>0))
length(find(Continuec>0))
%length(find(Deriv>0))
% 
% D=@(test) coef3rdstep(6)+test*coef3rdstep(7)+test.^2/10*coef3rdstep(8)+test.^3/100*coef3rdstep(9)+test.^4/1000*coef3rdstep(10);
% space1=9:0.5:15;
% mat1=D(space1);
% scatter(space1,mat1);




%% use LOGW_EVCT, etc, variables for incumbents %%
% 
% PrE=normcdf([ones(size(LOGW_E_VCT,1),1),LOGW_E_VCT,XQEVCT,XS_EVCT.*PartyE_VCT,log(TenureE_VCT+1)]*theta,0,1);
% PrE3step=PrE;
% PrE3step(IND6CT,:)=[];
% PrE3step(DEL_Pen,:)=[];
% E_PrimaryN=[ones(size(LOGW_E_VCT,1),1),LOGW_E_VCT,XQEVCT,XS_EVCT.*PartyE_VCT,log(TenureE_VCT+1)]*theta2;
% E_PrimaryN3step=E_PrimaryN;
% E_PrimaryN3step(IND6CT,:)=[];
% E_PrimaryN3step(DEL_Pen,:)=[];
% 
% 

BX1=abs(1/mintheta2ndstage(3))*([LOGD_E_VCT,LOGD_E_VC,XS_EVCT_.*[SameE_VCT,PartyE_VCT],log(TenureE_VCT+1)]*Est2(1:5)+XQEVCT2+q_C_E_VCT);
%BX13step=BX1;
%BX13step(IND6CT,:)=[];
%BX13step(DEL_Pen,:)=[];
% thetain2=Sndstage_kettei_long;
% % Estim3=[   -1.9713322e+01
% %   -1.1704758e-02
% %    1.6123112e+03
% %   -7.4047131e+033]

datasetVC=[LOGTot_E_VC,LOGW_NXT_E_VC,LOGD_E_VC,q_C_E_VCT,BX1,Continuec];
%datasetVC(Continuec<max(Continue)|(LOGW_NXT_E_VC==0),:)=[];
%datasetVC((LOGW_NXT_E_VC==0),:)=[];
datasetVC(((LOGW_NXT_E_VC==0)&(Continuec>100)),:)=[];
datasetVC(((LOGW_NXT_E_VC==0)|(q_C_E_VCT>0.2)),:)=[];
%datasetVC(Continuec>100,:)=[];
%%
para0=mintheta2ndstage([1,4,5,6]);
 options=optimset('MaxIter',200000,'MaxFunEvals',1000000,'Display','iter');
[mintheta,SRR]=fminsearch(@(theta3) Minimize3PS_long(Est2,mintheta2ndstage,theta3,datasetVC),para0,options);


est3rdstage=mintheta;
%save est3rdstage2 est3rdstage
%%
%Check results
%load est3rdstage2
alpha=1/2;
beta=2;
ben2=abs(mintheta2ndstage(2));
cost1=abs(mintheta2ndstage(1));
costc=abs(est3rdstage(1));
a=est3rdstage(2);
b=est3rdstage(3);
c=est3rdstage(4);
aa=mintheta2ndstage(4);
bb=mintheta2ndstage(5);
cc=mintheta2ndstage(6);
%Cost function
costchal=@(test,test2)costc*(alpha/beta)*ben2*(1+a*test+b*test.^2+c*test.^3).*test2.^beta;
costinc=@(test,test2)cost1*(alpha/beta)*ben2*(1+aa*test+bb*test.^2+cc*test.^3).*test2.^beta;

%Specification 1
% space1=0.001:0.01:0.201;
% matchal1=costchal(space1,10);
% matinc1=costinc(space1,10);
%Specification 2
space1=0:0.01:0.2;
space2=0.1:0.01:0.5;
space3=-0.2:0.01:0.2;
matchal1=costchal(space2,10);
matinc1=costinc(space1,10);

scatter(space3,matchal1)
hold on
scatter(space1,matinc1);
%%

load Best3step.txt
mintheta=Best3step(1:4,1);
Vandq_c=Vandq_c(thetain1step,thetain2,mintheta);
[Continue_IC,Continue_INC]=V_Iandq_I(thetain1step,thetain2);
Continue_IC=[Continue_IC,XQEV2];
Continue_INC=[Continue_INC,XQEV2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:size(PrE3step)
    [result(i,1), r(i,1)]=fminunc(@(r) cutoff(r,E_PrimaryN3step(i,1),PrE3step(i,1)),2); % result() is the parameter of the negative binomial from which draw the # of challengers
end
save result.mat result


load result
result=abs(result);
partition_t=5;
partition_v=4;
QQ=NaN(partition_t+1,partition_v+1,floor(size(result,1)/((partition_t+1)*(partition_v+1)))+3);
t_i=quantile(result,[1:partition_t]'/(partition_t+1));
t_i=[min(result);t_i;max(result)];
s_is=[];
beta_Vq_=[];
for i=1:partition_t+1
    I=find( (result>=t_i(i,1)).*(result<=t_i(i+1,1))==1);
    minV=Vandq_c(I,:);
    PrEI=PrE3step(I,1);
    s_i_=quantile(PrEI,[1:partition_v]'/(partition_v+1));
    s_is_=[quantile(PrEI,1/(2*partition_v));s_i_;quantile(PrEI,1-1/(2*partition_v))];
    s_i_=[min(PrEI)-1;s_i_;max(PrEI)+1];
    for j=1:partition_v+1
        J=find( (PrEI>s_i_(j,1)).*(PrEI<s_i_(j+1,1))==1);
        minVV=minV(J,:);
        [C_(j,1),wx]=min(minVV(:,2));
        C_(j,2)=minVV(wx,1);
        %C_(j,2)=min(minVV(:,1));
        QQ(i,j,1:size(minVV,1))=sort(minVV(:,2));
    end
    CQ(i,:,:)=sortrows(C_,1);
    C__=[ones(partition_v+1,1),C_(:,1)];
    beta_cq(i,1:2)=(inv(C__'*C__)*C__'*C_(:,2))';
    s_is=[s_is,s_is_];
end
t_is=quantile(result,[1:partition_t]'/(partition_t+1));
t_is=[quantile(result, 1/(2*partition_t));t_is;quantile(result, 1-1/(2*partition_t))];

save beta_cq.mat beta_cq QQ t_is s_is Continue_IC Continue_INC


% 
% s_j=quantile(minuu*beta_Vq(1:5,1),[1:6]'/7);
% J=find((minuu*beta_Vq(1:5,1)<s_j(1,1))==1);
% minVJ=minV(J,:);
% [underbar_V(1,1),ii]=min(minVJ(:,1));
% underbar_q(1,1)=minVJ(ii,2);
% distr=sort(minVJ(:,2));
% for j=2:6
%     J=find((minuu*beta_Vq(1:5,1)<s_j(j,1)).*(minuu*beta_Vq(1:5,1)>s_j(j-1,1))==1);
%     minVJ=minV(J,:);
%     [unberbar_V(1,j),ii]=min(minVJ(:,1));
%     underbar_q(1,j)=minVJ(ii,2);
%     if size(minVJ,1)<size(distr,1)
%         minVJ=[minVJ;[NaN,NaN]];
%     end
%     distr=[distr,sort(minVJ(:,2))];
% end
% 
% for i=2:6
%     I=find((result<t_i(i,1)).*(result>t_i(i-1,1))==1);
%     minuu=uu(I,1:5);
%     minV=Vandq_c(I,:);   
%     s_j=quantile(minuu*beta_Vq(1:5,1),[1:6]'/7);
%     J=find((minuu*beta_Vq(1:5,1)<s_j(1,1))==1);
%     minVJ=minV(J,:);
%     [underbar_V(i,1),ii]=min(minVJ(:,1));
%     underbar_q(i,1)=minVJ(ii,2);
%     if size(minVJ,1)<size(distr,1)
%         minVJ=[minVJ;[NaN,NaN]];
%     end
%     distr=[distr,sort(minVJ(:,2))];
%     
%     for j=2:6
%         J=find((minuu*beta_Vq(1:5,1)<s_j(j,1)).*(minuu*beta_Vq(1:5,1)>s_j(j-1,1))==1);
%         minVJ=minV(J,:);
%         [unberbar_V(i,j),ii]=min(minVJ(:,1));
%         underbar_q(i,j)=minVJ(ii,2);
%         if size(minVJ,1)<size(distr,1)
%             minVJ=[minVJ;NaN(size(distr,1)-size(minVJ,1),2)];
%         end
%         if size(minVJ,1)>size(distr,1)
%             distr=[distr;NaN(size(minVJ,1)-size(distr,1),size(distr,2))];
%         end
%     distr=[distr,sort(minVJ(:,2))];
%     end
%     J=find((minuu*beta_Vq(1:5,1)>s_j(6,1))==1);
%     minVJ=minV(J,:);
%     [underbar_V(i,7),ii]=min(minVJ(:,1));
%     underbar_q(i,7)=minVJ(ii,2);
%     if size(minVJ,1)<size(distr,1)
%         minVJ=[minVJ;NaN(size(distr,1)-size(minVJ,1),2)];
%     end
%     distr=[distr,sort(minVJ(:,2))];
% end
% 
% I=find((result>t_i(6,1))==1);
% minuu=uu(I,1:5);
% minV=Vandq_c(I,:);
% s_j=quantile(minuu*beta_Vq(1:5,1),[1:6]'/7);
% J=find((minuu*beta_Vq(1:5,1)<s_j(1,1))==1);
% minVJ=minV(J,:);
% [underbar_V(7,1),ii]=min(minVJ(:,1));
% underbar_q(7,1)=minVJ(ii,2);
% if size(minVJ,1)<size(distr,1)
%    minVJ=[minVJ;NaN(size(distr,1)-size(minVJ,1),2)];
% end
% distr=[distr,sort(minVJ(:,2))];
% for j=2:6
%     J=find((minuu*beta_Vq(1:5,1)<s_j(j,1)).*(minuu*beta_Vq(1:5,1)>s_j(j-1,1))==1);
%     minVJ=minV(J,:);
%     [unberbar_V(7,j),ii]=min(minVJ(:,1));
%     underbar_q(7,j)=minVJ(ii,2);
%     if size(minVJ,1)<size(distr,1)
%         minVJ=[minVJ;NaN(size(distr,1)-size(minVJ,1),2)];
%     end
%     distr=[distr,sort(minVJ(:,2))];
% end
% J=find((minuu*beta_Vq(1:5,1)>s_j(6,1))==1);
% minVJ=minV(J,:);
% [underbar_V(7,7),ii]=min(minVJ(:,1));
% underbar_q(7,j)=minVJ(ii,2);
% if size(minVJ,1)<size(distr,1)
%         minVJ=[minVJ;NaN(size(distr,1)-size(minVJ,1),2)];
% end
% distr=[distr,sort(minVJ(:,2))];
% 
% 
% 
% 
% 
% 
% minV=Vandq_c(I,1);
% underbar_V(21,1)=min(minV);
% save Vandq_c.mat Vandq_c t_i underbar_V
% 
% [a,b,c,d]=kde2d(Vandq_c,128);
% bwidth=0.2;
% gridsize=100;
% range=max(Vandq_c(:,1))-min(Vandq_c(:,1));
% q_cgrid=[1:100]*(0.8/100)-0.5;
% q_cgrid=q_cgrid';
% f=[];
%     for i=1:gridsize
%         Vgrid(i,1)=min(Vandq_c(:,1))+(range/(gridsize+1))*i;
%         Vandq_c_i=Vandq_c((((Vandq_c(:,1))>Vgrid(i,1)-bwidth)&((Vandq_c(:,1))<Vgrid(i,1)+bwidth)),2);
%         f_i=ksdensity(Vandq_c_i,q_cgrid);
%         f=[f,f_i];
%     end
%     
% %%%% Temp %%%%%%%%
% clear VS
% Cap=2000;
% chp=0.1;
% PrE1=normcdf([1,mean(LOGW_E_VCT),chp,mean(XS_EVCT.*PartyE_VCT),mean(log(TenureE_VCT+1))]*theta,0,1);
% E_PrimaryN1=[1,mean(LOGW_E_VCT),chp,mean(XS_EVCT.*PartyE_VCT),mean(log(TenureE_VCT+1))]*theta2;
% N_E1=monotone_exp(PrE1,E_PrimaryN1);
% F_cutoff1=(1-PrE1)^(1/N_E1);
% U=rand(Cap,N_E1);
% X=max(U,[],2);
% DX=find(X<F_cutoff1);
% X(DX,:)=[];
% size(X)
% X=max(U,[],2);
% DX=find(X<F_cutoff1);
% X(DX,:)=[];
% X=qc_LB+(qc_UB-qc_LB)*icdf('beta',X,Trdstage_kettei_short(3,1),Trdstage_kettei_short(4,1));
% Z=X;
% 
% 
% epsi=randn(Cap,1)*Sndstage_kettei_short(5,1);
% for j=1:size(X,1)
%     [val,id]=min(abs(XQEVCT-q_C_E_VCT+XS_EVCT.*PartyE_VCT-(chp-X(j,1)+mean(XS_EVCT.*PartyE_VCT))));
%     VS(j,1)=B_I*LOGD_E_VCT(id,1)+B_C*LOGD_E_VC(id,1)+mean(XS_EVCT.*PartyE_VCT)-X(j,1)+XQEVCT(id,1)+epsi(j,1);
% end
% DX=(VS>0);
% XA=X;
% XA(DX,:)=[];
% 
% for j=1:size(XA,1)
%     [val,id]=min(abs(XA(j)-q_c3step));
%     k=ceil(rand(1,1)*N);
%     incad(1:5,:,j)=Cc(:,:,k,id);
%     incad(6,:,j)=XA(j);
% end
% X1=squeeze(incad(6,1,:));
% DX=find(incad(5,1,:)==0);
% X1(DX,:)=[];
% X3=squeeze(incad(6,3,:));
% DX=find(incad(5,3,:)==0);
% X3(DX,:)=[];
% X5=squeeze(incad(6,5,:));
% DX=find(incad(5,5,:)==0);
% X5(DX,:)=[];
% X7=squeeze(incad(6,7,:));
% DX=find(incad(5,7,:)==0);
% X7(DX,:)=[];
% XAN=quantile(XA,rand(Cap,1));
% X3N=quantile(X3,rand(Cap,1));
% X5N=quantile(X5,rand(Cap,1));
% X7N=quantile(X7,rand(Cap,1));
% subplot(5,1,1)
% hist(X)
% subplot(5,1,2)
% hist(XAN)
% subplot(5,1,3)
% hist(X3N)
% subplot(5,1,4)
% hist(X5N)
% subplot(5,1,5)
% hist(X7N)
