clear all
% **
% make sure you get the sign of q_e right. Also note that f_d_qe, f_w_qe have the PARTY of Challenger.
% **

rng(958)

%Globals

%Calculating C matrix
%global N
%global interest



%load ('./E_V_july8.mat')
E_V_july8=csvread('E_Vjuly8partisan.csv',1,0);

% load('./op_inc_july8.mat')
%OP_INC_july8=csvread('OP_INC_july8partisan.csv',1,0);

% load('./op_inc_iv_july8.mat')
%OP_INC_IV_july8=csvread('OP_INC_IV_july8partisan.csv',1,0);

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

load ('./retire.txt');

% load ('./theta0.txt')
load ('./stateevol.mat')
% 
% load ('./Fststage_kettei.txt')
 load ('./Sndstage_initial_long.txt')
 


% rand('state',1000);
% randn('state',10);

% dele1=find(sum(isnan(OP_INC_IV_july8),2)>0);
% OP_INC_IV_july8(dele1,:)=[]; %Drop NaN
% OP_INC_IV_july8(OP_INC_IV_july8(:,2)>2002,:)=[]; %Drop year 2004 and on
% deleA=find(OP_INC_IV_july8(:,16)==0);
% OP_INC_IV_july8(deleA,:)=[];%Drop if oppornent disburse=0
% OP_INC_IV_july8(102,:)=[];%Drop an Rtotd outlier
% deleA2=find(OP_INC_IV_july8(:,12)<1);
% OP_INC_IV_july8(deleA2,:)=[];%Drop if nonpositive saving
% deleA3=find(OP_INC_IV_july8(:,11)<1);
% OP_INC_IV_july8(deleA3,:)=[];%Drop if nonpositive begcash

% dele2=find(sum(isnan(OP_INC_july8),2)>0);
% OP_INC_july8(dele2,:)=[]; %Drop NaN
% OP_INC_july8(OP_INC_july8(:,2)>2002,:)=[]; %Drop year 2004 and on
% deleB=find(((OP_INC_july8(:,16)==0).*OP_INC_july8(:,8))==1);
% OP_INC_july8(deleB,:)=[];%Drop if oppornent disburse=0 and contested
% OP_INC_july8(119:121,:)=[];%Drop an Rtotd outlier
% deleB2=find(OP_INC_july8(:,12)<1);
% OP_INC_july8(deleB2,:)=[];%Drop if nonpositive saving
% deleB3=find(OP_INC_july8(:,11)<1);
% OP_INC_july8(deleB3,:)=[];%Drop if nonpositive begcash

dele3=find(sum(isnan(E_V_july8),2)>0); 
E_V_july8(dele3,:)=[]; %Drop NaN
E_V_july8(E_V_july8(:,2)>2002,:)=[]; %Drop year 2004 and on
deleC=find((E_V_july8(:,16)==0).*E_V_july8(:,8)==1);
E_V_july8(deleC,:)=[]; %Drop if oppornent disburse=0 and contested
E_V_july8(174:177,:)=[]; %Drop an Rtotd outlier
deleC2=find(E_V_july8(:,12)<5000);
E_V_july8(deleC2,:)=[];%Drop if nonpositive saving
deleC3=find(E_V_july8(:,11)<1);
E_V_july8(deleC3,:)=[];%Drop if nonpositive begcash
deleC4=find(E_V_july8(:,16)<5000&E_V_july8(:,8)==1);
E_V_july8(deleC4,:)=[];%Drop if nonpositive challenger spending

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
%% take knots to be between 0.88 to 1.025 with 8 knots (8 basis functions). (Almost all lie
% within this range) Let X_Knot be the matrix with 8 columns that contain
% the value of the B-Spline basis function evaluated at each of the 8
% basis functions.
% 
% mesh=quantile(LOGLOGD_E_V,[.125;.25;.375;.5;.625;.75;.875]);
% mesh=[min(LOGLOGD_E_V);mesh;max(LOGLOGD_E_V)];
% X_KnotEV1=(LOGLOGD_E_V<mesh(2,1)).*(1-(LOGLOGD_E_V-mesh(1,1))/(mesh(2,1)-mesh(1,1)));
% for i=1:6
%     PLUS=(LOGLOGD_E_V>=mesh(i+1,1)).*(LOGLOGD_E_V<mesh(i+2,1)).*((LOGLOGD_E_V-mesh(i+1,1))/(mesh(i+2,1)-mesh(i+1,1)))...
%         +(LOGLOGD_E_V>mesh(i+2,1)).*(LOGLOGD_E_V<mesh(i+3,1)).*(1-(LOGLOGD_E_V-mesh(i+2,1))/(mesh(i+3,1)-mesh(i+2,1)));
%     X_KnotEV1=[X_KnotEV1,PLUS];
% end
% PLUS=(LOGLOGD_E_V>=mesh(8,1)).*(LOGLOGD_E_V-mesh(8,1))/(mesh(9,1)-mesh(8,1));
% X_KnotEV1=[X_KnotEV1,PLUS];
% 
% mesh=quantile(LOGLOGTot_E_V,[.125;.25;.375;.5;.625;.75;.875]);
% mesh=[min(LOGLOGTot_E_V);mesh;max(LOGLOGTot_E_V)];
% X_KnotEV2=(LOGLOGTot_E_V<mesh(2,1)).*(1-(LOGLOGTot_E_V-mesh(1,1))/(mesh(2,1)-mesh(1,1)));
% for i=1:6
%     PLUS=(LOGLOGTot_E_V>=mesh(i+1,1)).*(LOGLOGTot_E_V<mesh(i+2,1)).*((LOGLOGTot_E_V-mesh(i+1,1))/(mesh(i+2,1)-mesh(i+1,1)))...
%         +(LOGLOGTot_E_V>mesh(i+2,1)).*(LOGLOGTot_E_V<mesh(i+3,1)).*(1-(LOGLOGTot_E_V-mesh(i+2,1))/(mesh(i+3,1)-mesh(i+2,1)));
%     X_KnotEV2=[X_KnotEV2,PLUS];
% end
% PLUS=(LOGLOGTot_E_V>=mesh(8,1)).*(LOGLOGTot_E_V-mesh(8,1))/(mesh(9,1)-mesh(8,1));
% X_KnotEV2=[X_KnotEV2,PLUS];
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

%%%%%%%%%%%%%%%%%%%%%%%% 
% Defining variables: done
%%%%%%%%%%%%%%%%%%%%%%%%




%Estimated parameters in the first stage
coefentry=Est1(1:(length(Est1)/2));
%Saving and fundraising parameters when contested
%     E_VCTa(:,1)=[Est411(1:8);Est421(1:8);Est431(1:8)];
%     E_VCTa(:,2)=[Est412(1:8);Est422(1:8);Est432(1:8)];
%     E_VCTa(:,3)=[Est413(1:8);Est423(1:8);Est433(1:8)];
%     E_VCTt(:,1)=[Est411(9:16);Est421(9:16);Est431(9:16)];
%     E_VCTt(:,2)=[Est412(9:16);Est422(9:16);Est432(9:16)];
%     E_VCTt(:,3)=[Est413(9:16);Est423(9:16);Est433(9:16)];
%     gammaCT(:,1)=[Est411(17:24);Est421(17:24);Est431(17:24)];
%     gammaCT(:,2)=[Est412(17:24);Est422(17:24);Est432(17:24)];
%     gammaCT(:,3)=[Est413(17:24);Est423(17:24);Est433(17:24)];
    E_VCTa(:,1)=[Est411(1:8);Est421(1:8)];
    E_VCTa(:,2)=[Est412(1:8);Est422(1:8)];
    E_VCTa(:,3)=[Est413(1:8);Est423(1:8)];
    E_VCTt(:,1)=[Est411(9:16);Est421(9:16)];
    E_VCTt(:,2)=[Est412(9:16);Est422(9:16)];
    E_VCTt(:,3)=[Est413(9:16);Est423(9:16)];
    gammaCT(:,1)=[Est411(17:24);Est421(17:24)];
    gammaCT(:,2)=[Est412(17:24);Est422(17:24)];
    gammaCT(:,3)=[Est413(17:24);Est423(17:24)];
%Saving and fundraising parameters conditional on contested and winning
    E_VCTa(:,4)=[Est4311(1:8);Est4321(1:8)];
    E_VCTa(:,5)=[Est4312(1:8);Est4322(1:8)];
    E_VCTa(:,6)=[Est4313(1:8);Est4323(1:8)];
    E_VCTt(:,4)=[Est4311(9:16);Est4321(9:16)];
    E_VCTt(:,5)=[Est4312(9:16);Est4322(9:16)];
    E_VCTt(:,6)=[Est4313(9:16);Est4323(9:16)];
    gammaCT(:,4)=[Est4311(17:24);Est4321(17:24)];
    gammaCT(:,5)=[Est4312(17:24);Est4322(17:24)];
    gammaCT(:,6)=[Est4313(17:24);Est4323(17:24)];


coefspend(:,1)=Est51;
coeffund(:,1)=Est52;
coefsave(:,1)=Est53;
coefprobwin=Est3;

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
   
   
% iterate=1;
% % XS2=[Unemployment,White]*thetaS2;
% % XSEV2=XSEV_*thetaS2;
% thetaQ(1,1)=1;
% %     for ii=2:size(X_Knot1,2)
% %         thetaQ(ii,1)=1-sum(abs(Est1(3:(ii+1),1)));
% %     end
% thetaQ(2:9,1)=thetain1step(3:10,1);
% thetaQ2(1,1)=thetain1step(24,1);
% %     for ii=2:size(X_Knot1,2)
% %         thetaQ2(ii,1)=thetaQ2(1,1)-sum(abs(thetain(4:(ii+2),1)));
% %     end
% thetaQ2(2:9,1)=thetain1step(25:32,1);
% 
% XQEV=X_KnotEV1*thetaQ;  % q_I
% XQEVCT=XQEV;
% XQEVCT(E_VContestFUL,:)=[];
% XQEVNCT=XQEV;
% XQEVNCT(E_VNContestFUL,:)=[];
% 
% XQEV2=X_KnotEV1*thetaQ2;  % q_I
% XQEVCT2=XQEV2;
% XQEVCT2(E_VContestFUL,:)=[];
% XQEVNCT2=XQEV2;
% XQEVNCT2(E_VNContestFUL,:)=[];
% save XQEVNCT2.mat XQEVNCT2
% save XQEVCT2.mat XQEVCT2

%Result from the first stage till here.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



Sofarbest=10^8;
bestiter=0;
iterate=1;
N=100;                      % N is the number of simulations.
NumSim=50;              %NumSim is the number of simulations in the 1st step to obtain f(.,q_e)

T=10;                      % T is the number of periods that we move the simlation forward.
interest=0.1;               %interest on money saved.
Entry=rand(T,N,length(E_V_july8));             %Simulation draw for computing the continuation value E_V
Winrnd=rand(T,N,length(E_V_july8));            %Simulation draw for computing the continuation value E_V (for challenger FOC)

Retirernd=(rand(T-1,N,length(E_V_july8))<repmat(retire',[1 N, length(E_V_july8)]));
Retirernd(T,:,:)=1;
RetirerndC=(rand(T-1,N,length(E_V_july8))<repmat(retire',[1 N, length(E_V_july8)]));
RetirerndC(T,:,:)=1;
RetirerndR=(rand(T-1,N,length(E_V_july8))<repmat(retire',[1 N, length(E_V_july8)]));
RetirerndR(T,:,:)=1;


for i=1:N
    for j=1:length(E_V_july8)
        Ret(i,j)=find(Retirernd(:,i,j)==1,1,'first');
    end
end
for i=1:N
    for j=1:length(E_V_july8)
        RetC(i,j)=find(RetirerndC(:,i,j)==1,1,'first');
    end
end
for i=1:N
    for j=1:length(E_V_july8)
        RetR(i,j)=find(RetirerndR(:,i,j)==1,1,'first');
    end
end


% minmax=min(max(Ret,[],1))
% 
% minmax=min(max(RetC,[],1))
% 
% minmax=min(max(RetR,[],1))


dF_gamma_ct=rand(T,N,length(E_V_july8));
dF_total_ct=rand(T,N,length(E_V_july8));
dF_gammasv_ct=rand(T,N,length(E_V_july8));
dF_totalsv_ct=rand(T,N,length(E_V_july8));
dF_nxt_nxt_ct=rand(T,N,length(E_V_july8));

% dF_gamma_nct=rand(T,N,length(E_V_july8));
% dF_total_nct=rand(T,N,length(E_V_july8));
% dF_nxt_nxt_nct=rand(T,N,length(E_V_july8));

dF_gamma_ctC=rand(T,N,length(E_V_july8));
dF_total_ctC=rand(T,N,length(E_V_july8));
dF_nxt_nxt_ctC=rand(T,N,length(E_V_july8));

% dF_gamma_nctC=rand(T,N,length(E_V_july8));
% dF_total_nctC=rand(T,N,length(E_V_july8));
% dF_nxt_nxt_nctC=rand(T,N,length(E_V_july8));

dF_gamma_ctR=rand(T,N,length(E_V_july8));
dF_total_ctR=rand(T,N,length(E_V_july8));
dF_nxt_nxt_ctR=rand(T,N,length(E_V_july8));

% dF_gamma_nctR=rand(T,N,length(E_V_july8));
% dF_total_nctR=rand(T,N,length(E_V_july8));
% dF_nxt_nxt_ncR=rand(T,N,length(E_V_july8));



%Draw evolution of X
%Parameters governing evolution of partisanship already adjusted to
%two-year window.
Shockump=sqrt(epsump)*randn(T,N,length(E_V_july8));
Shockpartisan=sqrt(epspartisan)*randn(T,N,length(E_V_july8));

Shockump2=sqrt(epsump)*randn(T,N,length(E_V_july8));
Shockpartisan2=sqrt(epspartisan)*randn(T,N,length(E_V_july8));

Shockump3=sqrt(epsump)*randn(T,N,length(E_V_july8));
Shockpartisan3=sqrt(epspartisan)*randn(T,N,length(E_V_july8));

%%
%Draw evolution of president incumbency. i.e. probability of president
%party switching
incwin=0.75;
newele=0.5;

%Midone considers a case:
%Presidential election takes place today, in which incumbent is running for reelection. 
%This implies the next election (initial element of the matrix) corresponds
%to mid-term election where president might be either new/old.
midone=find(YEARE_V==1984|YEARE_V==1992|YEARE_V==1996|YEARE_V==2004);


%Draws for winning.
presseqmidone=rand(T,N,length(midone));

%I need to keep track of whether the current president can run for
%reelection or not.
presyearmidone=ones(T,N,length(midone));

%Elections overlapping with presidential election (even periods): previous election being
%midterm. So no change in presidential incumbency from previous period.
presseqmidone([2,4,6,8,10],:,:)=0;
%The next election from now=midterm election. In the previous election (today),
%incumbent runs for reelection. 1 if incumbent loses, 0 otherwise.
presseqmidone(1,:,:)=(presseqmidone(1,:,:)>incwin);
presyearmidone(1,:,:)=1+(1-presseqmidone(1,:,:));
%At election i, if incumbent year=2,  new election. If year=1,
%election with incumbency.
for k=1:length(midone)
    for n=1:N
        for i=[1,3,5,7]
            if presyearmidone(i,n,k)==1
                presseqmidone((i+2),n,k)=(presseqmidone((i+2),n,k)>incwin);
                presyearmidone((1+2),n,k)=presyearmidone(i,n,k)+(1-presseqmidone(i+2,n,k));
            else
                presseqmidone((i+2),n,k)=(presseqmidone((i+2),n,k)>newele);
                presyearmidone((1+2),n,k)=1;
            end
        end
    end
end


%Midtwo considers a case:
%Presidential election takes place today, in which incumbent cannot run for reelection. 
%This implies the next election (initial element of the matrix) corresponds
%to mid-term election where president must be new.
midtwo=find(YEARE_V==1988|YEARE_V==2000);


%Draws for winning.
presseqmidtwo=rand(T,N,length(midtwo));

%I need to keep track of whether the current president can run for
%reelection or not.
presyearmidtwo=ones(T,N,length(midtwo));

%Midterm election: No change in incumbency
presseqmidtwo([2,4,6,8,10],:,:)=0;
%Today:incumbent cannot run for reelection. 1 if party changes.
presseqmidtwo(1,:,:)=(presseqmidtwo(1,:,:)>newele);
presyearmidtwo(1,:,:)=1;
%At election i, if incumbent year=2,  new election. If year=1,
%election with incumbency.
for k=1:length(midtwo)
    for n=1:N
        for i=[1,3,5,7]
            if presyearmidtwo(i,n,k)==1
                presseqmidtwo((i+2),n,k)=(presseqmidtwo((i+2),n,k)>incwin);
                presyearmidtwo((1+2),n,k)=presyearmidtwo(i,n,k)+(1-presseqmidtwo(i+2,n,k));
            else
                presseqmidtwo((i+2),n,k)=(presseqmidtwo((i+2),n,k)>newele);
                presyearmidtwo((1+2),n,k)=1;
            end
        end
    end
end

%Begone considers a case:
%Today:midterm election.  
%and Next election has presidential, where incumbent can run for reelection.
begone=find(YEARE_V==1990|YEARE_V==1994|YEARE_V==2002);

%Draws for winning.
presseqbegone=rand(T,N,length(begone));

%I need to keep track of whether the current president can run for
%reelection or not.
presyearbegone=ones(T,N,length(begone));

%If previous election is midterm election: No change in incumbency
presseqbegone([1,3,5,7,9],:,:)=0;
%2nd election=incumbent runs for reelection at 1st election. 1 if incumbent loses, 0 otherwise.
presseqbegone(2,:,:)=(presseqbegone(2,:,:)>incwin);
presyearbegone(2,:,:)=1+(1-presseqbegone(2,:,:));
%At election i, if incumbent year=2,  new election. If year=1,
%election with incumbency.
for k=1:length(begone)
    for n=1:N
        for i=[2,4,6,8]
            if presyearbegone(i,n,k)==1
                presseqbegone((i+2),n,k)=(presseqbegone((i+2),n,k)>incwin);
                presyearbegone((1+2),n,k)=presyearbegone(i,n,k)+(1-presseqbegone(i+2,n,k));
            else
                presseqbegone((i+2),n,k)=(presseqbegone((i+2),n,k)>newele);
                presyearbegone((1+2),n,k)=1;
            end
        end
    end
end


%Begtwo considers a case:
%Today:midterm election.  
%and Next election has presidential, where incumbent cannot run for reelection.
begtwo=find(YEARE_V==1986|YEARE_V==1998);

%Draws for winning.
presseqbegtwo=rand(T,N,length(begtwo));

%I need to keep track of whether the current president can run for
%reelection or not.
presyearbegtwo=ones(T,N,length(begtwo));

%begterm election: No change in incumbency
presseqbegtwo([1,3,5,7,9],:,:)=0;
%2nd election=incumbent cannot run for reelection. 1 if party changes.
presseqbegtwo(2,:,:)=(presseqbegtwo(2,:,:)>newele);
presyearbegtwo(2,:,:)=1;
%At election i, if incumbent year=2,  new election. If year=1,
%election with incumbency.
for k=1:length(begtwo)
    for n=1:N
        for i=[2,4,6,8]
            if presyearbegtwo(i,n,k)==1
                presseqbegtwo((i+2),n,k)=(presseqbegtwo((i+2),n,k)>incwin);
                presyearbegtwo((1+2),n,k)=presyearbegtwo(i,n,k)+(1-presseqbegtwo(i+2,n,k));
            else
                presseqbegtwo((i+2),n,k)=(presseqbegtwo((i+2),n,k)>newele);
                presyearbegtwo((1+2),n,k)=1;
            end
        end
    end
end

%Finally, aggregate them all to make one big matrix.
presseq=zeros(T,N,length(E_V_july8));
presseq(:,:,midone)=presseqmidone;
presseq(:,:,midtwo)=presseqmidtwo;
presseq(:,:,begone)=presseqbegone;
presseq(:,:,begtwo)=presseqbegtwo;



%%

%Continue1=zeros(length(NCE_V),N);
%DContinue1=zeros(length(NCE_V),N);
% C=zeros(5,T,N,length(NCE_V));
% %DC=zeros(5,T,N,length(NCE_V));
% parfor i=1:length(NCE_V)
%          C(:,:,:,i)=Actions(i,XSEV_(i,:)',SameE_V(i,1), XQEV2(i,1),X_KnotEV1(i,:),RTotDE_V(i,:), TenureE_V(i,1),LOGW_NXT_E_V(i,1),PartyE_V(i,1),PresdumE_V(i,1),MidtermE_V(i,1),...
%          Shockpartisan(:,:,i), Shockump(:,:,i), presseq(:,:,i),Entry(:,:,i), coefentry, E_VCTa, E_VCTt, gammaCT, coefspend, coeffund,coefsave, dF_gamma_ct(:,:,i), dF_total_ct(:,:,i),...
%        dF_gammasv_ct(:,:,i),dF_totalsv_ct(:,:,i),Winrnd(:,:,i),coefprobwin,Ret(:,i)',Betapartisan,Betaump,N,interest);
%     squeeze(C(:,:,:,i));
%     i
% %     DC(:,:,:,i)=Actions(i,XSEV_(i,:)',XQEV(i,1),RTotDE_V(i,1),TenureE_V(i,1),LOGW_NXT_E_V(i,1)+Delt,...
% %         Shockwh(:,:,i), Shockump(:,:,i), Entry(:,:,i), theta, E_VCTa, E_VCTt, gammaCT, E_VNCTa, E_VNCTt, gammaNCT, dF_gamma_ct(:,:,i), dF_total_ct(:,:,i),...
% %         dF_nxt_nxt_ct(:,:,i),dF_gamma_nct(:,:,i),dF_total_nct(:,:,i),dF_nxt_nxt_nct(:,:,i), Winrnd(:,:,i),thetawin,Ret(:,i)',Betawh,Betaump,thetaS,PartyE_V(i,1));
% %         if (sum(squeeze(DC(5,1,:,i)-C(5,1,:,i)))<0)
% %             i
% %             C(:,:,:,i)
% %             DC(:,:,:,i)
% % 
% %         end
% %     for k=1:N
% %         for j=1:10
% %             if C(1,j,k,i)==1 %Contest
% %                 Continue1(i,k)=Continue1(i,k)+0.9*(ben1*C(2,j,k,i).^2-cost1*C(3,j,k,i).^2+C(5,j,k,i));
% %             else
% %                 Continue1(i,k)=Continue1(i,k)+0.9*(ben2*C(2,j,k,i).^2-cost2*C(3,j,k,i).^2+C(5,j,k,i));
% %             end
% %             if DC(1,j,k,i)==1 %Contest
% %                 DContinue1(i,k)=DContinue1(i,k)+0.9*(ben1*DC(2,j,k,i).^2-cost1*DC(3,j,k,i).^2+DC(5,j,k,i));
% %             else
% %                 DContinue1(i,k)=DContinue1(i,k)+0.9*(ben2*DC(2,j,k,i).^2-cost2*DC(3,j,k,i).^2+DC(5,j,k,i));
% %             end
% %         end
% %     end
% %     Continue(i,1)=Continuation(XSEV_(i,:)',XQEV(i,1),TenureE_V(i,1),LOGW_NXT_E_V(i,1),...
% %         Shockwh(:,:,i), Shockump(:,:,i), Entry(:,:,i), theta, E_VCTa, E_VCTt, gammaCT, E_VNCTa, E_VNCTt, gammaNCT, dF_gamma_ct(:,:,i), dF_total_ct(:,:,i),...
% %         dF_nxt_nxt_ct(:,:,i),dF_gamma_nct(:,:,i),dF_total_nct(:,:,i),dF_nxt_nxt_nct(:,:,i), Winrnd(:,:,i),thetawin,Ret(:,i)',Betawh,Betaump,cost1,ben1,cost2,ben2,thetaS,PartyE_V(i,1));
% %     DContinue(i,1)=Continuation(XSEV_(i,:)',XQEV(i,1),TenureE_V(i,1),LOGW_NXT_E_V(i,1)+Delt,...
% %         Shockwh(:,:,i), Shockump(:,:,i), Entry(:,:,i), theta, E_VCTa, E_VCTt, gammaCT, E_VNCTa, E_VNCTt, gammaNCT, dF_gamma_ct(:,:,i), dF_total_ct(:,:,i),...
% %         dF_nxt_nxt_ct(:,:,i),dF_gamma_nct(:,:,i),dF_total_nct(:,:,i),dF_n
% %         xt_nxt_nct(:,:,i), Winrnd(:,:,i),thetawin,Ret(:,i)',Betawh,Betaump,cost1,ben1,cost2,ben2,thetaS,PartyE_V(i,1));
% end
% 
%   save C.mat C
  %save DC.mat DC
%%
 load('./C5000-1-5000.mat');
 C=C(:,:,1:100,:);
load('residstd.mat');
% load('./DC.mat');

% LB(1)=-Inf; % normalized to 1, vacuous.
% LB(2)=-Inf; % ben1, coefficient on benefit
% LB(3)=-1; % alpha, benefit of spending
% LB(4)=1; % beta, cost or fundraising
% LB(5)=-1; % sig
% UB(1)=Inf;
% UB(2)=Inf;
% UB(3)=1;
% UB(4)=Inf;
% UB(5)=1;
% Aineq=[0,0,1,-1,0];
% Bineq=0;
% LB(1)=-Inf; % normalized to 1, vacuous.
% LB(2)=-Inf; % ben1, coefficient on benefit
% LB(3)=0; % alpha, benefit of spending
% LB(4)=0; % beta, cost or fundraising
% LB(5)=-Inf; % sig
% UB=[];

%Compute probability of winning from step 1-3
probwin=[ones(length(NCE_VCT),1),LOGW_E_VCT,XS_EVCT_(:,1).*SameE_VCT,XS_EVCT_(:,1).^2.*SameE_VCT,XS_EVCT_(:,2).*PartyE_VCT,LOGW_E_VCT.^2,X_KnotE_VCT,log(TenureE_VCT+1),log(TenureE_VCT+1).^2,...
    X_KnotE_VCT.*(XS_EVCT_(:,1).*SameE_VCT*ones(1,8)),X_KnotE_VCT.*(XS_EVCT_(:,2).*PartyE_VCT*ones(1,8)),X_KnotE_VCT.*(log(TenureE_VCT+1)*ones(1,8)),...;
  LOGW_E_VCT.^3]*Est3;


%Input for the function evaluated
datasetV=[NCE_V, LOGTot_NCE_V,XQEV2,LOGW_NXT_E_V, SameE_V,PartyE_V,TenureE_V, XSEV_,PresdumE_V,MidtermE_V];
datasetVCT=[E_VNContestFUL,XQEVCT2,LOGTot_NCE_VCT, NCE_VCT,LOGW_NXT_E_VCT,LOGTotal_E_VCT,VSEVCT,XS_EVCT_,TenureE_VCT, PartyE_VCT,SameE_VCT, LOGD_E_VCT,LOGD_E_VC,XQEVCT2];
datasetVNCT=[E_VContestFUL,LOGW_NXT_E_VNCT,LOGTotal_E_VNCT,NCE_VNCT,LOGTot_NCE_VNCT,XQEVNCT2];


%%
mintheta=[];
SRR=[];
 options=optimset('MaxIter',6000,'MaxFunEvals',1000000,'Display','iter');
 iter=11;  
 Sofarbest=10^8;
   minthetadist=zeros(5,iter);
   initthetadist=zeros(5,iter);
   srrdist=zeros(1,iter);
   tic
 parfor i=1:iter
 %load inittheta2step
 %theta0([1,2,3,5])=inittheta2step;
  %  theta0(4)=rand(1,1);
  theta0=rand(5,1);
  theta0(1)=theta0(1)*20;
  theta0(2)=theta0(2)*20;
  theta0(3)=theta0(3)*20;
  theta0(5)=theta0(5)*0.01;
  %theta0(6)=theta0(3)*10*rand(1,1);     
  
funvalue=Minimize2S_new(Est2,theta0,probwin,datasetV,datasetVCT,datasetVNCT,N,C,T);
funvalueend=0;
thetaup=theta0;
sss=1;
while abs(funvalue-funvalueend)>0.01
   sss=sss+1;
   funvalue=Minimize2S_new(Est2,thetaup,probwin,datasetV,datasetVCT,datasetVNCT,N,C,T);
    [mintheta,SRR]=fminsearch(@(theta2) Minimize2S_new(Est2,theta2,probwin,datasetV,datasetVCT,datasetVNCT,N,C,T,residstd),thetaup,options);
   funvalueend=Minimize2S_new(Est2,mintheta,probwin,datasetV,datasetVCT,datasetVNCT,N,C,T);
      if mod(sss,30)==1
      thetaup=mintheta+0.5*(rand(size(mintheta,1),1)-0.5).*mintheta;
      else
        thetaup=mintheta;
     end
end
 minthetadist(:,i)=mintheta;
 srrdist(i)=SRR;
 initthetadist(:,i)=theta0;
%theta2ndstep=aux;
%save minimizedtheta.txt mintheta SRR -ASCII
% 
% 
% Qualc(:,1)=E_V_july8(:,31);
% Qualc(:,2)=E_V_july8(:,32);
% Qualc(:,3)=E_V_july8(:,33);
% Qualc(:,4)=E_V_july8(:,34);
% Qualc(E_VContestFUL,:)=[];
% load q_C_E_VCT.mat
% Qualc(:,5)=q_C_E_VCT;
% save Qualc.txt Qualc -ASCII
% save QualI.txt XQEV2 -ASCII
 end
 toc
%%
%Check winning prob and distribution
thetain=Est2;


theta2=mintheta;
%LOGD_E_VC(LOGD_E_VC<9.21)=9.21;


vdelta=0.90;
% thetaS=thetain(1:2,1);
B_I=thetain(1,1);
B_C=thetain(2,1);
B_T=thetain(5,1);
thetaS2=thetain(3:4,1);
% for k=1:13
%     ACDa(:,k)=thetain(24+12*(k-1):24+12*k-1,1);
% end
% for k=1:13
%     ACDt(:,k)=thetain(180+12*(k-1):180+12*k-1,1);
% end
% for k=1:13
%     gammaACD(:,k)=thetain(336+12*(k-1):336+12*k-1,1);
% end


%h=1;
cost1=abs(theta2(1,1));  %% coefficient on the cost function of incumbent, contested
ben1=abs(theta2(2,1));   %% coefficient on the benefit function of incumbent, contested
ben2=ben1;
%ben2=abs(theta2(3,1));   %% coefficient on the benefit function of incumbent, uncontested
% cost_c=theta2(5,1); %% coefficient on the cost function of challenger
% ben_c=theta2(6,1);  %% coefficient on the benefit funciton of challenger.
% RES=theta2(7,1); %% reservation value

alpha=abs(theta2(4,1));
% alpha=cdf('norm',theta2(3,1),0,1);
% beta=1+abs(theta2(4,1));
beta=abs(theta2(3,1));
cost2=1;%abs(theta2(6,1));  %% coefficient on the cost function of the incumbent, uncontested >normalized to 1. Can be normalized
          % because we make C_I(x)=(x^alpha)*exp(f(q)) and we treat f(q)
          % non parametrically. If we set f(q)=C+f(q) then
          % C_I(x)=C*(x^alpha)*exp(f(q))
sig=abs(theta2(5,1));




%Continuation(ST,q_I,Ten,w_I,epswh, epsump, E_VCTa, E_VCTt, gamma,dF_gamma_ct, dF_total_ct, dF_nxt_nxt,Winrnd,thetawin,Ret,Betawh,Betaump,cost,ben,thetaS,Party)
%%%%%%%%          E_V         %%%%%%%%%
%Contcontest=payoff if contested
costpara=permute(repmat((ones(N,1)*((exp(LOGTot_NCE_V(:,1))./exp(NCE_V(:,1))).*...
                    ((max(0,NCE_V(:,1)).^(alpha-1))./(max(0,LOGTot_NCE_V(:,1)).^(beta-1)))).'),[1 1 10]),[3 1 2]);
Contcontest=zeros(3,T,N,length(NCE_V));
Contcontest(1,:,:,:)=ben1*(max(0,C(2,:,:,:))).^alpha.*C(1,:,:,:); %Benefit
Contcontest(2,:,:,:)=-cost1*(alpha/beta)*ben2*costpara.*squeeze(C(3,:,:,:)).^beta.*squeeze(C(1,:,:,:)); %Cost
Contcontest(3,:,:,:)=C(5,:,:,:).*C(1,:,:,:); %Return from winning

%Discount
for j=1:T
    Contcontest(:,j,:,:)=((vdelta)^(j-1))*Contcontest(:,j,:,:);
end
%Sum over benefit+cost+winning and over t=1~T
Contcontestpath=squeeze(sum(sum(Contcontest,1),2));

%Contuncontest=payoff if uncontested
Contuncontest=zeros(3,T,N,length(NCE_V));
Contuncontest(1,:,:,:)=ben2*(max(0,C(2,:,:,:))).^alpha.*(1-C(1,:,:,:)); %Benefit
Contuncontest(2,:,:,:)=-cost2*(alpha/beta)*ben2*costpara.*squeeze(C(3,:,:,:)).^beta.*squeeze(1-C(1,:,:,:)); %Cost
Contuncontest(3,:,:,:)=C(5,:,:,:).*(1-C(1,:,:,:)); %Return from winning

%Discount
for j=1:T
    Contuncontest(:,j,:,:)=((vdelta)^(j-1))*Contuncontest(:,j,:,:);
end

%Sum over benefit+cost+winning and over t=1~T
Contuncontestpath=squeeze(sum(sum(Contuncontest,1),2));

%Sum over contest and uncontest
Continue1=(Contcontestpath+Contuncontestpath).';

%This gives us sample size* number of simulation matrix of continuation
%payoff


% Given State, Tenure, warchest, compute incumbent continuation value %%
%Continue1=zeros(length(NCE_V),N);
%DContinue1=zeros(length(NCE_V),N);
% for i=1:length(NCE_V)  %% (all elements of E_V)
%     for k=1:N % simulation
%         for j=1:(min(find(C(5,:,k,i)==0))) % number of periods
%             if C(1,j,k,i)==1 %Contest
%                 Continue1(i,k)=Continue1(i,k)+((vdelta)^(j-1))*(ben1*max(0,C(2,j,k,i))^alpha-cost1*(alpha/beta)*ben2*(exp(LOGTot_NCE_V(i,1))/exp(NCE_V(i,1)))*...
%                     ((max(0,NCE_V(i,1))^(alpha-1))/(max(0,LOGTot_NCE_V(i,1))^(beta-1)))*C(3,j,k,i)^beta+C(5,j,k,i));
%             else
%                 Continue1(i,k)=Continue1(i,k)+((vdelta)^(j-1))*(ben2*max(0,C(2,j,k,i))^alpha-cost2*(alpha/beta)*ben2*(exp(LOGTot_NCE_V(i,1))/exp(NCE_V(i,1)))*...
%                     ((max(0,NCE_V(i,1))^(alpha-1))/(max(0,LOGTot_NCE_V(i,1))^(beta-1)))*C(3,j,k,i)^beta+C(5,j,k,i));
%             end
% %             if DC(1,j,k,i)==1 %Contest
% %                 DContinue1(i,k)=DContinue1(i,k)+((0.9)^(j-1))*(ben1*max(0,DC(2,j,k,i))^(1/2)-cost1*(1/4)*(ben2/cost2)*RTotDE_V(i,1)*DC(3,j,k,i)^2+DC(5,j,k,i));
% %             else
% %                 DContinue1(i,k)=DContinue1(i,k)+((0.9)^(j-1))*(ben2*max(0,DC(2,j,k,i))^(1/2)-ben2*(1/4)*RTotDE_V(i,1)*DC(3,j,k,i)^2+DC(5,j,k,i));
% %             end
%         end
%     end
% end
Continue=mean(Continue1,2);
Continuation1=Continue;
Continuation1(E_VContestFUL,:)=[];           %% Continuation1 is the Continuation value for periods in which incumebent is contested.
Continuation2=Continue;
Continuation2(E_VNContestFUL,:)=[];          %%  Continuation1 is the Continuation value for periods in which incumebent is uncontested.

%Regress Continue on State variables to find the derivative. %
%Regressand=[ones(length(NCE_V),1),LOGW_NXT_E_V,LOGW_NXT_E_V.^2,XQEV2,TenureE_V,XSEV_(:,1).*SameE_V,XSEV_(:,2).*PartyE_V,PresdumE_V,MidtermE_V];
Regressand=[ones(length(NCE_V),1),LOGW_NXT_E_V,LOGW_NXT_E_V.^2/10,LOGW_NXT_E_V.^3/100,LOGW_NXT_E_V.^4/1000,XQEV2,TenureE_V,XSEV_(:,1).*SameE_V,XSEV_(:,2).*PartyE_V,PresdumE_V,MidtermE_V];
%Outlier=(LOGW_NXT_E_V<quantile(LOGW_NXT_E_V,.05));
%Outlier=[];
%Continue_san_OL=Continue;
%Continue_san_OL(Outlier,:)=[];
%Regressand_san_OL=Regressand;
%Regressand_san_OL(Outlier,:)=[];
coef=Regressand\Continue;
%coef=(inv(Regressand_san_OL'*Regressand_san_OL))*Regressand_san_OL'*Continue_san_OL;

%Use estimated derivative in computing FOC.
%Deriv=[zeros(length(Continuation1),1),ones(length(Continuation1),1),2*LOGW_NXT_E_VCT,zeros(length(Continuation1),1),zeros(length(Continuation1),5)]*coef;
%DerivNCT=[zeros(length(Continuation2),1),ones(length(Continuation2),1),2*LOGW_NXT_E_VNCT,zeros(length(Continuation2),1),zeros(length(Continuation2),5)]*coef;
Deriv=[zeros(length(Continuation1),1),ones(length(Continuation1),1),2*LOGW_NXT_E_VCT/10,3*LOGW_NXT_E_VCT.^2/100,4*LOGW_NXT_E_VCT.^3/1000,zeros(length(Continuation1),1),zeros(length(Continuation1),5)]*coef;
DerivNCT=[zeros(length(Continuation2),1),ones(length(Continuation2),1),2*LOGW_NXT_E_VNCT/10,3*LOGW_NXT_E_VNCT.^2/100,4*LOGW_NXT_E_VNCT.^3/1000,zeros(length(Continuation2),1),zeros(length(Continuation2),5)]*coef;
%
Deriv=max(Deriv,0.00001);
DerivNCT=max(DerivNCT,0.00001);

OUT=((beta*cost1*(alpha/beta)*ben2*(exp(LOGTot_NCE_VCT(:,1))./exp(NCE_VCT(:,1))).*...
    ((max(0,NCE_VCT(:,1)).^(alpha-1))./(max(0,LOGTot_NCE_VCT(:,1)).^(beta-1))).*LOGTotal_E_VCT.^(beta-1)).*(1./exp(LOGTotal_E_VCT(:,1))))./((vdelta*Deriv).*(1./exp(LOGW_NXT_E_VCT)));
DE=exp(LOGTotal_E_VCT(:,1)).*((vdelta*Deriv).*(1./exp(LOGW_NXT_E_VCT)));
SRR13=mean(abs(OUT-probwin),1);
std13=std(OUT-probwin);
SRR14=(OUT>1);
SRR15=(OUT<min(probwin));

Pen=max(OUT-1,0)-min(OUT,0);
Pen2=Pen;
Pen2(IND6CT,:)=[];
DEL_Pen=[find(Pen2<quantile(Pen2,.025));find(Pen2>quantile(Pen2,.975))];      %Cannot invert some observations:
%DEL2=find(Pen2>0);
%Pen2(DEL_Pen,:)=[];
OUT_m=OUT;
OUT_m(DEL_Pen,:)=[];
%ako=1/2+(1/pi)*atan((OUT-0.5)*h);
ako=min(max(OUT,0.000001),0.999999);
BX1=norminv(ako);
% BX1=norminv(max(0.01,min(0.99,(beta*(alpha/beta)*ben2*...
%     ((max(0,NCE_VCT(:,1)).^(alpha-1))./(max(0,LOGTot_NCE_VCT(:,1)).^(beta-1))).*LOGTotal_E_VCT.^(beta-1))./(vdelta*Deriv))));  %%BX1=(1/sig)*(-0.5+B_I*d_I-B_C*d_C+...)=(K) in the paper.
BX1sig=BX1;
BX1sig(IND6CT,:)=[];

q_C_E_VCT=(-1)*(sig*BX1-B_I*LOGD_E_VCT-B_C*LOGD_E_VC-[XS_EVCT_(:,1).*SameE_VCT,XS_EVCT_(:,2).*PartyE_VCT]*thetaS2-B_T*log(TenureE_VCT+1)-XQEVCT2);


% XXX=Pen;


VV=VSEVCT-0.5-B_I*LOGD_E_VCT-B_C*LOGD_E_VC-[XS_EVCT_(:,1).*SameE_VCT,XS_EVCT_(:,2).*PartyE_VCT]*thetaS2-XQEVCT2+q_C_E_VCT-B_T*log(TenureE_VCT+1);
VV(IND6CT,:)=[];
VV(DEL_Pen,:)=[];
LOGD_E_VCT_m=LOGD_E_VCT;
LOGD_E_VCT_m(IND6CT,:)=[];
LOGD_E_VCT_m(DEL_Pen,:)=[];
LOGD_E_VC_m=LOGD_E_VC;
LOGD_E_VC_m(IND6CT,:)=[];
LOGD_E_VC_m(DEL_Pen,:)=[];
XS_m=[XS_EVCT_(:,1).*SameE_VCT,XS_EVCT_(:,2).*PartyE_VCT];
XS_m(IND6CT,:)=[];
XS_m(DEL_Pen,:)=[];
XQEVCT_m=XQEVCT2;
XQEVCT_m(IND6CT,:)=[];
XQEVCT_m(DEL_Pen,:)=[];
q_C_E_VCT_m=q_C_E_VCT;
q_C_E_VCT_m(IND6CT,:)=[];
q_C_E_VCT_m(DEL_Pen,:)=[];
TenureE_VCT_m=TenureE_VCT;
TenureE_VCT_m(IND6CT,:)=[];
TenureE_VCT_m(DEL_Pen,:)=[];
SRR11_1=mean(VV.^2);     %To estimate sigma.
std11_1=std(VV);
SRR11_2=mean(VV.*LOGD_E_VCT_m)^2;
std11_2=std(VV.*LOGD_E_VCT_m);
SRR11_3=mean(VV.*LOGD_E_VC_m)^2;
std11_3=std(VV.*LOGD_E_VC_m);
SRR11_4=mean(VV.*XS_m(:,1))^2;
std11_4=std(VV.*XS_m(:,1));
SRR11_5=mean(VV.*XS_m(:,2))^2;
std11_5=std(VV.*XS_m(:,2));
SRR11_6=mean(VV.*XQEVCT_m)^2;
std11_6=std(VV.*XQEVCT_m);
SRR11_7=mean(VV.*q_C_E_VCT_m)^2;
std11_7=std(VV.*q_C_E_VCT_m);
SRR11_8=mean(VV.*log(TenureE_VCT_m+1))^2;
std11_8=std(VV.*log(TenureE_VCT_m+1));
SRR11=SRR11_1/std11_1+SRR11_2/std11_2+SRR11_3/std11_3+SRR11_4/std11_4...
    +SRR11_5/std11_5+SRR11_6/std11_6+SRR11_7/std11_7+SRR11_8/std11_8;
%VV(DEL2,:)=[];
%SRR11=(1/290)*(sum(VV))^2;
%SRR11P=(1/290)*sum(VV.^2);

SRR9=-log(sig)-(1/length(VV))*sum((VV.^2)/(2*sig^2));
SRR9=-SRR9;
%SRR9P=(1/290)*sum(VV.^2-sig^2);
%SRR9P_1=mean(((VV-mean(VV)).^2-sig^2))^2;
%std9P_1=std(((VV-mean(VV)).^2-sig^2));
%SRR9P_2=mean(sqrt(VV.^2)-sig)^2;
%std9P_2=std(sqrt(VV.^2)-sig);
% foc of incumbent, contested periods %%
FOC11=(beta*cost1*(alpha/beta)*ben2*(exp(LOGTot_NCE_VCT(:,1))./exp(NCE_VCT(:,1))).*...
    ((max(0,NCE_VCT(:,1)).^(alpha-1))./(max(0,LOGTot_NCE_VCT(:,1)).^(beta-1))).*LOGTotal_E_VCT.^(beta-1)).*(1./exp(LOGTotal_E_VCT(:,1)))...      %% d/dI C_I(total)
    -(B_I/(sig*exp(LOGD_E_VCT(:,1))))*normpdf(BX1).*(1+vdelta*Continuation1)-alpha*ben1*(1./exp(LOGD_E_VCT(:,1))).*(max(0,(LOGD_E_VCT))).^(alpha-1); %% d/d_I(P_2)(B+delta*EV+(d/dI)H_I())

FOC11(IND6CT,:)=[];
FOC11(DEL_Pen,:)=[];
%FOC11(DEL2,:)=[];
SRR10P=mean(FOC11.^2,1);
std10P=std(FOC11);
% FOC12=cost1*(1/4)*(ben2/cost2)*RTotDE_VCT.*((LOGTotal_E_VCT+Delt).^2-LOGTotal_E_VCT.^2)...      %% d/dI C_I(total)
%     -(B_I/sig)*normpdf(BX1).*(1+0.9*Continuation1)-ben1*((LOGD_E_VCT+Delt).^2-(LOGD_E_VCT).^2); %% d/d_I(P_2)(B+delta*EV+(d/dI)H_I())
% FOC12(IND6CT,:)=[];
% FOC12delta=cost1*(1/4)*(ben2/cost2)*RTotDE_VCT.*((LOGTotal_E_VCT+Delt).^2-LOGTotal_E_VCT.^2)... 
%     -(B_I/sig)*normpdf(BX1b1).*(1+0.9*Continuation1b1)-(ben1+incr)*(max(0,(LOGD_E_VCT+Delt)).^(1/2)-max(0,LOGD_E_VCT).^(1/2));
% FOC12delta(IND6CT,:)=[];
% SRR11=((1/incr)*((1/290)*sum(FOC12delta.^2,1)-(1/290)*sum(FOC12.^2,1)))^2;
% SRR11P=(1/290)*sum(FOC12.^2,1);

% %% foc of incumbent, contested periods %%


% %% foc of incumbent, uncontested periods %%
FOC21=vdelta*DerivNCT.*(1./exp(LOGW_NXT_E_VNCT))-... 
(beta*cost2*(alpha/beta)*ben2*(exp(LOGTot_NCE_VNCT(:,1))./exp(NCE_VNCT(:,1))).*...
    ((max(0,NCE_VNCT(:,1)).^(alpha-1))./(max(0,LOGTot_NCE_VNCT(:,1)).^(beta-1))).*LOGTotal_E_VNCT.^(beta-1)).*(1./exp(LOGTotal_E_VNCT(:,1)));   %% delta*(d/dw_I)E_V-C_I'
FOC21(IND6NCT,:)=[];

% FOC21delta=0.9*(DContinuation2b2-Continuation2b2)-(ben2+incr)*RTotDE_VNCT.*((LOGTotal_E_VNCT+Delt).^2-LOGTotal_E_VNCT.^2);   %% delta*(d/dw_I)E_V-C_I'
% FOC21delta(IND6NCT,:)=[];
% SRR12=((1/incr)*((1/290)*sum(FOC21delta.^2,1)-(1/290)*sum(FOC21.^2,1)))^2;
SRR12P=mean(FOC21.^2,1);
std12P=std(FOC21);

% SRR9
% SRR10P
% % SRR11P
% SRR12P

SRR2step=SRR9+10^7*SRR10P/std10P+SRR11+10^4*SRR12P/std12P+sum(SRR15)+sum(SRR14)+0.0000*SRR13/std13;%+;%;%
%%
figure(1)
hist(OUT,20)
figure(2)
hist(OUT(OUT<2),20)
%%
hist(DE(OUT<2),20)

%%
size(OUT(OUT>1))
%%
hist(XQEV2,20)
%%
hist(q_C_E_VCT,20)
%%
hist(q_C_E_VCT(Deriv<0.001),20)
%%
hist(q_C_E_VCT(LOGW_NXT_E_VCT<10),20)
%%

bins=linspace(-0.5,0.5,100);
[hist1,scale1]=hist(q_C_E_VCT,bins);
[hist2,scale2]=hist(XQEV2,bins);
bar(hist2)
hold on;
bar(hist1,'r')

%%
figure(1)
hist(LOGW_NXT_E_VCT,20)
figure(2)
hist(LOGW_NXT_E_VCT(q_C_E_VCT>0.1),20)
%%
figure(1)
hist(LOGW_E_VCT,20)
figure(2)
hist(LOGW_E_VCT(q_C_E_VCT>0.1),20)
%%
figure(1)
hist(LOGD_E_VC,20)
figure(2)
hist(LOGD_E_VC(q_C_E_VCT>0.1),20)
%%
figure(1)
hist(LOGD_E_VCT,20)
figure(2)
hist(LOGD_E_VCT(q_C_E_VCT>0.1),20)