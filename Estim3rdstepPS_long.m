
clear
% **
% make sure you get the sign of q_e right. Also note that f_d_qe, f_w_qe have the PARTY of Challenger.
% **

global LOGW_I
global Contest
global Party
global Entry
global Winrnd
global Tenure
global E_VContestFUL
global E_VNContestFUL
global NCE_V
global NCE_VCT
global NCE_VNCT
global PartyE_V
global PartyE_VCT
global XS_EVCT_
global XS_EVNCT_
global XS_EVCTwnxt_
global XSEV_
global TenureE_V
global TenureE_VCT
global TenureE_VNCT
global Win_E_VCT
global LOGW_NXT_E_VCT
global LOGW_NXT_E_VNCT
global LOGTotal_E_VCT
global LOGTotal_E_VNCT
global LOGTot_NCE_V
global LOGW_E_VCT
global LOGW_E_VNCT
global LOGW_E_VCTwnxt
global LOGD_E_VCT
global LOGD_E_VNCT
global N
global Sofarbest
global bestiter
global TenureE_VCTwnxt
global LOGW_NXT_E_V
global LOGW_NXT_E_VCTwnxt
global LOGW_NXT_E_VC
global LOGD_E_VC
global LOGTot_E_VC
global NumSim;
global IND6CTto7
global IND6CT
global IND6NCT
global RTotDE_V
global RTotDE_VCT
global XQEV2
global RTotDE_VNCT
global VSEVCT;
global PrE3step
global q_c3step
global LOGD_E_VC3step
global LOGTot_E_VC3step
global rtotd3step
global max_N
global LOGW_E_VCT3step
global LOGW_NXT3step
global XS3step_
global Tenure3step
global Party3step
global TenureE_VCT3step
global interest


% load ('./coef.mat')
% load ('./BX1.mat') 





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
load('./ q_C_E_VCT.mat');
%load ('./presseq.mat')
load ('./DEL_Pen.mat')
load ('./Continue.mat')

%Others
load ('./retire.txt');
load ('./stateevol.mat')



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


%%%%%%%%%%%%%%%%%%%%%%%% 
% Defining variables: done
%%%%%%%%%%%%%%%%%%%%%%%%
%%

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

%     
% 
% XS_EVCT2=XS_EVCT_*thetaS2;
% XS_EVCT=XS_EVCT_*thetaS;
% % XSEV=XSEV_*thetaS;
% XQ2=X_Knot1*thetaQ2;
% XQEV2=X_KnotEV1*thetaQ2;  % q_I
% XQEVCT2=XQEV2;
% XQEVCT2(E_VContestFUL,:)=[];
% XQEV=X_KnotEV1*thetaQ;  % q_I2
% XQEVCT=XQEV;
% XQEVCT(E_VContestFUL,:)=[];

%Results from first stage up here


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Sofarbest=10^8;
bestiter=0;
%Delt=0.5;
interest=0.1;
dfSim=100;              %dfSim is the number of simulations used to evaluate the conditional expectation of the entrant's quality.
N=20;                      % N is the number of simulations.
NumSim=50;              %NumSim is the number of simulations in the 1st step to obtain f(.,q_e)
max_N=max(OP_INC_july8(:,36));               % max_N is the maximum number of possible entrants in the Primary.
T=10;                      % T is the number of periods that we move the simlation forward.



%If create Ccmatrix, we need these shocks

% Entry=rand(T,N,length(E_V_july8));             %Simulation draw for computing the continuation value E_V
% Winrnd=rand(T,N,length(E_V_july8));
% 
% Retirernd=(rand(T-1,N,length(E_V_july8))<repmat(retire',[1 N, length(E_V_july8)]));
% Retirernd(T,:,:)=1;
% 
% 
% for i=1:N
%     for j=1:length(E_V_july8)
%         Ret(i,j)=find(Retirernd(:,i,j)==1,1,'first');
%     end
% end
% 
% 
% 
% dF_gamma_ct=rand(T,N,length(E_V_july8));
% dF_total_ct=rand(T,N,length(E_V_july8));
% dF_nxt_nxt_ct=rand(T,N,length(E_V_july8));
% 
% 
% dF_gamma_ctC=rand(T,N,length(E_V_july8));
% dF_total_ctC=rand(T,N,length(E_V_july8));
% dF_nxt_nxt_ctC=rand(T,N,length(E_V_july8));
% 
% 
% dF_gamma_ctR=rand(T,N,length(E_V_july8));
% dF_total_ctR=rand(T,N,length(E_V_july8));
% dF_nxt_nxt_ctR=rand(T,N,length(E_V_july8));
% 
% 
% %Shock of state variables
% Shockump=sqrt(epsump)*randn(T,N,length(E_V_july8));
% Shockpartisan=sqrt(epspartisan)*randn(T,N,length(E_V_july8));
% 
% Shockump2=sqrt(epsump)*randn(T,N,length(E_V_july8));
% Shockpartisan2=sqrt(epspartisan)*randn(T,N,length(E_V_july8));
% 
% Shockump3=sqrt(epsump)*randn(T,N,length(E_V_july8));
% Shockpartisan3=sqrt(epspartisan)*randn(T,N,length(E_V_july8));






q_c3step=q_C_E_VCT;
%q_c3step(IND6CT,:)=[];    %% Challenger quality.
q_c3step(DEL_Pen,:)=[]; % DEL_Pen outlier. cannot use.
rtotd3step=zeros(length(q_c3step),1);
for i=1:length(q_c3step)
    [minq,argminq]=min(abs(q_c3step(i,1)-XQEV2));
%     [result,r]=fmincon(@(rtotd) f_qinverse(rtotd),0.03,[],[],[],[],[],[],@(x) f_constraint(x,q_c3step(i,1),thetaQ2));   %%invert f_q function to get d/tot.
    rtotd3step(i,1)=RTotDE_V(argminq,1);
%    rtotd3step(i,1)=result;
end
% x_knot3step=(rtotd3step<mesh(2,1)).*(1-(rtotd3step-mesh(1,1))/(mesh(2,1)-mesh(1,1)));
% for i=0:6
%     PLUS=(rtotd3step>=mesh(i+1,1)).*(rtotd3step<mesh(i+2,1)).*((rtotd3step-mesh(i+1,1))/(mesh(i+2,1)-mesh(i+1,1)))...
%         +(rtotd3step>=mesh(i+2,1)).*(rtotd3step<mesh(i+3,1)).*(1-(rtotd3step-mesh(i+2,1))/(mesh(i+3,1)-mesh(i+2,1)));
%     x_knot3step=[x_knot3step,PLUS];
% end
% PLUS=(rtotd3step>=mesh(8,1)).*(rtotd3step-mesh(8,1))/(mesh(9,1)-mesh(8,1));
% x_knot3step=[x_knot3step,PLUS];
% q_c3step0=x_knot3step*thetaQ;
% save rtotd3step.mat rtotd3step
% load('./rtotd3step.mat');
[thetae,SRR]=fminsearch(@(x) entry_probit(x,Contest, LOGW_I,XS_, XQ2, Party,Tenure ),[theta(1:3,1);thetaS*theta(4,1);theta(5,1)]);
save thetae.mat thetae
XS3step_=XS_EVCT_;
%XS3step_(IND6CT,:)=[];
XS3step_(DEL_Pen,:)=[];
Party3step=(-1)*PartyE_VCT;  %% Challenger party=(-1)*Incumbent party
%Party3step(IND6CT,:)=[];
Party3step(DEL_Pen,:)=[];
Tenure3step=zeros(size(Party3step,1),1);               %% Challenger starts with Tenure==0;
LOGW_NXT3step=LOGW_NXT_E_VC;
%LOGW_NXT3step(IND6CT,:)=[];                                    %% challenger warchest.
LOGW_NXT3step(DEL_Pen,:)=[];
LOGW_E_VCT3step=LOGW_E_VCT;
%LOGW_E_VCT3step(IND6CT,:)=[];
LOGW_E_VCT3step(DEL_Pen,:)=[];
TenureE_VCT3step=TenureE_VCT;
%TenureE_VCT3step(IND6CT,:)=[];
TenureE_VCT3step(DEL_Pen,:)=[];

LOGD_E_VC3step=LOGD_E_VC;
%LOGD_E_VC3step(IND6CT,:)=[];            %challenger spending
LOGD_E_VC3step(DEL_Pen,:)=[];
LOGTot_E_VC3step=LOGTot_E_VC;
%LOGTot_E_VC3step(IND6CT,:)=[];              %challenger tot
LOGTot_E_VC3step(DEL_Pen,:)=[];
VSEVCT3step=VSEVCT;
%VSEVCT3step(IND6CT,:)=[];
VSEVCT3step(DEL_Pen,:)=[];
IND6CTto7=find(LOGW_NXT3step==0);        %IND6CTto7 indexes those for which war chest next period==0 for challenger.
% % 
% Cc=zeros(5,T,N,length(Party3step));
% % DCc=zeros(5,T,N,length(Party3step));
% for i=1:length(Party3step)
%     Cc(:,:,:,i)=Actions3(i,XS3step_(i,:)',thetae,q_c3step(i,1),Tenure3step(i,1),LOGW_NXT3step(i,1),...
%         Shockwh(:,:,i), Shockump(:,:,i), Entry(:,:,i), E_VCTa, E_VCTt, gammaCT, gammaNCT, tNCT, wNCT, dF_gamma_ct(:,:,i), dF_total_ct(:,:,i),...
%         dF_nxt_nxt_ct(:,:,i), Winrnd(:,:,i),thetawin,Ret(:,i)',Betawh,Betaump,thetaS,thetaS2,Party3step(i,1));
% %     DCc(:,:,:,i)=Actions(i,XS3step_(i,:)',q_c3step(i,1),rtotd3step(i,1),Tenure3step(i,1),LOGW_NXT3step(i,1)+Delt,...
% %         Shockwh(:,:,i), Shockump(:,:,i), Entry(:,:,i), theta, E_VCTa, E_VCTt, gammaCT, E_VNCTa, E_VNCTt, gammaNCT, dF_gamma_ct(:,:,i), dF_total_ct(:,:,i),...
% %         dF_nxt_nxt_ct(:,:,i),dF_gamma_nct(:,:,i),dF_total_nct(:,:,i),dF_nxt_nxt_nct(:,:,i), Winrnd(:,:,i),thetawin,Ret(:,i)',Betawh,Betaump,thetaS,Party3step(i,1));
% end
% save Cc.mat Cc
% % save DCc.mat DCc

load('./Cc.mat');
% load('./DCc.mat');
load('./C.mat');
% % load('./DC.mat');
% %% compute PrE %%
% %% probability of entry ex-ante %%

%%

%Inpute continuation payoff
Regressand=[ones(length(NCE_V),1),LOGW_NXT_E_V,LOGW_NXT_E_V.^2/10,LOGW_NXT_E_V.^3/100,LOGW_NXT_E_V.^4/1000,...
    XQEV2,LOGW_NXT_E_V.*XQEV2,LOGW_NXT_E_V.^2/10.*XQEV2,LOGW_NXT_E_V.^3/100.*XQEV2,LOGW_NXT_E_V.^4/1000.*XQEV2,...
    TenureE_V,TenureE_V.^2,XSEV_(:,1).*SameE_V,XSEV_(:,1).^2.*SameE_V,XSEV_(:,2).*PartyE_V,XSEV_(:,2).^2.*PartyE_V,...
    TenureE_V.*XSEV_(:,1).*SameE_V,TenureE_V.*XSEV_(:,2).*PartyE_V,PresdumE_V,MidtermE_V];
coef=Regressand\Continue;



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
% %BX1=abs(1/Sndstage_kettei_long(5,1))*(B_I*LOGD_E_VCT+B_C*LOGD_E_VC+XS_EVCT2.*PartyE_VCT+B_T*log(TenureE_VCT+1)+XQEVCT2+q_C_E_VCT);
% BX13step=BX1;
% BX13step(IND6CT,:)=[];
% BX13step(DEL_Pen,:)=[];
% thetain2=Sndstage_kettei_long;
% % Estim3=[   -1.9713322e+01
% %   -1.1704758e-02
% %    1.6123112e+03
% %   -7.4047131e+033]

for i=1:5
thetain3=Trdstage_initial(1:4,i);
[mintheta,SRR]=fminsearch(@(theta3) Minimize3PS_long(thetain1step,thetain2,theta3),thetain3);
end

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
