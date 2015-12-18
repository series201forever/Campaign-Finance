clear all
% **
% make sure you get the sign of q_e right. Also note that f_d_qe, f_w_qe have the PARTY of Challenger.
% **

rng(958)

global iterate
global results
global LOGW_I
global Contest
global Party
global dF_gamma_ct
global dF_total_ct
global dF_nxt_nxt_ct
global dF_gamma_nct
global dF_total_nct
global dF_nxt_nxt_nct
global dF_gamma_ctC
global dF_total_ctC
global dF_nxt_nxt_ctC
global dF_gamma_nctC
global dF_total_nctC
global dF_nxt_nxt_nctC
global dF_gamma_ctR
global dF_total_ctR
global dF_nxt_nxt_ctR
global dF_gamma_nctR
global dF_total_nctR
global dF_nxt_nxt_nctR
global dF_dgammaC
global dF_wgammaC
global dF_vtgammaC
global Entry
global Winrnd
global Ret
global RetC
global RetR
global Betawh
global Betaump
global epswh
global epsump
global Shockump
global Shockwh
global epswh2
global epsump2
global Shockump2
global Shockwh2
global epswh3
global epsump3
global Shockump3
global Shockwh3
global Samplesize
global Contest
global LOGW_I
global LOGD_NC
global LOGW_NC
global LOGWNXT_NC
global Tenure
global Party
global White
global White_NC
global Unemployment
global Unemployment_NC
global SamplesizeAC
global LOGD_IAC
global LOGW_IAC
global VSAC
global NCAC
global LOGD_CAC
global PartyAC
global TenureAC
global WhiteAC
global White_NCAC
global UnemploymentAC
global Unemployment_NCAC
global E_VContestFUL
global E_VNContestFUL
global NCE_V
global NCE_VCT
global NCE_VNCT
global PartyE_V
global PartyE_VCT
global PartyE_VNCT
global XSNCEVCT_
global XSNCEVNCT_
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
global LOGW_E_V
global LOGW_E_VCT
global LOGW_E_VNCT
global LOGW_E_VCTwnxt
global LOGD_E_VCT
global LOGD_E_VNCT

global Sofarbest
global bestiter
global NCE_VCTwnxt
global TenureE_VCTwnxt
global PartyE_VCTwnxt
global XSNCEVCTwnxt_
global LOGW_NXT_E_V
global LOGW_NXT_E_VCTwnxt
global LOGW_NXT_E_VC
global LOGD_E_VC
global LOGTot_E_VC
global LOGTot_NCE_VCT
global YearAC
global YEARE_VCT
global YEARE_VNCT
global YEARE_VCTwnxt
global LOGW_CNXTAC
global NumSim;
global Year
global IND6CTto7
global IND6CT
global IND7CT
global YEARE_VCT
global YEARE_VNCT
global YEARE_VCTwnxt
global IND6NCT
global C
global LOGTot_NCE_VNCT
% global DC
% global RTotDE_V
% global RTotDE_VCT
global XQEV2
global XQEVCT2
global XQEVNCT2
global LOGTot_NCE_V
global coef
% global RTotDE_VNCT
global Tenure_INCCT
global VSEVCT;
global XXX
global incr


%Global used to construct C.mat
global N
global interest

%load ('./E_V_july8.mat')
E_V_july8=csvread('E_Vjuly8partisan.csv',1,0);

% load('./op_inc_july8.mat')
OP_INC_july8=csvread('OP_INC_july8partisan.csv',1,0);

% load('./op_inc_iv_july8.mat')
OP_INC_IV_july8=csvread('OP_INC_IV_july8partisan.csv',1,0);

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
% load ('./Sndstage_initial_long.txt')

% rand('state',1000);
% randn('state',10);

dele1=find(sum(isnan(OP_INC_IV_july8),2)>0);
OP_INC_IV_july8(dele1,:)=[]; %Drop NaN
OP_INC_IV_july8(OP_INC_IV_july8(:,2)>2002,:)=[]; %Drop year 2004 and on
deleA=find(OP_INC_IV_july8(:,16)==0);
OP_INC_IV_july8(deleA,:)=[];%Drop if oppornent disburse=0
OP_INC_IV_july8(102,:)=[];%Drop an Rtotd outlier
deleA2=find(OP_INC_IV_july8(:,12)<1);
deleA3=find(OP_INC_IV_july8(:,11)<1);
OP_INC_IV_july8(deleA2,:)=[];%Drop if nonpositive saving
OP_INC_IV_july8(deleA3,:)=[];%Drop if nonpositive begcash

dele2=find(sum(isnan(OP_INC_july8),2)>0);
OP_INC_july8(dele2,:)=[]; %Drop NaN
OP_INC_july8(OP_INC_july8(:,2)>2002,:)=[]; %Drop year 2004 and on
deleB=find(((OP_INC_july8(:,16)==0).*OP_INC_july8(:,8))==1);
OP_INC_july8(deleB,:)=[];%Drop if oppornent disburse=0 and contested
OP_INC_july8(119:121,:)=[];%Drop an Rtotd outlier
deleB2=find(OP_INC_july8(:,12)<1);
deleB3=find(OP_INC_july8(:,11)<1);
OP_INC_july8(deleB2,:)=[];%Drop if nonpositive saving
OP_INC_july8(deleB3,:)=[];%Drop if nonpositive begcash

dele3=find(sum(isnan(E_V_july8),2)>0); 
E_V_july8(dele3,:)=[]; %Drop NaN
E_V_july8(E_V_july8(:,2)>2002,:)=[]; %Drop year 2004 and on
deleC=find((E_V_july8(:,16)==0).*E_V_july8(:,8)==1);
E_V_july8(deleC,:)=[]; %Drop if oppornent disburse=0 and contested
E_V_july8(174:177,:)=[]; %Drop an Rtotd outlier
deleC2=find(E_V_july8(:,12)<1);
deleC3=find(E_V_july8(:,11)<1);
E_V_july8(deleC2,:)=[];%Drop if nonpositive saving
E_V_july8(deleC3,:)=[];%Drop if nonpositive begcash

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Defining Variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
 Estimation_OP=OP_INC_july8; %Sample to estimate entry probability
% 
% 
 Samplesize=length(Estimation_OP);
% %Year=Estimation_OP(:,2); % year of election.
% Contest=Estimation_OP(:,8); % contest
% Primary_N=Estimation_OP(:,36);
% 
% %LOGD_I=log(max(ones(Samplesize,1),Estimation_OP(:,10))); % spend
% LOGW_I=log(max(ones(Samplesize,1),Estimation_OP(:,11))); % w_I
% %LOGW_INXT=log(max(ones(Samplesize,1),Estimation_OP(:,12))); % W_I'
% %VS=Estimation_OP(:,13); %V_I
% %IDN=Estimation_OP(:,21); %Candidate ID
 LOGTot_NC=log(max(ones(Samplesize,1),Estimation_OP(:,31))); % Tot_INC
 LOGD_NC=log(max(ones(Samplesize,1),Estimation_OP(:,32))); % D_INC
% %LOGW_NC=log(max(ones(Samplesize,1),Estimation_OP(:,33))); % W_INC
% %LOGWNXT_NC=log(max(ones(Samplesize,1),Estimation_OP(:,34))); %W'_INC
 RTotD_NC=(max(0,LOGD_NC).^(-1/2))./LOGTot_NC;  %% ratio of spending over total. (modify)
% RTotD_NC=RTotD_NC.*(exp(LOGTot_NC)./exp(LOGD_NC));
% % LOGLOGD_NC=log(max(ones(Samplesize,1),LOGD_NC));
% % LOGLOGTot_NC=log(max(ones(Samplesize,1),LOGTot_NC));
% %Tenure_NC=Estimation_OP(:,35); % Tenure_INC
% Tenure=Estimation_OP(:,15);
% Party=3*ones(Samplesize,1)-2*Estimation_OP(:,3);  %%Party=1  if candidate i is Democrat and Party=-1 if candidate i is Republican.
% Presparty=3*ones(Samplesize,1)-2*Estimation_OP(:,39);
% Same=Party==Presparty;
% Dif=Party~=Presparty;
% Same=Same-Dif;%Same party=1 if candidate party=President party, otherwise -1
% % White=Estimation_OP(:,21);
% % Black=Estimation_OP(:,22);
% % Other=Estimation_OP(:,23);
% % White_NC=Estimation_OP(:,27);
% % Black_NC=Estimation_OP(:,28);
% % Other_NC=Estimation_OP(:,29);
% Unemployment=Estimation_OP(:,24);
% % Unemployment_NC=Estimation_OP(:,30);
% Partisan=Estimation_OP(:,63);
% 
% Unempsame=Unemployment.*Same;
% Partdemo=Partisan.*Party;
% 
% %%B-Spline for q_I(RTotD_NC), where RTotD_NC=(LOGD_NC./LOGTot_NC)  etc.%%
% %% take knots to be between 0.88 to 1.025 with 8 knots (8 basis functions). (Almost all lie
% %% within this range) Let X_Knot be the matrix with 8 columns that contain
% %% the value of the B-Spline basis function evaluated at each of the 8
% %% basis functions.
% % mesh=quantile(LOGLOGD_NC,[.125;.25;.375;.5;.625;.75;.875]);
% % mesh=[min(LOGLOGD_NC);mesh;max(LOGLOGD_NC)];
% % X_Knot1=(LOGLOGD_NC<mesh(2,1)).*(1-(LOGLOGD_NC-mesh(1,1))/(mesh(2,1)-mesh(1,1)));
% % for i=1:6
% %     PLUS=(LOGLOGD_NC>=mesh(i+1,1)).*(LOGLOGD_NC<mesh(i+2,1)).*((LOGLOGD_NC-mesh(i+1,1))/(mesh(i+2,1)-mesh(i+1,1)))...
% %         +(LOGLOGD_NC>mesh(i+2,1)).*(LOGLOGD_NC<mesh(i+3,1)).*(1-(LOGLOGD_NC-mesh(i+2,1))/(mesh(i+3,1)-mesh(i+2,1)));
% %     X_Knot1=[X_Knot1,PLUS];
% % end
% % PLUS=(LOGLOGD_NC>=mesh(8,1)).*(LOGLOGD_NC-mesh(8,1))/(mesh(9,1)-mesh(8,1));
% % X_Knot1=[X_Knot1,PLUS];
% % 
% % mesh=quantile(LOGLOGTot_NC,[.125;.25;.375;.5;.625;.75;.875]);
% % mesh=[min(LOGLOGTot_NC);mesh;max(LOGLOGTot_NC)];
% % X_Knot2=(LOGLOGTot_NC<mesh(2,1)).*(1-(LOGLOGTot_NC-mesh(1,1))/(mesh(2,1)-mesh(1,1)));
% % for i=1:6
% %     PLUS=(LOGLOGTot_NC>=mesh(i+1,1)).*(LOGLOGTot_NC<mesh(i+2,1)).*((LOGLOGTot_NC-mesh(i+1,1))/(mesh(i+2,1)-mesh(i+1,1)))...
% %         +(LOGLOGTot_NC>mesh(i+2,1)).*(LOGLOGTot_NC<mesh(i+3,1)).*(1-(LOGLOGTot_NC-mesh(i+2,1))/(mesh(i+3,1)-mesh(i+2,1)));
% %     X_Knot2=[X_Knot2,PLUS];
% % end
% % PLUS=(LOGLOGTot_NC>=mesh(8,1)).*(LOGLOGTot_NC-mesh(8,1))/(mesh(9,1)-mesh(8,1));
% % X_Knot2=[X_Knot2,PLUS];
% 
% fineness=9;
% mesh=quantile(RTotD_NC,linspace(0,1,fineness).');
% X_Knot1=(RTotD_NC<mesh(2,1)).*(1-(RTotD_NC-mesh(1,1))/(mesh(2,1)-mesh(1,1)));
% for i=0:(numel(mesh)-3)
%     PLUS=(RTotD_NC>=mesh(i+1,1)).*(RTotD_NC<mesh(i+2,1)).*((RTotD_NC-mesh(i+1,1))/(mesh(i+2,1)-mesh(i+1,1)))...
%         +(RTotD_NC>=mesh(i+2,1)).*(RTotD_NC<mesh(i+3,1)).*(1-(RTotD_NC-mesh(i+2,1))/(mesh(i+3,1)-mesh(i+2,1)));
%     X_Knot1=[X_Knot1,PLUS];
% end
% PLUS=(RTotD_NC>=mesh((numel(mesh)-1),1)).*(RTotD_NC-mesh((numel(mesh)-1),1))/(mesh(numel(mesh),1)-mesh((numel(mesh)-1),1));
% X_Knot1=[X_Knot1,PLUS];
% 
% %First column:redundant
% X_Knot1=X_Knot1(:,2:fineness);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%      Ai & Chen     %%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% AC=OP_INC_IV_july8;%sample to estimate vote share
% 
% %Variables used in Number=2, estimation of vote share
% SamplesizeAC=length(OP_INC_IV_july8);
% LOGD_IAC=log(max(ones(SamplesizeAC,1),AC(:,10))); % spend
% LOGW_IAC=log(max(ones(SamplesizeAC,1),AC(:,11))); % w_I
% %LOGW_INXTAC=log(max(ones(SamplesizeAC,1),AC(:,12))); % W_I'
% VSAC=AC(:,13); %V_I
% LOGD_NCAC=log(max(ones(SamplesizeAC,1),AC(:,32))); % D_INC
% %LOGW_NCAC=log(max(ones(SamplesizeAC,1),AC(:,33))); % W_INC
% LOGTot_NCAC=log(max(ones(SamplesizeAC,1),AC(:,31))); % Tot_INC
% % LOGLOGD_NCAC=log(max(ones(SamplesizeAC,1),LOGD_NCAC));
% % LOGLOGTot_NCAC=log(max(ones(SamplesizeAC,1),LOGTot_NCAC));
% RTotD_NCAC=(max(0,LOGD_NCAC).^(-1/2))./LOGTot_NCAC; %% D_INC./Tot_INC
% RTotD_NCAC=RTotD_NCAC.*(exp(LOGTot_NCAC)./exp(LOGD_NCAC));
% %LOGWNXT_NCAC=log(max(ones(SamplesizeAC,1),AC(:,34))); %W'_INC
% %Tenure_NCAC=AC(:,35); % Tenure_INC
% %NCAC=[LOGD_NCAC,LOGW_NCAC,LOGWNXT_NCAC];
% LOGD_CAC=log(max(ones(SamplesizeAC,1),AC(:,16))); %Sepnding of Challenger
% % LOGW_CAC=log(max(ones(SamplesizeAC,1),AC(:,17))); %Warchest of Challenger
% % LOGW_CNXTAC=log(max(ones(SamplesizeAC,1),AC(:,18))); %Savings of Challenger
% % YearAC=AC(:,2); % Current year.
% 
% TenureAC=AC(:,15);
% PartyAC=3*ones(SamplesizeAC,1)-2*AC(:,3);  %%PartyAC=1  if candidate i is Democrat and PartyAC=-1 if candidate i is Republican.
% PrespartyAC=3*ones(SamplesizeAC,1)-2*AC(:,39);
% SameAC=PartyAC==PrespartyAC; %Same party if the two match
% DifAC=PartyAC~=PrespartyAC;
% SameAC=SameAC-DifAC;
% 
% % WhiteAC=AC(:,21);
% % BlackAC=AC(:,22);
% % OtherAC=AC(:,23);
% % White_NCAC=AC(:,27);
% % Black_NCAC=AC(:,28);
% % Other_NCAC=AC(:,29);
% UnemploymentAC=AC(:,24);
% UnemploymentsqAC=UnemploymentAC.^2;
% %Unemployment_NCAC=AC(:,30);
% PartisanAC=AC(:,63);
% 
% 
% UnempsameAC=UnemploymentAC.*SameAC;
% UnempsqsameAC=UnemploymentsqAC.*SameAC;
% PartdemoAC=PartisanAC.*PartyAC;
% 
% 
% %Variables used in number=3, estimation of winning probability
% Win_AC=(AC(:,13)>=0.5);  % Dummy for whether incumbent won the election in (t).
% 
% 
% 
% %%B-Spline for q_I(RTotD_NCAC), where RTotD_NCAC=(LOGD_NCAC./LOGTot_NCAC)  etc.%%
% %% take knots to be between 0.88 to 1.025 with 8 knots (8 basis functions). (Almost all lie
% %% within this range) Let X_Knot be the matrix with 8 columns that contain
% %% the value of the B-Spline basis function evaluated at each of the 8
% %% basis functions.
% % 
% % mesh=quantile(LOGLOGD_NCAC,[.125;.25;.375;.5;.625;.75;.875]);
% % mesh=[min(LOGLOGD_NCAC);mesh;max(LOGLOGD_NCAC)];
% % X_KnotAC1=(LOGLOGD_NCAC<mesh(2,1)).*(1-(LOGLOGD_NCAC-mesh(1,1))/(mesh(2,1)-mesh(1,1)));
% % for i=1:6
% %     PLUS=(LOGLOGD_NCAC>=mesh(i+1,1)).*(LOGLOGD_NCAC<mesh(i+2,1)).*((LOGLOGD_NCAC-mesh(i+1,1))/(mesh(i+2,1)-mesh(i+1,1)))...
% %         +(LOGLOGD_NCAC>mesh(i+2,1)).*(LOGLOGD_NCAC<mesh(i+3,1)).*(1-(LOGLOGD_NCAC-mesh(i+2,1))/(mesh(i+3,1)-mesh(i+2,1)));
% %     X_KnotAC1=[X_KnotAC1,PLUS];
% % end
% % PLUS=(LOGLOGD_NCAC>=mesh(8,1)).*(LOGLOGD_NCAC-mesh(8,1))/(mesh(9,1)-mesh(8,1));
% % X_KnotAC1=[X_KnotAC1,PLUS];
% % 
% % mesh=quantile(LOGLOGTot_NCAC,[.125;.25;.375;.5;.625;.75;.875]);
% % mesh=[min(LOGLOGTot_NCAC);mesh;max(LOGLOGTot_NCAC)];
% % X_KnotAC2=(LOGLOGTot_NCAC<mesh(2,1)).*(1-(LOGLOGTot_NCAC-mesh(1,1))/(mesh(2,1)-mesh(1,1)));
% % for i=1:6
% %     PLUS=(LOGLOGTot_NCAC>=mesh(i+1,1)).*(LOGLOGTot_NCAC<mesh(i+2,1)).*((LOGLOGTot_NCAC-mesh(i+1,1))/(mesh(i+2,1)-mesh(i+1,1)))...
% %         +(LOGLOGTot_NCAC>mesh(i+2,1)).*(LOGLOGTot_NCAC<mesh(i+3,1)).*(1-(LOGLOGTot_NCAC-mesh(i+2,1))/(mesh(i+3,1)-mesh(i+2,1)));
% %     X_KnotAC2=[X_KnotAC2,PLUS];
% % end
% % PLUS=(LOGLOGTot_NCAC>=mesh(8,1)).*(LOGLOGTot_NCAC-mesh(8,1))/(mesh(9,1)-mesh(8,1));
% % X_KnotAC2=[X_KnotAC2,PLUS];
% 
% 
% %Use the same quantile defined in OP_INC_july8. No need to re-define mesh
% %here.
% %mesh=quantile(RTotD_NCAC,linspace(0,1,fineness).');
% 
% X_KnotAC1=(RTotD_NCAC<mesh(2,1)).*(1-(RTotD_NCAC-mesh(1,1))/(mesh(2,1)-mesh(1,1)));
% for i=0:(numel(mesh)-3)
%     PLUS=(RTotD_NCAC>=mesh(i+1,1)).*(RTotD_NCAC<mesh(i+2,1)).*((RTotD_NCAC-mesh(i+1,1))/(mesh(i+2,1)-mesh(i+1,1)))...
%         +(RTotD_NCAC>=mesh(i+2,1)).*(RTotD_NCAC<mesh(i+3,1)).*(1-(RTotD_NCAC-mesh(i+2,1))/(mesh(i+3,1)-mesh(i+2,1)));
%     X_KnotAC1=[X_KnotAC1,PLUS];
% end
% PLUS=(RTotD_NCAC>=mesh((numel(mesh)-1),1)).*(RTotD_NCAC-mesh((numel(mesh)-1),1))/(mesh(numel(mesh),1)-mesh((numel(mesh)-1),1));
% X_KnotAC1=[X_KnotAC1,PLUS];
% 
% %First column: redundant
% X_KnotAC1=X_KnotAC1(:,2:fineness);

E_VContestFUL=find(E_V_july8(:,8)==0);           %% E_VContestFUL is the elements in which contest==0 :277 distinct districts
E_VNContestFUL=find(E_V_july8(:,8)==1);          %% E_VNContestFUL is the elements in which contest==1 :297 distinct districts

Win_E_V=(E_V_july8(:,13)>=0.5);  % Dummy for whether incumbent won the election in (t).
Win_E_VCT=Win_E_V;
Win_E_VCT(E_VContestFUL,:)=[];
IND5=find(Win_E_VCT==0);  % index for elections with incumbent loss (among those with entry).
                          % Number of distinct districts within contest=0 is 27.
                          % Number of distinct districts within contest=1 is 295.




PartyE_V=3*ones(length(E_V_july8),1)-2*E_V_july8(:,3);  %%PartyE_V=1  if candidate i is Democrat and PartyE_V=-1 if candidate i is Republican.
% PartyE_VCT=PartyE_V;
% PartyE_VCT(E_VContestFUL,:)=[];     %PartyE_VCT is the party indicator when there is entry next period.
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

%%B-Spline for q_I(RTotDE_V), where RTotDE_V=(NCE_VCT(:,1)./LOGTot_NCE_V)  etc.%%
%% take knots to be between 0.88 to 1.025 with 8 knots (8 basis functions). (Almost all lie
%% within this range) Let X_Knot be the matrix with 8 columns that contain
%% the value of the B-Spline basis function evaluated at each of the 8
%% basis functions.
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

X_KnotE_VNCT=X_KnotEV1;
X_KnotE_VNCT(E_VNContestFUL,:)=[];

%%%%%%%%%%%%%%%%%%%%%%%% 
% Defining variables: done
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%               q_e               %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%thetaS 2 X 1 [unemp,white]  mean of unemp is about 0.07, white is about
%0.81
%thetaQ 9 X 1 [ones,d_NC,w_NC,w'_NXT,d_NC^2,w_NC^2,w'_NXT^2,party.*s_NC,tenure]
%theta  7 X 1 [ones,LOGW_I,XQ,XS.*Party,LOGW_I.^2,XQ.^2,log(T
%             enure.+1)]  
%             mean(LOGD_NCE_V) is 12.4 mean(LOGW_NCE_V) is 11.0,
%             mean(LOGWNXT_NCE_V) is 11.7  mean(LOGD_NXT_E_V) is 12.61
% [QC1,QC2]   [PRENTAC,PRENTAC.^2]
% E_VCTa 15 X 1
% E_VCTt 15 X 1
% gammaCT 15 X 1
% thetawin 7 X 1 [ones,LOGW_NXT_E_VCT,XQEVCT,XS_NXT_EVCT.*PartyE_VCT,LOGW_NXT_E_VCT.^2,XQEVCT.^2,...
%    log(TenureE_VCT.+1)]
% thetadNC 7 X 1 LOGD_NXT_NCT=[1,XQNCT,XQNCT.^2,LOGW_NXT_NCT,LOGW_NXT_NCT.^2,Party_NCT.*XS_NXT_NCT,log(1+TenureNCT)]*thetadNC...
% mean(LOGD_NXT_NCT)==12.30 mean(LOGW_NXT_NCT)==11.25
% thetatNC 7 X 1 LOGTotal_NXT_NCT=[1,XQNCT,XQNCT.^2,LOGW_NXT_NCT,LOGW_NXT_NCT.^2,Party_NCT.*XS_NXT_NCT,log(1+TenureNCT)]*thetatNC
% mean(LOGTotal_NXT_NCT)=12.56  mean(LOGW_NXT_NCT)=11.2534
% thetanxtNC 7 X 1 LOGWNXT_NXT_NCT=[1,XQNCT,XQNCT.^2,LOGW_NXT_NCT,LOGW_NXT_NCT.^2,Party_NCT.*XS_NXT_NCT,log(1+TenureNCT)]*thetanxtNC
% mean(LOGWNXT_NXT_E_VNCT)=11.72, mean(LOGW_NXT_E_VNCT)==11.25
% win=cdfnorm([1,LOGW_NXT_E_VCT,XQEVCT,XS_NXT_EVCT.*PartyE_VCT,LOGW_NXT_E_VCT.^2,XQEVCT.^2,...
%    log(TenureE_VCT+1)]*thetawin)


% % bdds: thetaS(1) in [-3,3], thetaS(2) in [-0.2,0.2]
% % mean(unemployment)*thetaS(1)+mean(white)*thetaS(2) in [-0.35,0.35]
% thetaS=[0.5,-0.05]';
% % bdds: thetaQ(1) in [0.3,0.9] thetaQ(2) in [-0.03,0.03] thetaQ(3) in
% % [-0.03,0.03] thetaQ(4) in [-0.03,0.03] thetaQ(5) in [-0.002,0.002]
% % thetaQ(6) in [-0.002,0.002] thetaQ(7) in [-0.002,0.002] thetaQ(8) in
% % [-0.8,0.8] thetaQ(9) in [-0.2,0.3]
% % overall restriction
% % [1,d_NC,w_NC,w'_NXT,d_NC^2,w_NC^2,w'_NXT^2,party.*s_NC,tenure]*thetaQ in
% % [0.2,0.8]
% thetaQ=[0.65,0.02,0.01,-0.01,-0.0002,-0.0002,0.0001,0.01,0.1]';
% % thata (1) in [-1,2] thata(2) in [-0.3,0.2] theta(3) in [-2,2]
% % theta(4) in [-2,2] theta(5) in [-0.007,0.007] theta(6) in  [-2,2]
% % theta(7) in [-0.5,0.5]
% % overall restriction
% % [1,LOGW_I,XQ,XS.*Party,LOGW_I.^2,XQ.^2,log(Tenure.+1)]*theta in
% % [0,2.5];
% theta=[1.8,-0.015,-0.7,1,0.0001,0.2,-0.07]';
% % 
% E_VCTa=[0.5,0.3,0.1,0.01,0.2,0.5,0.3,0.1,0.01,0.2,0.5,0.3,0.1,0.01,0.2]';
% 
% % E_VCTt(i) in [mean(x_i)-3*std(x_i),mean(x_i)+3*std(x_i)]
% E_VCTt=[mean(LOGD_E_VCT),0.6,mean(LOGW_E_VCT),0,mean(log(1+TenureE_VCT)),mean(LOGTotal_E_VCT),0.6,mean(LOGW_E_VCT),0,log(5)...
%     ,mean(LOGW_NXT_E_VCTwnxt),0.6,mean(LOGW_E_VCT),0,mean(log(1+TenureE_VCT))]';
% 
% lbE_VCTt=E_VCTt-3*[std([LOGD_E_VCT,0.6*ones(length(LOGD_E_VCT),1),LOGW_E_VCT,zeros(length(LOGD_E_VCT),1),...
%     log(1+TenureE_VCT),LOGTotal_E_VCT,0.6*ones(length(LOGD_E_VCT),1),LOGW_E_VCT,zeros(length(LOGD_E_VCT),1),log(1+TenureE_VCT)])';...
%     std(LOGW_NXT_E_VCTwnxt);std([0.6*ones(length(LOGD_E_VCT),1),LOGW_E_VCT,zeros(length(LOGD_E_VCT),1),log(1+TenureE_VCT)])'];
% ubE_VCTt=E_VCTt+3*[std([LOGD_E_VCT,0.6*ones(length(LOGD_E_VCT),1),LOGW_E_VCT,zeros(length(LOGD_E_VCT),1)...
%     ,log(1+TenureE_VCT),LOGTotal_E_VCT,0.6*ones(length(LOGD_E_VCT),1),LOGW_E_VCT,zeros(length(LOGD_E_VCT),1),log(1+TenureE_VCT)])';...
%     std(LOGW_NXT_E_VCTwnxt);std([0.6*ones(length(LOGD_E_VCT),1),LOGW_E_VCT,zeros(length(LOGD_E_VCT),1),log(1+TenureE_VCT)])'];
% lbE_VCTt(2,1)=0.2;
% lbE_VCTt(4,1)=0.2;
% lbE_VCTt(7,1)=0.2;
% lbE_VCTt(9,1)=0.2;
% lbE_VCTt(12,1)=0.2;
% lbE_VCTt(14,1)=0.2;
% ubE_VCTt(2,1)=0.9;
% ubE_VCTt(4,1)=0.9;
% ubE_VCTt(7,1)=0.9;
% ubE_VCTt(9,1)=0.9;
% ubE_VCTt(12,1)=0.9;
% ubE_VCTt(14,1)=0.9;
% 
% 
% 
% % gammaCT(i) in [(1/5)*sigma(x_i),5*sigma(x_i)]
% gammaCT=[0.1,0.14,0.1,0.1,0.2,0.1,0.14,0.1,0.1,0.2,0.1,0.14,0.1,0.1,0.2]';
% 
% lbgammaCT=(1/5)*[std([LOGD_E_VCT,0.6*ones(length(LOGD_E_VCT),1),LOGW_E_VCT,zeros(length(LOGD_E_VCT),1),log(1+TenureE_VCT),...
%     LOGTotal_E_VCT,0.6*ones(length(LOGD_E_VCT),1),LOGW_E_VCT,zeros(length(LOGD_E_VCT),1),log(1+TenureE_VCT)])';...
%     std(LOGW_NXT_E_VCTwnxt);std([0.6*ones(length(LOGD_E_VCT),1),LOGW_E_VCT,zeros(length(LOGD_E_VCT),1),log(1+TenureE_VCT)])'];
% ubgammaCT=5*[std([LOGD_E_VCT,0.6*ones(length(LOGD_E_VCT),1),LOGW_E_VCT,zeros(length(LOGD_E_VCT),1),log(1+TenureE_VCT)...
%     ,LOGTotal_E_VCT,0.6*ones(length(LOGD_E_VCT),1),LOGW_E_VCT,zeros(length(LOGD_E_VCT),1),log(1+TenureE_VCT)])';...
%     std(LOGW_NXT_E_VCTwnxt);std([0.6*ones(length(LOGD_E_VCT),1),LOGW_E_VCT,zeros(length(LOGD_E_VCT),1),log(1+TenureE_VCT)])'];
% lbgammaCT(2,1)=0.05;
% lbgammaCT(4,1)=0.05;
% lbgammaCT(7,1)=0.05;
% lbgammaCT(9,1)=0.05;
% lbgammaCT(12,1)=0.05;
% lbgammaCT(14,1)=0.05;
% ubgammaCT(2,1)=0.9;
% ubgammaCT(4,1)=0.9;
% ubgammaCT(7,1)=0.9;
% ubgammaCT(9,1)=0.9;
% ubgammaCT(12,1)=0.9;
% ubgammaCT(14,1)=0.9;
% 
% % E_VNCTa=[0.5,0.3,0.1,0.01,0.2,0.5,0.3,0.1,0.01,0.2,0.5,0.3,0.1,0.01,0.2]';
% % E_VNCTt=[mean(LOGD_NXT_E_VNCT),0.6,mean(LOGW_NXT_E_VNCT),0,log(5),mean(LOGTotal_NXT_E_VNCT),0.6,mean(LOGW_NXT_E_VNCT),0,log(5)...
% %     ,mean(LOGWNXT_NXT_E_VNCT),0.6,mean(LOGW_NXT_E_VNCT),0,log(5)]';
% % gammaNCT=[0.1,0.14,0.1,0.1,0.2,0.1,0.14,0.1,0.1,0.2,0.1,0.14,0.1,0.1,0.2]';
% 
% % thetadNC(1) in [0,14] thetadNC(2) in [-5,5] thetadNC(3) in [-3,3]
% % thetadNC(4) in [-0.5,1.5] thetadNC(5) in [-0.04,0.04] thetadNC(6) in [-5,5]
% % thetadNC(7) in [-5,5]
% thetadNC=[7,1,-0.1,0.3,0.001,1,0.7]';
% 
% % thetatNC(1) in [0,18] thetatNC(2) in [-5,5] thetatNC(3) in [-3,3]
% % thetatNC(4) in [-1.5,0.5] thetaNC(5) in [-0.04,0.04] thetaNC(6) in [-5,5]
% % thetatNC(7) in [-5,5]
% thetatNC=[15,1,-0.1,-0.3,0.001,1,0.5]';
% 
% % thetanxtNC(1) in [0,18] thetanxtNC(2) in [-5,5] thetanxtNC(3) in [-3,3]
% % thetanxtNC(4) in [-0.5,1.5] thetanxtNC(5) in [-0.04,0.04] thetanxtNC(6) in [-5,5]
% % thetanxtNC(7) in [-5,5]
% thetanxtNC=[9,1,-0.1,0.4,0.001,1,0.15]';
% 
% % win=cdfnorm([1,LOGW_NXT_E_VCT,XQEVCT,XS_NXT_EVCT.*PartyE_VCT,LOGW_NXT_E_VCT.^2,XQEVCT.^2,...
% %    log(TenureE_VCT+1)]*thetawin)
% 
% % thetawin(1) in [-0.5,2] thetawin(2) in [-0.2,0.2] thetawin(3) in [-2,2]
% % thetawin(4) in [-2,2] thetawin(5) in [-0.01,0.01] thetawin(6) in [-2,2]
% % thetawin(7) in [-0.7,0.7]
% % overall restriction [1,LOGW_NXT_E_VCT,XQEVCT,XS_NXT_EVCT.*PartyE_VCT,LOGW_NXT_E_VCT.^2,XQEVCT.^2,...
% %    log(TenureE_VCT+1)]*thetawin in [0.5,0.99]
% thetawin=[0.15,0.030,0.9,1,-0.001,-0.1,0.3]';

% B_I in [-0.01,0.05]
B_I=0.01;
% B_C in [-0.05,0.01]
B_C=-0.015;
% Q_C1 in [-0.5,0.5]
Q_C1=0.06;
% Q_C2 in [-0.5,0.5]
Q_C2=-0.01;
% B_T in [-0.1,0.3]
B_T=0.04;
% overall restriction
% B_I*LOGD_IAC+B_C*LOGD_CAC+XSAC.*PartyAC+XQAC+[PRENTAC,PRENTAC.^2]*[Q_C1,Q
% _C2]'+B_T*log(TenureAC+1) in [0,1.0]

% cost1 in [-0.004,0]
cost1=-0.0013;
% ben1 in [0,0.004]
ben1=0.0005;
% cost2 in [-0.004,0]
cost2=-0.0015;
% ben2 in [0,0.004]
ben2=0.001;
% sig in [0.001,1]
sig=0.3;



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
        XQEVCTwnxt2=XQEVCT2;
        XQEVCTwnxt2(IND5,:)=[];
    else
        XQEV2=X_KnotEV1*thetaQ2;  % q_I
        XQEVCT2=XQEV2;
        XQEVCT2(E_VContestFUL,:)=[];
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
N=10;                      % N is the number of simulations.
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
 C=zeros(5,T,N,length(NCE_V));
% DC=zeros(5,T,N,length(NCE_V));
for i=1:length(NCE_V)
         C(:,:,:,i)=Actions(i,XSEV_(i,:)',SameE_V(i,1), XQEV2(i,1),X_KnotEV1(i,:),RTotDE_V(i,:), TenureE_V(i,1),LOGW_NXT_E_V(i,1),PartyE_V(i,1),PresdumE_V(i,1),MidtermE_V(i,1),...
         Shockpartisan(:,:,i), Shockump(:,:,i), presseq(:,:,i),Entry(:,:,i), coefentry, E_VCTa, E_VCTt, gammaCT, coefspend, coeffund,coefsave, dF_gamma_ct(:,:,i), dF_total_ct(:,:,i),...
       dF_gammasv_ct(:,:,i),dF_totalsv_ct(:,:,i),Winrnd(:,:,i),coefprobwin,Ret(:,i)',Betapartisan,Betaump);
    squeeze(C(:,:,:,i))
    i
%     DC(:,:,:,i)=Actions(i,XSEV_(i,:)',XQEV(i,1),RTotDE_V(i,1),TenureE_V(i,1),LOGW_NXT_E_V(i,1)+Delt,...
%         Shockwh(:,:,i), Shockump(:,:,i), Entry(:,:,i), theta, E_VCTa, E_VCTt, gammaCT, E_VNCTa, E_VNCTt, gammaNCT, dF_gamma_ct(:,:,i), dF_total_ct(:,:,i),...
%         dF_nxt_nxt_ct(:,:,i),dF_gamma_nct(:,:,i),dF_total_nct(:,:,i),dF_nxt_nxt_nct(:,:,i), Winrnd(:,:,i),thetawin,Ret(:,i)',Betawh,Betaump,thetaS,PartyE_V(i,1));
%         if (sum(squeeze(DC(5,1,:,i)-C(5,1,:,i)))<0)
%             i
%             C(:,:,:,i)
%             DC(:,:,:,i)
% 
%         end
%     for k=1:N
%         for j=1:10
%             if C(1,j,k,i)==1 %Contest
%                 Continue1(i,k)=Continue1(i,k)+0.9*(ben1*C(2,j,k,i).^2-cost1*C(3,j,k,i).^2+C(5,j,k,i));
%             else
%                 Continue1(i,k)=Continue1(i,k)+0.9*(ben2*C(2,j,k,i).^2-cost2*C(3,j,k,i).^2+C(5,j,k,i));
%             end
%             if DC(1,j,k,i)==1 %Contest
%                 DContinue1(i,k)=DContinue1(i,k)+0.9*(ben1*DC(2,j,k,i).^2-cost1*DC(3,j,k,i).^2+DC(5,j,k,i));
%             else
%                 DContinue1(i,k)=DContinue1(i,k)+0.9*(ben2*DC(2,j,k,i).^2-cost2*DC(3,j,k,i).^2+DC(5,j,k,i));
%             end
%         end
%     end
%     Continue(i,1)=Continuation(XSEV_(i,:)',XQEV(i,1),TenureE_V(i,1),LOGW_NXT_E_V(i,1),...
%         Shockwh(:,:,i), Shockump(:,:,i), Entry(:,:,i), theta, E_VCTa, E_VCTt, gammaCT, E_VNCTa, E_VNCTt, gammaNCT, dF_gamma_ct(:,:,i), dF_total_ct(:,:,i),...
%         dF_nxt_nxt_ct(:,:,i),dF_gamma_nct(:,:,i),dF_total_nct(:,:,i),dF_nxt_nxt_nct(:,:,i), Winrnd(:,:,i),thetawin,Ret(:,i)',Betawh,Betaump,cost1,ben1,cost2,ben2,thetaS,PartyE_V(i,1));
%     DContinue(i,1)=Continuation(XSEV_(i,:)',XQEV(i,1),TenureE_V(i,1),LOGW_NXT_E_V(i,1)+Delt,...
%         Shockwh(:,:,i), Shockump(:,:,i), Entry(:,:,i), theta, E_VCTa, E_VCTt, gammaCT, E_VNCTa, E_VNCTt, gammaNCT, dF_gamma_ct(:,:,i), dF_total_ct(:,:,i),...
%         dF_nxt_nxt_ct(:,:,i),dF_gamma_nct(:,:,i),dF_total_nct(:,:,i),dF_n
%         xt_nxt_nct(:,:,i), Winrnd(:,:,i),thetawin,Ret(:,i)',Betawh,Betaump,cost1,ben1,cost2,ben2,thetaS,PartyE_V(i,1));
end

  save C.mat C
  %save DC.mat DC
%%
load('./C.mat');
% load('./DC.mat');
mintheta=[];
SRR=[];
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

 
for k=1:2
    theta2initial=Sndstage_initial_long(1:5,k);
    k;
    [mta,SR]=fminsearch(@(theta2) Minimize2S_new(thetain1step,theta2),theta2initial)
    mintheta=[mintheta,mta]
    SRR=[SRR,SR]
end
mintheta
SRR
save minimizedtheta.txt mintheta SRR -ASCII


Qualc(:,1)=E_V_july8(:,31);
Qualc(:,2)=E_V_july8(:,32);
Qualc(:,3)=E_V_july8(:,33);
Qualc(:,4)=E_V_july8(:,34);
Qualc(E_VContestFUL,:)=[];
load q_C_E_VCT.mat
Qualc(:,5)=q_C_E_VCT;
save Qualc.txt Qualc -ASCII
save QualI.txt XQEV2 -ASCII