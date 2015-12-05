clear
% **
% make sure you get the sign of q_e right. Also note that f_d_qe, f_w_qe have the PARTY of Challenger.
% **

% Variables used in number=4
global E_VContestFUL
global TenureE_VCT
global TenureE_VCTwnxt

global Est2
global IND5

global LOGD_E_VCT
global LOGTotal_E_VCT
global LOGW_E_VCT
global SameE_VCT
global XS_EVCT_
global PartyE_VCT
global PresdumE_VCT
global MidtermE_VCT

global LOGW_NXT_E_VCTwnxt
global LOGW_E_VCTwnxt
global SameE_VCTwnxt
global XS_EVCTwnxt_
global PartyE_VCTwnxt
global PresdumE_VCTwnxt
global MidtermE_VCTwnxt

global RTotDE_VCT
global X_KnotEV1







%%%%%%%%%%%%%%%%
%Loading dataset
%%%%%%%%%%%%%%%%

%load ('./E_V_july8.mat')
E_V_july8=csvread('E_Vjuly8partisan.csv',1,0);

% load('./op_inc_july8.mat')
OP_INC_july8=csvread('OP_INC_july8partisan.csv',1,0);

% load('./op_inc_iv_july8.mat')
OP_INC_IV_july8=csvread('OP_INC_IV_july8partisan.csv',1,0);



%Data cleaning
dele1=find(sum(isnan(OP_INC_IV_july8),2)>0);
OP_INC_IV_july8(dele1,:)=[]; %Drop NaN
OP_INC_IV_july8(OP_INC_IV_july8(:,2)>2002,:)=[]; %Drop year 2004 and on
deleA=find(OP_INC_IV_july8(:,16)==0);
OP_INC_IV_july8(deleA,:)=[];%Drop if oppornent disburse=0
OP_INC_IV_july8(102,:)=[];%Drop an Rtotd outlier

dele2=find(sum(isnan(OP_INC_july8),2)>0);
OP_INC_july8(dele2,:)=[]; %Drop NaN
OP_INC_july8(OP_INC_july8(:,2)>2002,:)=[]; %Drop year 2004 and on
deleB=find(((OP_INC_july8(:,16)==0).*OP_INC_july8(:,8))==1);
OP_INC_july8(deleB,:)=[];%Drop if oppornent disburse=0 and contested
OP_INC_july8(119:121,:)=[];%Drop an Rtotd outlier

dele3=find(sum(isnan(E_V_july8),2)>0); 
E_V_july8(dele3,:)=[]; %Drop NaN
E_V_july8(E_V_july8(:,2)>2002,:)=[]; %Drop year 2004 and on
deleC=find((E_V_july8(:,16)==0).*E_V_july8(:,8)==1);
E_V_july8(deleC,:)=[]; %Drop if oppornent disburse=0 and contested
E_V_july8(174:177,:)=[]; %Drop an Rtotd outlier



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Defining Variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%
%OP_INC_july 8 to estimate entry probability
%%%%%%%%%%%%%%%

Estimation_OP=OP_INC_july8; %Sample to estimate entry probability


Samplesize=length(Estimation_OP);
%Year=Estimation_OP(:,2); % year of election.
Contest=Estimation_OP(:,8); % contest
Primary_N=Estimation_OP(:,36);

%LOGD_I=log(max(ones(Samplesize,1),Estimation_OP(:,10))); % spend
LOGW_I=log(max(ones(Samplesize,1),Estimation_OP(:,11))); % w_I
%LOGW_INXT=log(max(ones(Samplesize,1),Estimation_OP(:,12))); % W_I'
%VS=Estimation_OP(:,13); %V_I
%IDN=Estimation_OP(:,21); %Candidate ID
LOGTot_NC=log(max(ones(Samplesize,1),Estimation_OP(:,31))); % Tot_INC
LOGD_NC=log(max(ones(Samplesize,1),Estimation_OP(:,32))); % D_INC
%LOGW_NC=log(max(ones(Samplesize,1),Estimation_OP(:,33))); % W_INC
%LOGWNXT_NC=log(max(ones(Samplesize,1),Estimation_OP(:,34))); %W'_INC
RTotD_NC=(max(0,LOGD_NC).^(-1/2))./LOGTot_NC;  %% ratio of spending over total. (modify)
RTotD_NC=RTotD_NC.*(exp(LOGTot_NC)./exp(LOGD_NC));
% LOGLOGD_NC=log(max(ones(Samplesize,1),LOGD_NC));
% LOGLOGTot_NC=log(max(ones(Samplesize,1),LOGTot_NC));
%Tenure_NC=Estimation_OP(:,35); % Tenure_INC
Tenure=Estimation_OP(:,15);
Party=3*ones(Samplesize,1)-2*Estimation_OP(:,3);  %%Party=1  if candidate i is Democrat and Party=-1 if candidate i is Republican.
Presparty=3*ones(Samplesize,1)-2*Estimation_OP(:,39);
Same=Party==Presparty;
Dif=Party~=Presparty;
Same=Same-Dif;%Same party=1 if candidate party=President party, otherwise -1
% White=Estimation_OP(:,21);
% Black=Estimation_OP(:,22);
% Other=Estimation_OP(:,23);
% White_NC=Estimation_OP(:,27);
% Black_NC=Estimation_OP(:,28);
% Other_NC=Estimation_OP(:,29);
Unemployment=Estimation_OP(:,24);
% Unemployment_NC=Estimation_OP(:,30);
Partisan=Estimation_OP(:,63);

Unempsame=Unemployment.*Same;
Partdemo=Partisan.*Party;
Presdum=Estimation_OP(:,64);
Presdumsame=Presdum.*Same;
Unemppresdumsame=Unemployment.*Presdumsame;
Midterm=Estimation_OP(:,65);
%Midtermsame=Midterm.*Same;
%Unempmidtermsame=Unemployment.*Midtermsame;


%%B-Spline for q_I(LOGLOGD_NC, LOGLOGTot_NC)  %%
% take knots to be between 0.88 to 1.025 with 9 knots (8 basis functions). (Almost all lie
% within this range) Let X_Knot be the matrix with 8 columns that contain
% the value of the B-Spline basis function evaluated at each of the 9
% basis functions.
%---------------------------------------------------------

% mesh=quantile(LOGLOGD_NC,[.125;.25;.375;.5;.625;.75;.875]);
% mesh=[min(LOGLOGD_NC);mesh;max(LOGLOGD_NC)];
% X_Knot1=(LOGLOGD_NC<mesh(2,1)).*(1-(LOGLOGD_NC-mesh(1,1))/(mesh(2,1)-mesh(1,1)));
% for i=0:6
%     PLUS=(LOGLOGD_NC>=mesh(i+1,1)).*(LOGLOGD_NC<mesh(i+2,1)).*((LOGLOGD_NC-mesh(i+1,1))/(mesh(i+2,1)-mesh(i+1,1)))...
%         +(LOGLOGD_NC>mesh(i+2,1)).*(LOGLOGD_NC<mesh(i+3,1)).*(1-(LOGLOGD_NC-mesh(i+2,1))/(mesh(i+3,1)-mesh(i+2,1)));
%     X_Knot1=[X_Knot1,PLUS];
% end
% PLUS=(LOGLOGD_NC>=mesh(8,1)).*(LOGLOGD_NC-mesh(8,1))/(mesh(9,1)-mesh(8,1));
% X_Knot1=[X_Knot1,PLUS];
% 
% mesh=quantile(LOGLOGTot_NC,[.125;.25;.375;.5;.625;.75;.875]);
% mesh=[min(LOGLOGTot_NC);mesh;max(LOGLOGTot_NC)];
% X_Knot2=(LOGLOGTot_NC<mesh(2,1)).*(1-(LOGLOGTot_NC-mesh(1,1))/(mesh(2,1)-mesh(1,1)));
% for i=0:6
%     PLUS=(LOGLOGTot_NC>=mesh(i+1,1)).*(LOGLOGTot_NC<mesh(i+2,1)).*((LOGLOGTot_NC-mesh(i+1,1))/(mesh(i+2,1)-mesh(i+1,1)))...
%         +(LOGLOGTot_NC>mesh(i+2,1)).*(LOGLOGTot_NC<mesh(i+3,1)).*(1-(LOGLOGTot_NC-mesh(i+2,1))/(mesh(i+3,1)-mesh(i+2,1)));
%     X_Knot2=[X_Knot2,PLUS];
% end
% PLUS=(LOGLOGTot_NC>=mesh(8,1)).*(LOGLOGTot_NC-mesh(8,1))/(mesh(9,1)-mesh(8,1));
% X_Knot2=[X_Knot2,PLUS];
fineness=9;
mesh=quantile(RTotD_NC,linspace(0,1,fineness).');
X_Knot1=(RTotD_NC<mesh(2,1)).*(1-(RTotD_NC-mesh(1,1))/(mesh(2,1)-mesh(1,1)));
for i=0:(numel(mesh)-3)
    PLUS=(RTotD_NC>=mesh(i+1,1)).*(RTotD_NC<mesh(i+2,1)).*((RTotD_NC-mesh(i+1,1))/(mesh(i+2,1)-mesh(i+1,1)))...
        +(RTotD_NC>=mesh(i+2,1)).*(RTotD_NC<mesh(i+3,1)).*(1-(RTotD_NC-mesh(i+2,1))/(mesh(i+3,1)-mesh(i+2,1)));
    X_Knot1=[X_Knot1,PLUS];
end
PLUS=(RTotD_NC>=mesh((numel(mesh)-1),1)).*(RTotD_NC-mesh((numel(mesh)-1),1))/(mesh(numel(mesh),1)-mesh((numel(mesh)-1),1));
X_Knot1=[X_Knot1,PLUS];

%First column:redundant
X_Knot1=X_Knot1(:,2:fineness);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%      Ai & Chen     %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AC=OP_INC_IV_july8;%sample to estimate vote share

%Variables used in Number=2, estimation of vote share
SamplesizeAC=length(OP_INC_IV_july8);
LOGD_IAC=log(max(ones(SamplesizeAC,1),AC(:,10))); % spend
LOGW_IAC=log(max(ones(SamplesizeAC,1),AC(:,11))); % w_I
%LOGW_INXTAC=log(max(ones(SamplesizeAC,1),AC(:,12))); % W_I'
VSAC=AC(:,13); %V_I
LOGD_NCAC=log(max(ones(SamplesizeAC,1),AC(:,32))); % D_INC
%LOGW_NCAC=log(max(ones(SamplesizeAC,1),AC(:,33))); % W_INC
LOGTot_NCAC=log(max(ones(SamplesizeAC,1),AC(:,31))); % Tot_INC
% LOGLOGD_NCAC=log(max(ones(SamplesizeAC,1),LOGD_NCAC));
% LOGLOGTot_NCAC=log(max(ones(SamplesizeAC,1),LOGTot_NCAC));
RTotD_NCAC=(max(0,LOGD_NCAC).^(-1/2))./LOGTot_NCAC; %% D_INC./Tot_INC
RTotD_NCAC=RTotD_NCAC.*(exp(LOGTot_NCAC)./exp(LOGD_NCAC));
%LOGWNXT_NCAC=log(max(ones(SamplesizeAC,1),AC(:,34))); %W'_INC
%Tenure_NCAC=AC(:,35); % Tenure_INC
%NCAC=[LOGD_NCAC,LOGW_NCAC,LOGWNXT_NCAC];
LOGD_CAC=log(max(ones(SamplesizeAC,1),AC(:,16))); %Sepnding of Challenger
% LOGW_CAC=log(max(ones(SamplesizeAC,1),AC(:,17))); %Warchest of Challenger
% LOGW_CNXTAC=log(max(ones(SamplesizeAC,1),AC(:,18))); %Savings of Challenger
% YearAC=AC(:,2); % Current year.

TenureAC=AC(:,15);
PartyAC=3*ones(SamplesizeAC,1)-2*AC(:,3);  %%PartyAC=1  if candidate i is Democrat and PartyAC=-1 if candidate i is Republican.
PrespartyAC=3*ones(SamplesizeAC,1)-2*AC(:,39);
SameAC=PartyAC==PrespartyAC; %Same party if the two match
DifAC=PartyAC~=PrespartyAC;
SameAC=SameAC-DifAC;

% WhiteAC=AC(:,21);
% BlackAC=AC(:,22);
% OtherAC=AC(:,23);
% White_NCAC=AC(:,27);
% Black_NCAC=AC(:,28);
% Other_NCAC=AC(:,29);
UnemploymentAC=AC(:,24);
UnemploymentsqAC=UnemploymentAC.^2;
%Unemployment_NCAC=AC(:,30);
PartisanAC=AC(:,63);


UnempsameAC=UnemploymentAC.*SameAC;
UnempsqsameAC=UnemploymentsqAC.*SameAC;
PartdemoAC=PartisanAC.*PartyAC;

PresdumAC=AC(:,64);
PresdumsameAC=PresdumAC.*SameAC;
UnemppresdumsameAC=UnemploymentAC.*PresdumsameAC;
MidtermAC=AC(:,65);
%MidtermsameAC=MidtermAC.*SameAC;
%UnempmidtermsameAC=UnemploymentAC.*MidtermsameAC;

%Variables used in number=3, estimation of winning probability
Win_AC=(AC(:,13)>=0.5);  % Dummy for whether incumbent won the election in (t).


%Where do we use them?
% EPS=randn(SamplesizeAC,NumSim);
% EPS2=randn(SamplesizeAC,NumSim);
% EPS3=randn(SamplesizeAC,NumSim);


%%B-Spline for q_I(RTotD_NCAC), where RTotD_NCAC=(LOGD_NCAC./LOGTot_NCAC)  etc.%%
%% take knots to be between 0.88 to 1.025 with 8 knots (8 basis functions). (Almost all lie
% within this range) Let X_Knot be the matrix with 8 columns that contain
% the value of the B-Spline basis function evaluated at each of the 8
% basis functions.
% 
% mesh=quantile(LOGLOGD_NCAC,[.125;.25;.375;.5;.625;.75;.875]);
% mesh=[min(LOGLOGD_NCAC);mesh;max(LOGLOGD_NCAC)];
% X_KnotAC1=(LOGLOGD_NCAC<mesh(2,1)).*(1-(LOGLOGD_NCAC-mesh(1,1))/(mesh(2,1)-mesh(1,1)));
% for i=0:6
%     PLUS=(LOGLOGD_NCAC>=mesh(i+1,1)).*(LOGLOGD_NCAC<mesh(i+2,1)).*((LOGLOGD_NCAC-mesh(i+1,1))/(mesh(i+2,1)-mesh(i+1,1)))...
%         +(LOGLOGD_NCAC>mesh(i+2,1)).*(LOGLOGD_NCAC<mesh(i+3,1)).*(1-(LOGLOGD_NCAC-mesh(i+2,1))/(mesh(i+3,1)-mesh(i+2,1)));
%     X_KnotAC1=[X_KnotAC1,PLUS];
% end
% PLUS=(LOGLOGD_NCAC>=mesh(8,1)).*(LOGLOGD_NCAC-mesh(8,1))/(mesh(9,1)-mesh(8,1));
% X_KnotAC1=[X_KnotAC1,PLUS];
% 
% mesh=quantile(LOGLOGTot_NCAC,[.125;.25;.375;.5;.625;.75;.875]);
% mesh=[min(LOGLOGTot_NCAC);mesh;max(LOGLOGTot_NCAC)];
% X_KnotAC2=(LOGLOGTot_NCAC<mesh(2,1)).*(1-(LOGLOGTot_NCAC-mesh(1,1))/(mesh(2,1)-mesh(1,1)));
% for i=0:6
%     PLUS=(LOGLOGTot_NCAC>=mesh(i+1,1)).*(LOGLOGTot_NCAC<mesh(i+2,1)).*((LOGLOGTot_NCAC-mesh(i+1,1))/(mesh(i+2,1)-mesh(i+1,1)))...
%         +(LOGLOGTot_NCAC>mesh(i+2,1)).*(LOGLOGTot_NCAC<mesh(i+3,1)).*(1-(LOGLOGTot_NCAC-mesh(i+2,1))/(mesh(i+3,1)-mesh(i+2,1)));
%     X_KnotAC2=[X_KnotAC2,PLUS];
% end
% PLUS=(LOGLOGTot_NCAC>=mesh(8,1)).*(LOGLOGTot_NCAC-mesh(8,1))/(mesh(9,1)-mesh(8,1));
% X_KnotAC2=[X_KnotAC2,PLUS];

%Use the same quantile defined in OP_INC_july8. No need to re-define mesh
%here.
%mesh=quantile(RTotD_NCAC,linspace(0,1,fineness).');

X_KnotAC1=(RTotD_NCAC<mesh(2,1)).*(1-(RTotD_NCAC-mesh(1,1))/(mesh(2,1)-mesh(1,1)));
for i=0:(numel(mesh)-3)
    PLUS=(RTotD_NCAC>=mesh(i+1,1)).*(RTotD_NCAC<mesh(i+2,1)).*((RTotD_NCAC-mesh(i+1,1))/(mesh(i+2,1)-mesh(i+1,1)))...
        +(RTotD_NCAC>=mesh(i+2,1)).*(RTotD_NCAC<mesh(i+3,1)).*(1-(RTotD_NCAC-mesh(i+2,1))/(mesh(i+3,1)-mesh(i+2,1)));
    X_KnotAC1=[X_KnotAC1,PLUS];
end
PLUS=(RTotD_NCAC>=mesh((numel(mesh)-1),1)).*(RTotD_NCAC-mesh((numel(mesh)-1),1))/(mesh(numel(mesh),1)-mesh((numel(mesh)-1),1));
X_KnotAC1=[X_KnotAC1,PLUS];

%First column: redundant
X_KnotAC1=X_KnotAC1(:,2:fineness);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%               E_V               %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
PartyE_VCTwnxt=PartyE_VCT;
PartyE_VCTwnxt(IND5,:)=[];
PartyE_VNCT=PartyE_V;
PartyE_VNCT(E_VNContestFUL,:)=[];    %PartyE_VCT is the party indicator when there is NO entry next period.

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
% YEARE_V=E_V_july8(:,2);           % YEARE_V is the year of the election.
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

%%B-Spline for q_I(RTotDE_VCT), where RTotDE_VCT=(NCE_VCT(:,1)./LOGTotal_E_VCT)  etc.%%
%% take knots to be between 0.88 to 1.025 with 8 knots (8 basis functions). (Almost all lie
% within this range) Let X_Knot be the matrix with 8 columns that contain
% the value of the B-Spline basis function evaluated at each of the 8
% basis functions.


% 
% mesh=quantile(LOGLOGD_E_V,[.125;.25;.375;.5;.625;.75;.875]);
% mesh=[min(LOGLOGD_E_V);mesh;max(LOGLOGD_E_V)];
% X_KnotEV1=(LOGLOGD_E_V<mesh(2,1)).*(1-(LOGLOGD_E_V-mesh(1,1))/(mesh(2,1)-mesh(1,1)));
% for i=0:6
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
% for i=0:6
%     PLUS=(LOGLOGTot_E_V>=mesh(i+1,1)).*(LOGLOGTot_E_V<mesh(i+2,1)).*((LOGLOGTot_E_V-mesh(i+1,1))/(mesh(i+2,1)-mesh(i+1,1)))...
%         +(LOGLOGTot_E_V>mesh(i+2,1)).*(LOGLOGTot_E_V<mesh(i+3,1)).*(1-(LOGLOGTot_E_V-mesh(i+2,1))/(mesh(i+3,1)-mesh(i+2,1)));
%     X_KnotEV2=[X_KnotEV2,PLUS];
% end
% PLUS=(LOGLOGTot_E_V>=mesh(8,1)).*(LOGLOGTot_E_V-mesh(8,1))/(mesh(9,1)-mesh(8,1));
% X_KnotEV2=[X_KnotEV2,PLUS];

%mesh=quantile(RTotDE_V,linspace(0,1,fineness).');

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




%%
%%%%%%%%%%%%%%%%%
%Estimation of entry probability and Primary N
%Number=1
%%%%%%%%%%%%%%%%%



%Specification 1
 Xentry=[ones(Samplesize,1),LOGW_I,Unempsame,Partdemo,log(Tenure+1),X_Knot1];
 coefentry=Xentry\Contest;
 coefprimaryN=Xentry\Primary_N;

%Specification 2
%  Xentry2=[ones(Samplesize,1),LOGW_I,Unempsame,Partdemo,log(Tenure+1),Presdum,Presdumsame,Unemppresdumsame,Midterm,X_Knot1,X_Knot1.*(LOGW_I*ones(1,fineness-1))];
%  coefentry2=Xentry2\Contest;
%  coefprimaryN2=Xentry2\Primary_N;



%%%%%%%%%%
%Result of specification check: available on dropbox.

%%
%%%%%%%%%%%%%%%%%
%Estimation of Vote share equation
%Number=2
%%%%%%%%%%%%%%%%%


Eentry=[ones(SamplesizeAC,1),LOGW_IAC,UnempsameAC,PartdemoAC,log(TenureAC+1),X_KnotAC1]*coefentry;
Eentry2=min(1,[ones(SamplesizeAC,1),LOGW_IAC,UnempsameAC,PartdemoAC,log(TenureAC+1),X_KnotAC1]*coefentry);
EprimaryN=[ones(SamplesizeAC,1),LOGW_IAC,UnempsameAC,PartdemoAC,log(TenureAC+1),X_KnotAC1]*coefprimaryN;

Xvoteshare=[LOGD_IAC,LOGD_CAC,UnempsameAC,PartdemoAC,log(TenureAC+1),EprimaryN,Eentry,EprimaryN.*Eentry,RTotD_NCAC,RTotD_NCAC.^2,RTotD_NCAC.^3];
IVvoteshare=[ones(length(VSAC),1),LOGW_IAC,LOGW_IAC.^2,LOGW_IAC.*log(TenureAC+1),UnempsameAC,UnempsqsameAC,PartdemoAC,log(TenureAC+1),(log(TenureAC+1)).^2,X_KnotAC1,X_KnotAC1.*(LOGW_IAC*ones(1,fineness-1)),PresdumAC,PresdumsameAC,UnemppresdumsameAC,MidtermAC];
IVvar=IVvoteshare.'*IVvoteshare;

coefvoteshare=inv(Xvoteshare.'*IVvoteshare*inv(IVvar)*IVvoteshare.'*Xvoteshare)*(Xvoteshare.'*IVvoteshare*inv(IVvar)*IVvoteshare.'*(VSAC-0.5));


% Eentry2=[ones(SamplesizeAC,1),LOGW_IAC,UnempsameAC,PartdemoAC,log(TenureAC+1),PresdumAC,PresdumsameAC,UnemppresdumsameAC,MidtermAC,X_KnotAC1,X_KnotAC1.*(LOGW_IAC*ones(1,fineness-1))]*coefentry2;
% EprimaryN2=[ones(SamplesizeAC,1),LOGW_IAC,UnempsameAC,PartdemoAC,log(TenureAC+1),PresdumAC,PresdumsameAC,UnemppresdumsameAC,MidtermAC,X_KnotAC1,X_KnotAC1.*(LOGW_IAC*ones(1,fineness-1))]*coefprimaryN2;
% Xvoteshare2=[LOGD_IAC,LOGD_CAC,UnempsameAC,PartdemoAC,log(TenureAC+1),EprimaryN2,Eentry2,EprimaryN2.*Eentry2,(EprimaryN2).^2,X_KnotAC1];
% IVvoteshare2=[ones(length(VSAC),1),LOGW_IAC,LOGW_IAC.^2,LOGW_IAC.*log(TenureAC+1),UnempsameAC,UnempsqsameAC,PartdemoAC,log(TenureAC+1),(log(TenureAC+1)).^2,X_KnotAC1,X_KnotAC1.*(LOGW_IAC*ones(1,fineness-1)),PresdumAC,PresdumsameAC,UnemppresdumsameAC,MidtermAC];
% IVvar2=IVvoteshare2.'*IVvoteshare2;
% 
% coefvoteshare2=inv(Xvoteshare2.'*IVvoteshare2*inv(IVvar2)*IVvoteshare2.'*Xvoteshare2)*(Xvoteshare2.'*IVvoteshare2*inv(IVvar2)*IVvoteshare2.'*(VSAC-0.5));
% % 
%%
%%%%%%%%%%%%%%%%%%
%Estimation of winning probability
%Number=3
%%%%%%%%%%%%%%%%%%
Xprobwin=[ones(length(Win_AC),1),LOGW_IAC,UnempsameAC,UnempsqsameAC,PartdemoAC,LOGW_IAC.^2,X_KnotAC1,log(TenureAC+1),log(TenureAC+1).^2,...
    X_KnotAC1.*(UnempsameAC*ones(1,8)),X_KnotAC1.*(PartdemoAC*ones(1,8)),X_KnotAC1.*(log(TenureAC+1)*ones(1,8)),...;
  LOGW_IAC.^3];%X_KnotAC1.*(LOGW_IAC*ones(1,8)), ,log(TenureAC+1).*LOGW_IAC

%X_KnotAC1.*(LOGW_IAC*ones(1,8)), RTotD_NCAC.*LOGW_IAC,LOGW_IAC.^2, LOGW_IAC.^3,log(TenureAC+1).*LOGW_IAC,
    %Xprobwin=[ones(length(Win_AC),1),LOGW_IAC,UnempsameAC,PartdemoAC,LOGW_IAC.^2,X_KnotAC1,log(TenureAC+1)];X_KnotAC1.*(LOGW_IAC*ones(1,8)),
%LOGW_IAC.^3,
coefprobwin=Xprobwin\Win_AC;
predprobwin=min(1,Xprobwin*coefprobwin);

%If modifying Xprobwin, make sure the corresponding part on "actions.m" is
%modified as well.

%%
%Specification check for winning probability
errorprobwin=sum((Win_AC-Xprobwin*coefprobwin).^2)/(var(Win_AC)*length(Win_AC));
plotprobwin=sortrows([Win_AC,predprobwin]);
scatter(plotprobwin(:,1),plotprobwin(:,2));

%%
%Specification check2: if probability of winning is increasing in beginning
%cash(logW_IAC) conditional on entry
derivprobwin=@(x,logtenure,rtotd) [zeros(numel(x),1),ones(numel(x),1),zeros(numel(x),3),2*x,zeros(numel(x),10),...
    zeros(numel(x),24),3*x.^2]*coefprobwin;%rtotd*ones(numel(x),1),
xspace=linspace(5,14.3,1000).';

scatter(xspace,derivprobwin(xspace,max(log(TenureAC+1)),quantile(RTotD_NCAC,0.9)));

%%
%Conditional winning probability negative w.r.t begcash.
%Check if overall probability of winning increasing in begcash
derivprobwin2=[zeros(numel(LOGW_IAC),1),ones(numel(LOGW_IAC),1),zeros(numel(LOGW_IAC),3),2*LOGW_IAC,zeros(numel(LOGW_IAC),10),...
    zeros(numel(LOGW_IAC),24),3*LOGW_IAC.^2]*coefprobwin;%2*LOGW_IAC,,3*LOGW_IAC.^2

derivprobwinall1=(coefentry(2)).*predprobwin+Eentry2.*derivprobwin2-(coefentry(2))*ones(numel(LOGW_IAC),1);
hist(derivprobwinall1(derivprobwinall1>-0.01),50)

test1=find(derivprobwinall1<0);
test=[derivprobwinall1(test1),test1,RTotD_NCAC(test1),LOGW_IAC(test1)];

%%
%%%%%%%%%%%%%%%%%%
%Define variables to feed into number=4
%%%%%%%%%%%%%%%%%%
Est1=[coefentry;coefprimaryN];
Est2=coefvoteshare;
%Est1=[coefentry2;coefprimaryN2];
%Est2=coefvoteshare2;


%%
%%%%%%%%%%%%%%%%%%
%Estimation of distribution of spending, saving and fund-raising when
%contested
%Number=4
%%%%%%%%%%%%%%%%%%

  
  %%%%%%%%%%
  %4-1: spending
  %%%%%%%%%%
load('./inittheta1.mat');
Sofarbest1=10^8;
Sofarbest2=10^8;
Sofarbest3=10^8;
  OPTIONS=optimset('MaxIter',10000,'MaxFunEvals',1000000,'Display','iter');
  
%%
for j=1:500
    
%Find initial theta0
    thetaQ2=Est2(9:numel(Est2));
if numel(thetaQ2)==3
    XQEV2=[RTotDE_VCT,RTotDE_VCT.^2,RTotDE_VCT.^3]*thetaQ2; % q_I
    XQEVCT2=XQEV2;
    %XQEVCT2(E_VContestFUL,:)=[]; Already VCT
    XQEVCTwnxt2=XQEVCT2;
    XQEVCTwnxt2(IND5,:)=[];
else
    XQEV2=X_KnotEV1*thetaQ2;  % q_I
    XQEVCT2=XQEV2;
    XQEVCT2(E_VContestFUL,:)=[];
    XQEVCTwnxt2=XQEVCT2;
    XQEVCTwnxt2(IND5,:)=[];
end

theta0=zeros(24,1);    
theta0(1:8)=5*rand(8,1);
theta0(9)=mean(LOGD_E_VCT);
theta0(10)=mean(XQEVCT2);
theta0(11)=mean(LOGW_E_VCT);
theta0(12)=mean(XS_EVCT_(:,1));
theta0(13)=mean(XS_EVCT_(:,2).*PartyE_VCT);
theta0(14)=mean(SameE_VCT);
theta0(15)=mean(PresdumE_VCT);
theta0(16)=mean(MidtermE_VCT);
theta0(17:24)=rand(8,1);


%Saving for young

funvalue=Minimize_Apr2013(theta0,11);
funvalueend=0;
theta411=theta0;
while abs(funvalue-funvalueend)>0.001
    funvalue=Minimize_Apr2013(theta411,11);
    [mintheta11,SRR]=fminsearch(@(x) Minimize_Apr2013(x,11),theta411,OPTIONS);
    funvalueend=Minimize_Apr2013(mintheta11,11);
%      if mod(sss,30)==1
%      theta0=mintheta+0.5*(rand(size(mintheta,1),1)-0.5).*mintheta;
%      theta0([9,17])=mintheta([9,17])+ones(2,1)*0.3;
%      else
        theta411=mintheta11;
%     end
    
    if SRR<Sofarbest1
        Sofarbest1=SRR;
        %    bestiter=iterate;
        %save Est4.txt mintheta -ASCII
        aux11=mintheta11;
    end
end

Est411=aux11;
save Est411 Est411

funvalue=Minimize_Apr2013(theta0,12);
funvalueend=0;
theta412=theta0;
while abs(funvalue-funvalueend)>0.001
    funvalue=Minimize_Apr2013(theta412,12);
    [mintheta12,SRR]=fminsearch(@(x) Minimize_Apr2013(x,12),theta412,OPTIONS);
    funvalueend=Minimize_Apr2013(mintheta12,12);
%      if mod(sss,30)==1
%      theta0=mintheta+0.5*(rand(size(mintheta,1),1)-0.5).*mintheta;
%      theta0([9,17])=mintheta([9,17])+ones(2,1)*0.3;
%      else
        theta412=mintheta12;
%     end
    
    if SRR<Sofarbest2
        Sofarbest2=SRR;
        %    bestiter=iterate;
        %save Est4.txt mintheta -ASCII
        aux12=mintheta12;
    end
end

Est412=mintheta12;
save Est412 Est412

funvalue=Minimize_Apr2013(theta0,13);
funvalueend=0;
theta413=theta0;
while abs(funvalue-funvalueend)>0.001
    funvalue=Minimize_Apr2013(theta413,13);
    [mintheta13,SRR]=fminsearch(@(x) Minimize_Apr2013(x,13),theta413,OPTIONS);
    funvalueend=Minimize_Apr2013(mintheta13,13);
%      if mod(sss,30)==1
%      theta0=mintheta+0.5*(rand(size(mintheta,1),1)-0.5).*mintheta;
%      theta0([9,17])=mintheta([9,17])+ones(2,1)*0.3;
%      else
        theta413=mintheta13;
%     end
    
    if SRR<Sofarbest3
        Sofarbest3=SRR;
        %    bestiter=iterate;
        %save Est4.txt mintheta -ASCII
        aux13=mintheta13;
    end
end

Est413=mintheta13;
save Est413 Est413
end

% for sss=1:1000
%     [mintheta,SRR]=fminsearch(@(x) Minimize_Apr2013(x,1),theta0,OPTIONS);
%     theta0=mintheta;
%     
%     if SRR<Sofarbest
%         Sofarbest=SRR;
%         %    bestiter=iterate;
%         %save Est4.txt mintheta -ASCII
%         save aux1 mintheta
%     end
% end
% 
% load ('aux1.mat');
% Est41=mintheta;
% save Est41 Est41
%%
  %%%%%%%%%%
  %4-2: fundraising
  %%%%%%%%%%
  Sofarbest=10^8;
  OPTIONS=optimset('MaxIter',10000,'MaxFunEvals',1000000,'Display','iter');
  
load('./inittheta2.mat');
theta0=inittheta2;

for sss=1:1000
    [mintheta,SRR]=fminsearch(@(x) Minimize_Apr2013(x,2),theta0,OPTIONS);
    theta0=mintheta;
    
    if SRR<Sofarbest
        Sofarbest=SRR;
        %    bestiter=iterate;
        %save Est4.txt mintheta -ASCII
        save aux2 mintheta
    end
end

load ('aux2.mat');
Est42=mintheta;
save Est42 Est42
%%
  %%%%%%%%%%
  %4-3: saving
  %%%%%%%%%%
  load('inittheta431.mat');
    load('inittheta432.mat');
      load('inittheta433.mat');
  Sofarbest1=Minimize_Apr2013(inittheta431,31);
   Sofarbest2=Minimize_Apr2013(inittheta432,32);
    Sofarbest3=Minimize_Apr2013(inittheta433,33);
  OPTIONS=optimset('MaxIter',10000,'MaxFunEvals',1000000);

for j=1:500
    
%Find initial theta0
    thetaQ2=Est2(9:numel(Est2));
if numel(thetaQ2)==3
    XQEV2=[RTotDE_VCT,RTotDE_VCT.^2,RTotDE_VCT.^3]*thetaQ2; % q_I
    XQEVCT2=XQEV2;
    %XQEVCT2(E_VContestFUL,:)=[]; Already VCT
    XQEVCTwnxt2=XQEVCT2;
    XQEVCTwnxt2(IND5,:)=[];
else
    XQEV2=X_KnotEV1*thetaQ2;  % q_I
    XQEVCT2=XQEV2;
    XQEVCT2(E_VContestFUL,:)=[];
    XQEVCTwnxt2=XQEVCT2;
    XQEVCTwnxt2(IND5,:)=[];
end

theta0=zeros(24,1);    
theta0(1:8)=5*rand(8,1);
theta0(9)=mean(LOGW_NXT_E_VCTwnxt);
theta0(10)=mean(XQEVCTwnxt2);
theta0(11)=mean(LOGW_E_VCTwnxt);
theta0(12)=mean(XS_EVCTwnxt_(:,1));
theta0(13)=mean(XS_EVCTwnxt_(:,2).*PartyE_VCTwnxt);
theta0(14)=mean(SameE_VCTwnxt);
theta0(15)=mean(PresdumE_VCTwnxt);
theta0(16)=mean(MidtermE_VCTwnxt);
theta0(17:24)=rand(8,1);


%Saving for young

funvalue=Minimize_Apr2013(theta0,31);
funvalueend=0;
theta431=theta0;
while abs(funvalue-funvalueend)>0.001
    funvalue=Minimize_Apr2013(theta431,31);
    [mintheta31,SRR]=fminsearch(@(x) Minimize_Apr2013(x,31),theta431,OPTIONS);
    funvalueend=Minimize_Apr2013(mintheta31,31);
%      if mod(sss,30)==1
%      theta0=mintheta+0.5*(rand(size(mintheta,1),1)-0.5).*mintheta;
%      theta0([9,17])=mintheta([9,17])+ones(2,1)*0.3;
%      else
        theta431=mintheta31;
%     end
    
    if SRR<Sofarbest1
        Sofarbest1=SRR;
        %    bestiter=iterate;
        %save Est4.txt mintheta -ASCII
        save aux31 mintheta31
    end
end

load ('aux31.mat');
Est431=mintheta31;
save Est431 Est431

funvalue=Minimize_Apr2013(theta0,32);
funvalueend=0;
theta432=theta0;
while abs(funvalue-funvalueend)>0.001
    funvalue=Minimize_Apr2013(theta432,32);
    [mintheta32,SRR]=fminsearch(@(x) Minimize_Apr2013(x,32),theta432,OPTIONS);
    funvalueend=Minimize_Apr2013(mintheta32,32);
%      if mod(sss,30)==1
%      theta0=mintheta+0.5*(rand(size(mintheta,1),1)-0.5).*mintheta;
%      theta0([9,17])=mintheta([9,17])+ones(2,1)*0.3;
%      else
        theta432=mintheta32;
%     end
    
    if SRR<Sofarbest2
        Sofarbest2=SRR;
        %    bestiter=iterate;
        %save Est4.txt mintheta -ASCII
        save aux32 mintheta32
    end
end

load ('aux32.mat');
Est432=mintheta32;
save Est432 Est432

funvalue=Minimize_Apr2013(theta0,33);
funvalueend=0;
theta433=theta0;
while abs(funvalue-funvalueend)>0.001
    funvalue=Minimize_Apr2013(theta433,33);
    [mintheta33,SRR]=fminsearch(@(x) Minimize_Apr2013(x,33),theta433,OPTIONS);
    funvalueend=Minimize_Apr2013(mintheta33,33);
%      if mod(sss,30)==1
%      theta0=mintheta+0.5*(rand(size(mintheta,1),1)-0.5).*mintheta;
%      theta0([9,17])=mintheta([9,17])+ones(2,1)*0.3;
%      else
        theta433=mintheta33;
%     end
    
    if SRR<Sofarbest3
        Sofarbest3=SRR;
        %    bestiter=iterate;
        %save Est4.txt mintheta -ASCII
        save aux33 mintheta33
    end
end

load ('aux33.mat');
Est433=mintheta33;
save Est433 Est433
end


%%
%when not computing number=4
load('Est41.mat');
load('Est42.mat');
load('Est431.mat');
load('Est432.mat');
load('Est433.mat');

%Specification check: See if actions have right derivative
spendderiv=[Est41(1),Est41(3);Est41(9),Est41(11);Est41(17),Est41(19)];
fundderiv=[Est42(1),Est42(3);Est42(9),Est42(11);Est42(17),Est42(19)];
savederiv=[Est431(1),Est431(3);Est42(1),Est432(3);Est433(1),Est433(3)];

%Correlations on data
corr([LOGW_E_VCT(TenureE_VCT>quantile(TenureE_VCT,0.7)),LOGD_E_VCT(TenureE_VCT>quantile(TenureE_VCT,0.7)),LOGTotal_E_VCT(TenureE_VCT>quantile(TenureE_VCT,0.7)),LOGW_NXT_E_VCT(TenureE_VCT>quantile(TenureE_VCT,0.7))])
corr([LOGW_E_VCT(TenureE_VCT<quantile(TenureE_VCT,0.3)),LOGD_E_VCT(TenureE_VCT<quantile(TenureE_VCT,0.3)),LOGTotal_E_VCT(TenureE_VCT<quantile(TenureE_VCT,0.3)),LOGW_NXT_E_VCT(TenureE_VCT<quantile(TenureE_VCT,0.3))])

corr([LOGW_E_VCT(TenureE_VCT>quantile(TenureE_VCT,0.7)&RTotDE_VCT<quantile(RTotDE_VCT,0.3)),LOGD_E_VCT(TenureE_VCT>quantile(TenureE_VCT,0.7)&RTotDE_VCT<quantile(RTotDE_VCT,0.3)),LOGTotal_E_VCT(TenureE_VCT>quantile(TenureE_VCT,0.7)&RTotDE_VCT<quantile(RTotDE_VCT,0.3)),LOGW_NXT_E_VCT(TenureE_VCT>quantile(TenureE_VCT,0.7)&RTotDE_VCT<quantile(RTotDE_VCT,0.3))])
corr([LOGW_E_VCT(TenureE_VCT>quantile(TenureE_VCT,0.7)&RTotDE_VCT>quantile(RTotDE_VCT,0.8)),LOGD_E_VCT(TenureE_VCT>quantile(TenureE_VCT,0.7)&RTotDE_VCT>quantile(RTotDE_VCT,0.8)),LOGTotal_E_VCT(TenureE_VCT>quantile(TenureE_VCT,0.7)&RTotDE_VCT>quantile(RTotDE_VCT,0.8)),LOGW_NXT_E_VCT(TenureE_VCT>quantile(TenureE_VCT,0.7)&RTotDE_VCT>quantile(RTotDE_VCT,0.8))])

corr([LOGW_E_VCT,LOGD_E_VCT,LOGTotal_E_VCT,LOGW_NXT_E_VCT])

%%
%%%%%%%%%%%%%%%%%%
%Estimation of distribution of spending, saving and fund-raising when
%not contested
%Number=5
%%%%%%%%%%%%%%%%%%


%Each actions

%X_Knot version...Derivative not necessarily positive
% Xaction=[ones(size(LOGD_E_VNCT,1),1),LOGW_E_VNCT,XS_EVNCT_(:,1),XS_EVNCT_(:,2).*PartyE_VNCT,SameE_VNCT,SameE_VNCT.*XS_EVNCT_(:,1),TenureE_VNCT,X_KnotE_VNCT,...
%        LOGW_E_VNCT.^2,(XS_EVNCT_(:,1)).^2,(XS_EVNCT_(:,1)).^2.*SameE_VNCT,PartyE_VNCT.*(XS_EVNCT_(:,2)).^2,TenureE_VNCT.^2,...
%        X_KnotE_VNCT.*(LOGW_E_VNCT*ones(1,8)),X_KnotE_VNCT.*((XS_EVNCT_(:,1).*SameE_VNCT)*ones(1,8)),X_KnotE_VNCT.*((XS_EVNCT_(:,2).*PartyE_VNCT)*ones(1,8)),X_KnotE_VNCT.*((TenureE_VNCT)*ones(1,8)),...
%        PresdumE_VNCT,PresdumE_VNCT.*SameE_VNCT,PresdumE_VNCT.*SameE_VNCT.*XS_EVNCT_(:,1),MidtermE_VNCT,LOGW_E_VNCT.*TenureE_VNCT,LOGW_E_VNCT.^3,X_KnotE_VNCT.*(LOGW_E_VNCT.^2*ones(1,8))];

%Rtotd version... Derivative of spend+fund positive
Xaction=[ones(size(LOGW_E_VNCT,1),1),LOGW_E_VNCT,XS_EVNCT_(:,1),XS_EVNCT_(:,2).*PartyE_VNCT,SameE_VNCT,SameE_VNCT.*XS_EVNCT_(:,1),TenureE_VNCT,X_KnotE_VNCT,...
        LOGW_E_VNCT.^2,(XS_EVNCT_(:,1)).^2,(XS_EVNCT_(:,1)).^2.*SameE_VNCT,PartyE_VNCT.*(XS_EVNCT_(:,2)).^2,TenureE_VNCT.^2,...
        RTotDE_VNCT.*LOGW_E_VNCT,X_KnotE_VNCT.*((XS_EVNCT_(:,1).*SameE_VNCT)*ones(1,8)),X_KnotE_VNCT.*((XS_EVNCT_(:,2).*PartyE_VNCT)*ones(1,8)),X_KnotE_VNCT.*((TenureE_VNCT)*ones(1,8)),...
        PresdumE_VNCT,PresdumE_VNCT.*SameE_VNCT,PresdumE_VNCT.*SameE_VNCT.*XS_EVNCT_(:,1),MidtermE_VNCT,LOGW_E_VNCT.*TenureE_VNCT,RTotDE_VNCT.^2.*LOGW_E_VNCT];%

Xaction2=[ones(size(LOGW_E_VNCT,1),1),LOGW_E_VNCT,XS_EVNCT_(:,1),XS_EVNCT_(:,2).*PartyE_VNCT,SameE_VNCT,SameE_VNCT.*XS_EVNCT_(:,1),TenureE_VNCT,X_KnotE_VNCT,...
        LOGW_E_VNCT.^2,(XS_EVNCT_(:,1)).^2,(XS_EVNCT_(:,1)).^2.*SameE_VNCT,PartyE_VNCT.*(XS_EVNCT_(:,2)).^2,TenureE_VNCT.^2,...
        RTotDE_VNCT.*LOGW_E_VNCT,X_KnotE_VNCT.*((XS_EVNCT_(:,1).*SameE_VNCT)*ones(1,8)),X_KnotE_VNCT.*((XS_EVNCT_(:,2).*PartyE_VNCT)*ones(1,8)),X_KnotE_VNCT.*((TenureE_VNCT)*ones(1,8)),...
        PresdumE_VNCT,PresdumE_VNCT.*SameE_VNCT,PresdumE_VNCT.*SameE_VNCT.*XS_EVNCT_(:,1),MidtermE_VNCT,LOGW_E_VNCT.*TenureE_VNCT,RTotDE_VNCT.^2.*LOGW_E_VNCT];%,RTotDE_VNCT.*LOGW_E_VNCT.^2,RTotDE_VNCT.^2.*LOGW_E_VNCT.^2

Xaction3=[ones(size(LOGW_E_VNCT,1),1),LOGW_E_VNCT,XS_EVNCT_(:,1),XS_EVNCT_(:,2).*PartyE_VNCT,SameE_VNCT,SameE_VNCT.*XS_EVNCT_(:,1),TenureE_VNCT,X_KnotE_VNCT,...
        LOGW_E_VNCT.^2,(XS_EVNCT_(:,1)).^2,(XS_EVNCT_(:,1)).^2.*SameE_VNCT,PartyE_VNCT.*(XS_EVNCT_(:,2)).^2,TenureE_VNCT.^2,...
        RTotDE_VNCT.*LOGW_E_VNCT,X_KnotE_VNCT.*((XS_EVNCT_(:,1).*SameE_VNCT)*ones(1,8)),X_KnotE_VNCT.*((XS_EVNCT_(:,2).*PartyE_VNCT)*ones(1,8)),X_KnotE_VNCT.*((TenureE_VNCT)*ones(1,8)),...
        PresdumE_VNCT,PresdumE_VNCT.*SameE_VNCT,PresdumE_VNCT.*SameE_VNCT.*XS_EVNCT_(:,1),MidtermE_VNCT,RTotDE_VNCT.^2.*LOGW_E_VNCT];

%     Xaction(LOGW_E_VNCT<1,:)=[]; %Drop outliers   
%     Xaction2(LOGW_E_VNCT<1,:)=[]; %Drop outliers
% if size(LOGD_E_VNCT)==size(LOGW_E_VNCT)
% LOGD_E_VNCT(LOGW_E_VNCT<1,:)=[];
% LOGTotal_E_VNCT(LOGW_E_VNCT<1,:)=[];
% LOGW_NXT_E_VNCT(LOGW_E_VNCT<1,:)=[];
% end
coefspend=Xaction\LOGD_E_VNCT;
coeffund=Xaction2\LOGTotal_E_VNCT;
coefsave=Xaction3\LOGW_NXT_E_VNCT;

%If modifying Xaction, make sure the corresponding part on "actions.m" is
%modified as well.

%%
%Specification check for actions
predspend=Xaction*coefspend;
predfund=Xaction2*coeffund;
predsave=Xaction3*coefsave;

errorspend=(sum((LOGD_E_VNCT-predspend).^2))/(var(LOGD_E_VNCT)*length(LOGD_E_VNCT));
errorfund=(sum((LOGTotal_E_VNCT-predfund).^2))/(var(LOGTotal_E_VNCT)*length(LOGTotal_E_VNCT));
errorsave=(sum((LOGW_NXT_E_VNCT-predsave).^2))/(var(LOGW_NXT_E_VNCT)*length(LOGW_NXT_E_VNCT));
%[exp(LOGD_E_VNCT),exp(predspend)];
scatter(LOGD_E_VNCT,predspend);
%%
%Specification check: see if actions are increasing in begcash
derivactions=@(x,tenure,rtotd,coef) [zeros(numel(x),1),ones(numel(x),1),zeros(numel(x),13),2*x,zeros(numel(x),4),ones(numel(x),1)*rtotd,...
    zeros(numel(x),28),tenure*ones(numel(x),1),rtotd.^2*ones(numel(x),1)]*coef;
derivactions2=@(x,tenure,rtotd,coef) [zeros(numel(x),1),ones(numel(x),1),zeros(numel(x),13),2*x,zeros(numel(x),4),ones(numel(x),1)*rtotd,...
    zeros(numel(x),28),tenure*ones(numel(x),1),rtotd.^2*ones(numel(x),1)]*coef;%,rtotd*2*x,rtotd.^2*2*x
derivactions3=@(x,tenure,rtotd,coef) [zeros(numel(x),1),ones(numel(x),1),zeros(numel(x),13),2*x,zeros(numel(x),4),ones(numel(x),1)*rtotd,...
    zeros(numel(x),28),rtotd.^2*ones(numel(x),1)]*coef;
xspace=linspace(5,17,1000).';

%X_Knot version
%scatter(xspace,derivactions(xspace,min(TenureE_VNCT),[max(X_KnotE_VNCT).*[0,1,0,0,0,0,0,0]],coefspend)+derivactions(xspace,min(TenureE_VNCT),[max(X_KnotE_VNCT).*[0,1,0,0,0,0,0,0]],coefsave));

%Rtotd version
figure(1)
scatter(xspace,derivactions(xspace,quantile(TenureE_VNCT,0.95),quantile(RTotDE_VNCT,0.05),coefspend));
figure(2)
scatter(xspace,derivactions2(xspace,quantile(TenureE_VNCT,0.95),quantile(RTotDE_VNCT,0.05),coeffund));
figure(3)
scatter(xspace,derivactions3(xspace,quantile(TenureE_VNCT,0.95),quantile(RTotDE_VNCT,0.05),coefsave));
figure(4)
scatter(xspace,derivactions(xspace,quantile(TenureE_VNCT,0.95),quantile(RTotDE_VNCT,0.05),coefspend)+derivactions3(xspace,quantile(TenureE_VNCT,0.95),quantile(RTotDE_VNCT,0.05),coefsave)-derivactions2(xspace,quantile(TenureE_VNCT,0.95),quantile(RTotDE_VNCT,0.05),coeffund));

%Correlations on data
% corr([LOGW_E_VNCT(TenureE_VNCT>quantile(TenureE_VNCT,0.7)),LOGD_E_VNCT(TenureE_VNCT>quantile(TenureE_VNCT,0.7)),LOGTotal_E_VNCT(TenureE_VNCT>quantile(TenureE_VNCT,0.7)),LOGW_NXT_E_VNCT(TenureE_VNCT>quantile(TenureE_VNCT,0.7))])
% corr([LOGW_E_VNCT(TenureE_VNCT<quantile(TenureE_VNCT,0.3)),LOGD_E_VNCT(TenureE_VNCT<quantile(TenureE_VNCT,0.3)),LOGTotal_E_VNCT(TenureE_VNCT<quantile(TenureE_VNCT,0.3)),LOGW_NXT_E_VNCT(TenureE_VNCT<quantile(TenureE_VNCT,0.3))])
% 
% corr([LOGW_E_VNCT(TenureE_VNCT>quantile(TenureE_VNCT,0.7)&RTotDE_VNCT<quantile(RTotDE_VNCT,0.3)),LOGD_E_VNCT(TenureE_VNCT>quantile(TenureE_VNCT,0.7)&RTotDE_VNCT<quantile(RTotDE_VNCT,0.3)),LOGTotal_E_VNCT(TenureE_VNCT>quantile(TenureE_VNCT,0.7)&RTotDE_VNCT<quantile(RTotDE_VNCT,0.3)),LOGW_NXT_E_VNCT(TenureE_VNCT>quantile(TenureE_VNCT,0.7)&RTotDE_VNCT<quantile(RTotDE_VNCT,0.3))])
% corr([LOGW_E_VNCT(TenureE_VNCT>quantile(TenureE_VNCT,0.7)&RTotDE_VNCT>quantile(RTotDE_VNCT,0.8)),LOGD_E_VNCT(TenureE_VNCT>quantile(TenureE_VNCT,0.7)&RTotDE_VNCT>quantile(RTotDE_VNCT,0.8)),LOGTotal_E_VNCT(TenureE_VNCT>quantile(TenureE_VNCT,0.7)&RTotDE_VNCT>quantile(RTotDE_VNCT,0.8)),LOGW_NXT_E_VNCT(TenureE_VNCT>quantile(TenureE_VNCT,0.7)&RTotDE_VNCT>quantile(RTotDE_VNCT,0.8))])

%%
%Define coefficients
Est3=coefprobwin;
Est51=coefspend;
Est52=coeffund;
Est53=coefsave;



%minthetaall=[Est1(1:10,1);Est2(12,1);Est1(11:20,1);Est2(1:11,1);Est2(13:17,1);Est4(1:135,1);Est5(1:33,1);Est3(1:7,1)];
%minthetaall2=[Est1;Est2;Est3;Est41;Est42;Est43;Est51;Est52;Est53]; %Est1,Est2length differ depending on specification
%mintheta1=[Est1;Est2];
%mintheta2=[Est3;Est41;Est42;Est43;Est51;Est52;Est53];

save Est1 Est1
save Est2 Est2
save Est3 Est3
save Est41 Est41
save Est42 Est42
save Est431 Est431
save Est432 Est432
save Est433 Est433
save Est51 Est51
save Est52 Est52
save Est53 Est53
%save minimizedtheta.txt minthetaall2 -ASCII