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

global LOGW_NXT_E_VCTwnxt
global LOGW_E_VCTwnxt
global SameE_VCTwnxt
global XS_EVCTwnxt_
global PartyE_VCTwnxt

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

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%               q_e               %%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% %thetaS 2 X 1 [unemp,white]  mean of unemp is about 0.07, white is about
% %0.81
% %thetaQ 9 X 1 [ones,d_NC,w_NC,w'_NXT,d_NC^2,w_NC^2,w'_NXT^2,party.*s_NC,tenure]
% %theta  7 X 1 [ones,LOGW_I,XQ,XS.*Party,LOGW_I.^2,XQ.^2,log(T
% %             enure.+1)]  
% %             mean(LOGD_NCE_V) is 12.4 mean(LOGW_NCE_V) is 11.0,
% %             mean(LOGWNXT_NCE_V) is 11.7  mean(LOGD_NXT_E_V) is 12.61
% % [QC1,QC2]   [PRENTAC,PRENTAC.^2]
% % E_VCTa 15 X 1
% % E_VCTt 15 X 1
% % gammaCT 15 X 1
% % thetawin 7 X 1 [ones,LOGW_NXT_E_VCT,XQEVCT,XS_NXT_EVCT.*PartyE_VCT,LOGW_NXT_E_VCT.^2,XQEVCT.^2,...
% %    log(TenureE_VCT.+1)]
% % thetadNC 7 X 1 LOGD_NXT_NCT=[1,XQNCT,XQNCT.^2,LOGW_NXT_NCT,LOGW_NXT_NCT.^2,Party_NCT.*XS_NXT_NCT,log(1+TenureNCT)]*thetadNC...
% % mean(LOGD_NXT_NCT)==12.30 mean(LOGW_NXT_NCT)==11.25
% % thetatNC 7 X 1 LOGTotal_NXT_NCT=[1,XQNCT,XQNCT.^2,LOGW_NXT_NCT,LOGW_NXT_NCT.^2,Party_NCT.*XS_NXT_NCT,log(1+TenureNCT)]*thetatNC
% % mean(LOGTotal_NXT_NCT)=12.56  mean(LOGW_NXT_NCT)=11.2534
% % thetanxtNC 7 X 1 LOGWNXT_NXT_NCT=[1,XQNCT,XQNCT.^2,LOGW_NXT_NCT,LOGW_NXT_NCT.^2,Party_NCT.*XS_NXT_NCT,log(1+TenureNCT)]*thetanxtNC
% % mean(LOGWNXT_NXT_E_VNCT)=11.72, mean(LOGW_NXT_E_VNCT)==11.25
% % win=cdfnorm([1,LOGW_NXT_E_VCT,XQEVCT,XS_NXT_EVCT.*PartyE_VCT,LOGW_NXT_E_VCT.^2,XQEVCT.^2,...
% %    log(TenureE_VCT+1)]*thetawin)
% 
% 
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
% 
% % B_I in [-0.01,0.05]
% B_I=0.01;
% % B_C in [-0.05,0.01]
% B_C=-0.015;
% % Q_C1 in [-0.5,0.5]
% Q_C1=0.06;
% % Q_C2 in [-0.5,0.5]
% Q_C2=-0.01;
% % B_T in [-0.1,0.3]
% B_T=0.04;
% % overall restriction
% % B_I*LOGD_IAC+B_C*LOGD_CAC+XSAC.*PartyAC+XQAC+[PRENTAC,PRENTAC.^2]*[Q_C1,Q
% % _C2]'+B_T*log(TenureAC+1) in [0,1.0]
% 
% % cost1 in [-0.004,0]
% cost1=-0.0013;
% % ben1 in [0,0.004]
% ben1=0.0005;
% % cost2 in [-0.004,0]
% cost2=-0.0015;
% % ben2 in [0,0.004]
% ben2=0.001;
% % sig in [0.001,1]
% sig=0.3;
% 
% % theta0=[thetaS;thetaQ;theta;B_I;B_C;Q_C1;Q_C2;B_T;E_VCTa;E_VCTt;gammaCT;...
% %     thetadNC;thetatNC;thetanxtNC;thetawin;cost1;ben1;sig;cost2;ben2];
% 
% % theta0(1:23,1)=initialvalue(1:23,1);
% % theta0(24:42,1)=initialvalue(24:42,1);  %ACDa
% % theta0(43:61,1)=initialvalue(43:61,1);  %ACDt
% % theta0(62:80,1)=initialvalue(62:80,1);  %gammaACD
% % theta0(81:86,1)=initialvalue(81:86,1); %ACdqea
% % theta0(87:92,1)=initialvalue(87:92,1); %ACdqet
% % theta0(93:98,1)=initialvalue(93:98,1); % gammadqe
% % theta0(99:113,1)=initialvalue(99:113,1); % E_VCTa
% % theta0(114:128,1)=initialvalue(114:128,1); %E_VCTt
% % theta0(129:143,1)=initialvalue(129:143,1); %gammaCT
% % theta0(144:158,1)=initialvalue(144:158,1); % E_VNCTa
% % theta0(159:173,1)=initialvalue(159:173,1); % E_VNCTt
% % theta0(174:188,1)=initialvalue(174:188,1); % gammaNCT
% % theta0(189:194,1)=initialvalue(189:194,1);  % ACwqea
% % theta0(195:200,1)=initialvalue(195:200,1);  % ACwqet
% % theta0(201:206,1)=initialvalue(201:206,1);  % gammawqe
% % theta0(207:213,1)=initialvalue(207:213,1);  %ACvqea
% % theta0(214:220,1)=initialvalue(214:220,1);   %ACvqet
% % theta0(221:227,1)=initialvalue(221:227,1);  %gammavqe
% % theta0(228:234,1)=initialvalue(228:234,1);  % thetawin
% % theta0(235,1)=initialvalue(235,1);  % beta_1
% % theta0(236,1)=initialvalue(236,1);  % beta_2
% % theta0(237,1)=initialvalue(237,1);  % sig
% 
% bddtheta0=zeros(length(theta0),2);
% bddtheta0(1:2,:)=[-3,3;-0.2,0.2]; % thetaS
% bddtheta0(3:11,:)=[0.3,0.9;-0.03,0.03;-0.03,0.03;-0.03,0.03;-0.002,0.002;-0.002,0.002;-0.002,0.002;-0.8,0.8;-0.2,0.3]; % thetaQ
% bddtheta0(12:18,:)=[-1,2;-0.3,0.2;-2,2;-2,2;-0.007,0.007;-2,2;-0.5,0.5]; %theta
% bddtheta0(19,:)=[-0.01,0.05]; %B_I
% bddtheta0(20,:)=[-0.05,0.01]; % B_C
% bddtheta0(21,:)=[-0.5,0.5]; %Q_C1
% bddtheta0(22,:)=[-0.5,0.5]; %Q_C2
% bddtheta0(23,:)=[-0.1,0.3]; % B_T
% bddtheta0(24:38,:)=[(-1)./lbgammaCT,1./lbgammaCT]; % E_VCTa
% bddtheta0(39:53,:)=[lbE_VCTt,ubE_VCTt]; % E_VCTt
% bddtheta0(54:68,:)=[lbgammaCT,ubgammaCT]; % gammaCT
% bddtheta0(69:75,:)=[0,14;-5,5;-3,3;-0.5,1.5;-0.04,0.04;-5,5;-5,5]; % thetadNC
% bddtheta0(76:82,:)=[0,18;-5,5;-3,3;-1.5,0.5;-0.04,0.04;-5,5;-5,5]; % thetatNC
% bddtheta0(83:89,:)=[0,18;-5,5;-3,3;-0.5,1.5;-0.04,0.04;-5,5;-5,5]; % thetanxtNC
% bddtheta0(90:96,:)=[-0.5,2;-0.2,0.2;-2,2;-2,2;-0.01,0.01;-2,2;-0.7,0.7]; % thetawin
% bddtheta0(98,:)=[-0.004,0]; %cost1
% bddtheta0(98,:)=[0,0.004]; %ben1
% bddtheta0(99,:)=[0.001,1]; %sig
% bddtheta0(100,:)=[-0.004,0]; %cost2
% bddtheta0(101,:)=[0,0.004]; %ben2
% Aineq=zeros(6,101);
% bineq=zeros(6,1);
% Aineq(1,1:2)=[mean(Unemployment),mean(White)];
% bineq(1,1)=0.35;
% Aineq(2,:)=-Aineq(1,:);
% bineq(1,1)=0.35;
% 
% Aineq(3,3:11)=[1,mean(LOGD_NC),mean(LOGW_NC),mean(LOGWNXT_NC),mean(LOGD_NC.^2),mean(LOGW_NC.^2)...
%     ,mean(LOGWNXT_NC.^2),0,mean(log(Tenure+1))];
% Aineq(4,:)=-Aineq(3,:);
% bineq(3,1)=0.8;
% bineq(4,1)=-0.2
% Aineq(5,90:96)=[1,mean(LOGW_NXT_E_VCT),0.5,0.1,mean(LOGW_NXT_E_VCT.^2),mean(0.25),mean(log(TenureE_VCT+1))];
% Aineq(6,:)=-Aineq(5,:);
% bineq(5,1)=0.99;
% bineq(6,1)=-0.5;
% 





%%
%%%%%%%%%%%%%%%%%
%Estimation of entry probability and Primary N
%Number=1
%%%%%%%%%%%%%%%%%


Xentry=[ones(Samplesize,1),LOGW_I,Unempsame,Partdemo,log(Tenure+1),X_Knot1];
coefentry=Xentry\Contest;
coefprimaryN=Xentry\Primary_N;

%  Xentry2=[ones(Samplesize,1),LOGW_I,Unempsame,Partdemo,log(Tenure+1),X_Knot1,X_Knot1.*(LOGW_I*ones(1,fineness-1))];
%  coefentry2=Xentry2\Contest;
%  coefprimaryN2=Xentry2\Primary_N;


%%
%%%%%%%%%%%%%%%%%
%Estimation of Vote share equation
%Number=2
%%%%%%%%%%%%%%%%%

Eentry=[ones(SamplesizeAC,1),LOGW_IAC,UnempsameAC,PartdemoAC,log(TenureAC+1),X_KnotAC1]*coefentry;
EprimaryN=[ones(SamplesizeAC,1),LOGW_IAC,UnempsameAC,PartdemoAC,log(TenureAC+1),X_KnotAC1]*coefprimaryN;

Xvoteshare=[LOGD_IAC,LOGD_CAC,UnempsameAC,PartdemoAC,log(TenureAC+1),EprimaryN,Eentry,EprimaryN.*Eentry,RTotD_NCAC,RTotD_NCAC.^2,RTotD_NCAC.^3];
IVvoteshare=[ones(length(VSAC),1),LOGW_IAC,LOGW_IAC.^2,LOGW_IAC.*log(TenureAC+1),UnempsameAC,UnempsqsameAC,PartdemoAC,log(TenureAC+1),(log(TenureAC+1)).^2,X_KnotAC1,X_KnotAC1.*(LOGW_IAC*ones(1,fineness-1))];
IVvar=IVvoteshare.'*IVvoteshare;

coefvoteshare=inv(Xvoteshare.'*IVvoteshare*inv(IVvar)*IVvoteshare.'*Xvoteshare)*(Xvoteshare.'*IVvoteshare*inv(IVvar)*IVvoteshare.'*(VSAC-0.5));


% Eentry2=[ones(SamplesizeAC,1),LOGW_IAC,UnempsameAC,PartdemoAC,log(TenureAC+1),X_KnotAC1,X_KnotAC1.*(LOGW_IAC*ones(1,fineness-1))]*coefentry2;
% EprimaryN2=[ones(SamplesizeAC,1),LOGW_IAC,UnempsameAC,PartdemoAC,log(TenureAC+1),X_KnotAC1,X_KnotAC1.*(LOGW_IAC*ones(1,fineness-1))]*coefprimaryN2;
% Xvoteshare2=[LOGD_IAC,LOGD_CAC,UnempsameAC,PartdemoAC,log(TenureAC+1),EprimaryN2,Eentry2,EprimaryN2.*Eentry2,(EprimaryN2).^2,X_KnotAC1];
% IVvoteshare2=[ones(length(VSAC),1),LOGW_IAC,LOGW_IAC.^2,LOGW_IAC.*log(TenureAC+1),UnempsameAC,UnempsqsameAC,PartdemoAC,log(TenureAC+1),(log(TenureAC+1)).^2,X_KnotAC1,X_KnotAC1.*(LOGW_IAC*ones(1,fineness-1))];
% IVvar2=IVvoteshare2.'*IVvoteshare2;
% 
% coefvoteshare2=inv(Xvoteshare2.'*IVvoteshare2*inv(IVvar2)*IVvoteshare2.'*Xvoteshare2)*(Xvoteshare2.'*IVvoteshare2*inv(IVvar2)*IVvoteshare2.'*(VSAC-0.5));
% % 
%%
%%%%%%%%%%%%%%%%%%
%Estimation of winning probability
%Number=3
%%%%%%%%%%%%%%%%%%
Xprobwin=[ones(length(Win_AC),1),LOGW_IAC,UnempsameAC,UnempsqsameAC,PartdemoAC,LOGW_IAC.^2,X_KnotAC1,log(TenureAC+1),log(TenureAC+1).^2,log(TenureAC+1).*LOGW_IAC,...
    X_KnotAC1.*(LOGW_IAC*ones(1,8)),X_KnotAC1.*(UnempsameAC*ones(1,8)),X_KnotAC1.*(PartdemoAC*ones(1,8)),X_KnotAC1.*(log(TenureAC+1)*ones(1,8))];
%Xprobwin=[ones(length(Win_AC),1),LOGW_IAC,UnempsameAC,PartdemoAC,LOGW_IAC.^2,X_KnotAC1,log(TenureAC+1)];

coefprobwin=Xprobwin\Win_AC;

%If modifying Xprobwin, make sure the corresponding part on "actions.m" is
%modified as well.

%%
%Specification check for winning probability
predprobwin=Xprobwin*coefprobwin;
errorprobwin=sum((Win_AC-Xprobwin*coefprobwin).^2)/(var(Win_AC)*length(Win_AC));
plotprobwin=sortrows([Win_AC,predprobwin]);
scatter(plotprobwin(:,1),plotprobwin(:,2));

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
Sofarbest=10^8;
%bestiter=0;
%iterate=1;
%Delt=0.5;

  OPTIONS=optimset('MaxIter',500000,'MaxFunEvals',1000000,'Display','iter');
load ('./inittheta.mat');
theta0=Est4;

for sss=1:10
    [mintheta,SRR]=fminsearch(@(x) Minimize_Apr2013(x,4),theta0,OPTIONS);
    theta0=mintheta+0.2*(rand(size(mintheta,1),1)-0.5).*mintheta;
    
    if SRR<Sofarbest
        Sofarbest=SRR;
        %    bestiter=iterate;
        %save Est4.txt mintheta SRR -ASCII
        save Est42 mintheta
    end
end

load ('Est42.mat');
Est4=mintheta;
%%
%%%%%%%%%%%%%%%%%%
%Estimation of distribution of spending, saving and fund-raising when
%not contested
%Number=5
%%%%%%%%%%%%%%%%%%


%Each actions
Xaction=[ones(size(LOGD_E_VNCT,1),1),LOGW_E_VNCT,XS_EVNCT_(:,1),XS_EVNCT_(:,2).*PartyE_VNCT,SameE_VNCT,SameE_VNCT.*XS_EVNCT_(:,1),TenureE_VNCT,X_KnotE_VNCT,...
        LOGW_E_VNCT.^2,(XS_EVNCT_(:,1)).^2,(XS_EVNCT_(:,1)).^2.*PartyE_VNCT,SameE_VNCT.*(XS_EVNCT_(:,2)).^2,TenureE_VNCT.^2,...
        X_KnotE_VNCT.*(LOGW_E_VNCT*ones(1,8)),X_KnotE_VNCT.*((XS_EVNCT_(:,1).*SameE_VNCT)*ones(1,8)),X_KnotE_VNCT.*((XS_EVNCT_(:,2).*PartyE_VNCT)*ones(1,8)),X_KnotE_VNCT.*((TenureE_VNCT)*ones(1,8))];

coefspend=Xaction\LOGD_E_VNCT;
coeffund=Xaction\LOGTotal_E_VNCT;
coefsave=Xaction\LOGW_NXT_E_VNCT;

%If modifying Xaction, make sure the corresponding part on "actions.m" is
%modified as well.

%%
%Specification check for actions
predspend=Xaction*coefspend;
predfund=Xaction*coeffund;
predsave=Xaction*coefsave;

errorspend=(sum((LOGD_E_VNCT-predspend).^2))/(var(LOGD_E_VNCT)*length(LOGD_E_VNCT));
errorfund=(sum((LOGTotal_E_VNCT-predfund).^2))/(var(LOGTotal_E_VNCT)*length(LOGTotal_E_VNCT));
errorsave=(sum((LOGW_NXT_E_VNCT-predsave).^2))/(var(LOGW_NXT_E_VNCT)*length(LOGW_NXT_E_VNCT));
[exp(LOGD_E_VNCT),exp(predspend)];
scatter(LOGD_E_VNCT,predspend);

%%
%Define coefficients
Est3=coefprobwin;
Est5=[coefspend;coeffund;coefsave];



%minthetaall=[Est1(1:10,1);Est2(12,1);Est1(11:20,1);Est2(1:11,1);Est2(13:17,1);Est4(1:135,1);Est5(1:33,1);Est3(1:7,1)];
minthetaall2=[Est1;Est2;Est3;Est4;Est5]; %Est1,Est2length differ depending on specification
mintheta1=[Est1;Est2];
mintheta2=[Est3;Est4;Est5];

save Est1 Est1
save Est2 Est2
save Est3 Est3
save Est4 Est4
save Est5 Est5

%save minimizedtheta.txt minthetaall2 -ASCII