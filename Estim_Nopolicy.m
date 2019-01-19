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
global LOGD_E_VCTwnxt
global LOGTotal_E_VCTwnxt
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
OP_INC_IV_july8(6,:)=[];%Drop an Rtotd outlier

 
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
E_V_july8(11:13,:)=[]; %Drop an Rtotd outlier
E_V_july8([850,1110,1111],:)=[]; %Drop an Rtotd outlier


%DO NOT DROP SAMPLES WITH NONPOSITIVE BEGCASH/NONPOSITIVE ENDCASH.
%ESTIMATION RESULTS CONSIDERABLY DIFFERENT.
%%


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
RTotD_NC=(max(0,LOGD_NC).^(-1/2))./LOGTot_NC;  %% ratio of spending over total.
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
%Unemployment=Unemployment-mean(Unemployment);

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


%%
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
%UnemploymentAC=UnemploymentAC-mean(UnemploymentAC);
UnemploymentsqAC=UnemploymentAC.^2;
%Unemployment_NCAC=AC(:,30);
PartisanAC=AC(:,63);


UnempsameAC=UnemploymentAC.*SameAC;
UnempsqsameAC=UnemploymentsqAC.*SameAC;
PartdemoAC=PartisanAC.*PartyAC;
PartsqdemoAC=PartisanAC.^2.*PartyAC;

PresdumAC=AC(:,64);
PresdumsameAC=PresdumAC.*SameAC;
UnemppresdumsameAC=UnemploymentAC.*PresdumsameAC;
MidtermAC=AC(:,65);
%MidtermsameAC=MidtermAC.*SameAC;
%UnempmidtermsameAC=UnemploymentAC.*MidtermsameAC;

%Variables used in number=3, estimation of winning probability
Win_AC=(AC(:,13)>=0.5);  % Dummy for whether incumbent won the election in (t).


%%B-Spline for q_I(RTotD_NCAC), where RTotD_NCAC=(LOGD_NCAC./LOGTot_NCAC)  etc.%%
% take knots to be between 0.88 to 1.025 with 8 knots (8 basis functions). (Almost all lie
% within this range) Let X_Knot be the matrix with 8 columns that contain
% the value of the B-Spline basis function evaluated at each of the 8
% basis functions.


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
%XSEV_(:,2)=XSEV_(:,2)-mean(XSEV_(:,2));
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
LOGD_E_VCT(E_VContestFUL,:)=[];   %log(Spend) limited to contested.
LOGD_E_VCTwnxt=LOGD_E_VCT;
LOGD_E_VCTwnxt(IND5,:)=[];
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
LOGTotal_E_VCTwnxt=LOGTotal_E_VCT;
LOGTotal_E_VCTwnxt(IND5,:)=[];
LOGTotal_E_VNCT=LOGTotal_E_V;
LOGTotal_E_VNCT(E_VNContestFUL,:)=[];
LOGTot_E_VC=log(max(ones(length(E_V_july8),1),E_V_july8(:,4)));  % log(realtotal_c)
LOGTot_E_VC(E_VContestFUL,:)=[];

%%B-Spline for q_I(RTotDE_VCT), where RTotDE_VCT=(NCE_VCT(:,1)./LOGTotal_E_VCT)  etc.%%
% take knots to be between 0.88 to 1.025 with 8 knots (8 basis functions). (Almost all lie
% within this range) Let X_Knot be the matrix with 8 columns that contain
% the value of the B-Spline basis function evaluated at each of the 8
% basis functions.



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
 Xentry=[ones(Samplesize,1),LOGW_I,LOGW_I.^2,Unempsame,Partdemo,log(Tenure+1),Midterm,Presdum,Midterm.*Presdum,X_Knot1,X_Knot1.*(repmat(LOGW_I,1,size(X_Knot1,2)))];
 coefentry=Xentry\Contest;
 residentry=Contest-Xentry*coefentry;
 sigma2entry=residentry.'*residentry/(Samplesize-size(Xentry,2));
 XX=Xentry.'*Xentry;
 covmatentry=inv(XX)*sigma2entry;
 seentry=sqrt(diag(covmatentry));
 
 coefprimaryN=Xentry\Primary_N;
 residprimaryn=Primary_N-Xentry*coefprimaryN;
 sigma2primaryn=residprimaryn.'*residprimaryn/(Samplesize-size(Xentry,2));
 XX=Xentry.'*Xentry;
 covmatprimaryn=inv(XX)*sigma2primaryn;
 seprimaryn=sqrt(diag(covmatprimaryn));
%Specification 2
%  Xentry=[ones(Samplesize,1),LOGW_I,LOGW_I.^2,Unemployment,Unemployment.^2,Same,Partdemo,log(Tenure+1),Midterm,Presdum,Midterm.*Presdum,X_Knot1,X_Knot1.*(repmat(LOGW_I,1,size(X_Knot1,2)))];
%  coefentry=Xentry\Contest;
%  residentry=Contest-Xentry*coefentry;
%  sigma2entry=residentry.'*residentry/(Samplesize-size(Xentry,2));
%  XX=Xentry.'*Xentry;
%  covmatentry=inv(XX)*sigma2entry;
%  seentry=sqrt(diag(covmatentry));
%  
%  coefprimaryN=Xentry\Primary_N;
%  residprimaryn=Primary_N-Xentry*coefprimaryN;
%  sigma2primaryn=residprimaryn.'*residprimaryn/(Samplesize-size(Xentry,2));
%  XX=Xentry.'*Xentry;
%  covmatprimaryn=inv(XX)*sigma2primaryn;
%  seprimaryn=sqrt(diag(covmatprimaryn));
% 
%   space2=min(Unemployment):0.0001:max(Unemployment);
%  regressor=[space2;space2.^2].';%*100;space2.^3*1000;space2.^4*10000;space2.^5*100000]';%;space2.^6*1000000].';%;space2.^5*10000;space2.^(1/2)
%  results=regressor*coefentry(4:5);
% % %results=regressor*coefprimaryN(length(coefprimaryN)-4:length(coefprimaryN));
%  scatter(space2,results)
% 
% space=min(LOGW_I):0.01:max(LOGW_I);
% result=[space;space.^2].'*coefentry(2:3);
% scatter(space,result)
%%%%%%%%%%
%Result of specification check: available on dropbox.

%%
%%%%%%%%%%%%%%%%%
%Estimation of Vote share equation
%Number=2
%%%%%%%%%%%%%%%%%

%Original specification
Eentry=[ones(SamplesizeAC,1),LOGW_IAC,LOGW_IAC.^2,UnempsameAC,PartdemoAC,log(TenureAC+1),MidtermAC,PresdumAC,MidtermAC.*PresdumAC,X_KnotAC1,X_KnotAC1.*(repmat(LOGW_IAC,1,size(X_KnotAC1,2)))]*coefentry;
EprimaryN=[ones(SamplesizeAC,1),LOGW_IAC,LOGW_IAC.^2,UnempsameAC,PartdemoAC,log(TenureAC+1),MidtermAC,PresdumAC,MidtermAC.*PresdumAC,X_KnotAC1,X_KnotAC1.*(repmat(LOGW_IAC,1,size(X_KnotAC1,2)))]*coefprimaryN;

Xvoteshare=[ones(length(VSAC),1),LOGD_IAC,LOGD_CAC,UnemploymentAC,UnemploymentsqAC,UnempsameAC,UnempsqsameAC,PartdemoAC,log(TenureAC+1),EprimaryN,Eentry,EprimaryN.*Eentry,EprimaryN.^2,Eentry.^2,RTotD_NCAC,RTotD_NCAC.^2*100,RTotD_NCAC.^3*1000,RTotD_NCAC.^4*10000,RTotD_NCAC.^5*100000];%,RTotD_NCAC.^6*1000000];%,,RTotD_NCAC.^(1/2)
IVvoteshare=[ones(length(VSAC),1),LOGW_IAC,LOGW_IAC.^2,LOGW_IAC.*log(TenureAC+1),LOGW_IAC.*UnemploymentAC,LOGW_IAC.*PartdemoAC,UnemploymentAC,UnemploymentsqAC,UnempsameAC,UnempsqsameAC,SameAC,PartdemoAC,PartsqdemoAC,log(TenureAC+1),(log(TenureAC+1)).^2,X_KnotAC1,X_KnotAC1.*(LOGW_IAC*ones(1,fineness-1)),PresdumAC,MidtermAC,PresdumAC.*MidtermAC];
IVvar=IVvoteshare.'*IVvoteshare;
IVmix=IVvoteshare*(IVvar\IVvoteshare.');




coefvoteshare=(Xvoteshare.'*IVmix*Xvoteshare)\(Xvoteshare.'*IVmix*(VSAC-0.5));

residual=VSAC-0.5-Xvoteshare*coefvoteshare;
residvar=sum(residual.^2)./(SamplesizeAC-numel(coefvoteshare));
covvoteshare=residvar*inv(Xvoteshare.'*IVmix*Xvoteshare);
sevoteshare=sqrt(diag(covvoteshare));
%Variant 1(to reduce coef of Logd_CAC and increase q_I)
%  Eentry=[ones(SamplesizeAC,1),LOGW_IAC,UnempsameAC,PartdemoAC,log(TenureAC+1),X_KnotAC1]*coefentry;
%  Eentry2=min(1,[ones(SamplesizeAC,1),LOGW_IAC,UnempsameAC,PartdemoAC,log(TenureAC+1),X_KnotAC1]*coefentry);
%  EprimaryN=[ones(SamplesizeAC,1),LOGW_IAC,UnempsameAC,PartdemoAC,log(TenureAC+1),X_KnotAC1]*coefprimaryN;
%  Xvoteshare=[LOGD_IAC,LOGD_CAC,UnempsameAC,PartdemoAC,log(TenureAC+1),EprimaryN,Eentry,EprimaryN.*Eentry,RTotD_NCAC,RTotD_NCAC.^2,RTotD_NCAC.^3];
%  IVvoteshare=[ones(length(VSAC),1),LOGW_IAC,LOGW_IAC.^2,LOGW_IAC.^3,LOGW_IAC.*log(TenureAC+1),UnempsameAC,UnempsqsameAC,PartdemoAC,log(TenureAC+1),(log(TenureAC+1)).^2,X_KnotAC1,X_KnotAC1.*(LOGW_IAC*ones(1,fineness-1)),X_KnotAC1.*(LOGW_IAC.^2*ones(1,fineness-1)),PresdumAC,X_KnotAC1.*(PresdumAC*ones(1,fineness-1)),PresdumsameAC,X_KnotAC1.*(PresdumsameAC*ones(1,fineness-1)),UnemppresdumsameAC,MidtermAC];
%  IVvar=IVvoteshare.'*IVvoteshare;
%  
%  coefvoteshare=inv(Xvoteshare.'*IVvoteshare*inv(IVvar)*IVvoteshare.'*Xvoteshare)*(Xvoteshare.'*IVvoteshare*inv(IVvar)*IVvoteshare.'*(VSAC-0.5));
%  residual=VSAC-0.5-Xvoteshare*coefvoteshare;
%  residstd=sqrt(sum(residual.^2)./(SamplesizeAC-numel(coefvoteshare)-1));

% Eentry2=[ones(SamplesizeAC,1),LOGW_IAC,UnempsameAC,PartdemoAC,log(TenureAC+1),PresdumAC,PresdumsameAC,UnemppresdumsameAC,MidtermAC,X_KnotAC1,X_KnotAC1.*(LOGW_IAC*ones(1,fineness-1))]*coefentry2;
% EprimaryN2=[ones(SamplesizeAC,1),LOGW_IAC,UnempsameAC,PartdemoAC,log(TenureAC+1),PresdumAC,PresdumsameAC,UnemppresdumsameAC,MidtermAC,X_KnotAC1,X_KnotAC1.*(LOGW_IAC*ones(1,fineness-1))]*coefprimaryN2;
% Xvoteshare2=[LOGD_IAC,LOGD_CAC,UnempsameAC,PartdemoAC,log(TenureAC+1),EprimaryN2,Eentry2,EprimaryN2.*Eentry2,(EprimaryN2).^2,X_KnotAC1];
% IVvoteshare2=[ones(length(VSAC),1),LOGW_IAC,LOGW_IAC.^2,LOGW_IAC.*log(TenureAC+1),UnempsameAC,UnempsqsameAC,PartdemoAC,log(TenureAC+1),(log(TenureAC+1)).^2,X_KnotAC1,X_KnotAC1.*(LOGW_IAC*ones(1,fineness-1)),PresdumAC,PresdumsameAC,UnemppresdumsameAC,MidtermAC];
% IVvar2=IVvoteshare2.'*IVvoteshare2;
% 
% coefvoteshare2=inv(Xvoteshare2.'*IVvoteshare2*inv(IVvar2)*IVvoteshare2.'*Xvoteshare2)*(Xvoteshare2.'*IVvoteshare2*inv(IVvar2)*IVvoteshare2.'*(VSAC-0.5));
% %
save residvar residvar
%%
%Specification check
space=(min(UnemploymentAC):0.0001:max(UnemploymentAC));
figure(1)
scatter(space,(coefvoteshare(4)-coefvoteshare(6))*space+(coefvoteshare(5)-coefvoteshare(7))*space.^2)
hold on
scatter(space,(coefvoteshare(4)+coefvoteshare(6))*space+(coefvoteshare(5)+coefvoteshare(7))*space.^2)
figure(2)
histogram(UnemploymentAC)

%%
space2=min(RTotD_NCAC):0.0001:0.07;
regressor=[space2;space2.^2*100;space2.^3*1000;space2.^4*10000;space2.^5*100000]';%;space2.^6*1000000].';%;space2.^5*10000;space2.^(1/2)
results=regressor*coefvoteshare(length(coefvoteshare)-size(regressor,2)+1:length(coefvoteshare));
scatter(space2,results)
%%
%%%%%%%%%%%%%%%%%%
%Estimation of winning probability
%Number=3
%%%%%%%%%%%%%%%%%%
Xprobwin=[ones(length(Win_AC),1),LOGW_IAC,UnemploymentAC,UnemploymentsqAC,PartdemoAC,LOGW_IAC.^2,X_KnotAC1,log(TenureAC+1),log(TenureAC+1).^2,...
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
expectedqc=[EprimaryN,Eentry,EprimaryN.*Eentry,EprimaryN.^2,Eentry.^2]*coefvoteshare(10:14);
histogram(expectedqc)
save expectedqc expectedqc
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

derivprobwinall1=(coefentry(2)).*predprobwin+Eentry.*derivprobwin2-(coefentry(2))*ones(numel(LOGW_IAC),1);
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
save Est1 Est1
save Est2 Est2

%%
%%%%%%%%%%%%%%%%%%
%Estimation of distribution of spending, saving and fund-raising when
%contested
%Number=4
%%%%%%%%%%%%%%%%%%

  
  %%%%%%%%%%
  %4-1: spending
  %%%%%%%%%%

%Find initial theta0
    thetaQ2=Est2(numel(Est2)-4:numel(Est2));
%if numel(thetaQ2)==2
aux=[RTotDE_VCT,RTotDE_VCT.^2*100,RTotDE_VCT.^3*1000,RTotDE_VCT.^4*10000,RTotDE_VCT.^5*100000]*thetaQ2;
    XQEV2=aux-mean(aux); % q_I
    XQEVCT2=XQEV2;
    %XQEVCT2(E_VContestFUL,:)=[]; Already VCT
    XQEVCTwnxt2=XQEVCT2;
    XQEVCTwnxt2(IND5,:)=[];
% else
%     XQEV2=X_KnotEV1*thetaQ2;  % q_I
%     XQEVCT2=XQEV2;
%     XQEVCT2(E_VContestFUL,:)=[];
%     XQEVCTwnxt2=XQEVCT2;
%     XQEVCTwnxt2(IND5,:)=[];
% end
%%

%Start from random points
seed=57;
rng(seed);

nn=500000;
theta0=zeros(24,nn);
Cutoff=[-1;3;7;100];
TAU=find((TenureE_VCT<=Cutoff(3,1))&(TenureE_VCT>Cutoff(2,1)));
% theta0(1:8,:)=20*rand(8,nn);
% theta0(9,:)=0.65*mean(LOGD_E_VCT(TAU))+0.65*mean(LOGD_E_VCT(TAU))*rand(1,nn);
% theta0(10,:)=0.65*mean(XQEVCT2(TAU))+0.65*mean(XQEVCT2(TAU))*rand(1,nn);
% theta0(11,:)=0.65*mean(LOGW_E_VCT(TAU))+0.65*mean(LOGW_E_VCT(TAU))*rand(1,nn);
% theta0(12,:)=0.65*mean(XS_EVCT_(TAU,1))+0.65*mean(XS_EVCT_(TAU,1))*rand(1,nn);
% theta0(13,:)=0.65*mean(XS_EVCT_(TAU,2).*PartyE_VCT(TAU))+0.65*mean(XS_EVCT_(TAU,2).*PartyE_VCT(TAU))*rand(1,nn);
% theta0(14,:)=0.65*mean(SameE_VCT(TAU))+0.65*mean(SameE_VCT(TAU))*rand(1,nn);
% theta0(15,:)=0.65*mean(PresdumE_VCT(TAU))+0.65*mean(PresdumE_VCT(TAU))*rand(1,nn);
% theta0(16,:)=0.65*mean(MidtermE_VCT(TAU))+0.65*mean(MidtermE_VCT(TAU))*rand(1,nn);
% theta0(17:24,:)=rand(8,nn);
theta0(1:8,:)=10*(rand(8,nn)-0.5);
theta0(9:16,:)=100*(rand(8,nn)-0.5);
theta0([1,3],:)=abs(theta0([1,3],:));
theta0(17:24,:)=10*rand(8,nn);

%Saving for young
funvalue11=zeros(1,nn);
funvalue12=zeros(1,nn);
funvalue13=zeros(1,nn);
funvalue21=zeros(1,nn);
funvalue22=zeros(1,nn);
funvalue23=zeros(1,nn);
% funvalue31=zeros(1,nn);
% funvalue32=zeros(1,nn);
% funvalue33=zeros(1,nn);
funvalue311=zeros(1,nn);
funvalue312=zeros(1,nn);
funvalue313=zeros(1,nn);
funvalue321=zeros(1,nn);
funvalue322=zeros(1,nn);
funvalue323=zeros(1,nn);
for i=1:nn
 i
funvalue11(i)=Minimize_Apr2013(theta0(:,i),11);
funvalue12(i)=Minimize_Apr2013(theta0(:,i),12);
funvalue13(i)=Minimize_Apr2013(theta0(:,i),13);
funvalue21(i)=Minimize_Apr2013(theta0(:,i),21);
funvalue22(i)=Minimize_Apr2013(theta0(:,i),22);
funvalue23(i)=Minimize_Apr2013(theta0(:,i),23);
% funvalue31(i)=Minimize_Apr2013(theta0(:,i),31);
% funvalue32(i)=Minimize_Apr2013(theta0(:,i),32);
% funvalue33(i)=Minimize_Apr2013(theta0(:,i),33);
funvalue311(i)=Minimize_Apr2013(theta0(:,i),311);
funvalue312(i)=Minimize_Apr2013(theta0(:,i),312);
funvalue313(i)=Minimize_Apr2013(theta0(:,i),313);
funvalue321(i)=Minimize_Apr2013(theta0(:,i),321);
funvalue322(i)=Minimize_Apr2013(theta0(:,i),322);
funvalue323(i)=Minimize_Apr2013(theta0(:,i),323);
end
id=funvalue11>0;
theta0=theta0(:,id);
funvalue11=funvalue11(id);
funvalue12=funvalue12(id);
funvalue13=funvalue13(id);
funvalue21=funvalue21(id);
funvalue22=funvalue22(id);
funvalue23=funvalue23(id);
% funvalue31=funvalue31(id);
% funvalue32=funvalue32(id);
% funvalue33=funvalue33(id);
funvalue311=funvalue311(id);
funvalue312=funvalue312(id);
funvalue313=funvalue313(id);
funvalue321=funvalue321(id);
funvalue322=funvalue322(id);
funvalue323=funvalue323(id);
save initfun2 funvalue11 funvalue12 funvalue13 funvalue21 funvalue22 funvalue23 funvalue311 funvalue312 funvalue313 funvalue321 funvalue322 funvalue323 theta0 seed
%%
 OPTIONS=optimset('MaxIter',6000,'MaxFunEvals',1000000,'Display','iter');
  
% load initfun2
% %Minimize for 411
% aux=sortrows([funvalue11(funvalue11<min(funvalue11)+1).',theta0(:,funvalue11<min(funvalue11)+1).']);
% theta0mat=aux(:,2:end).';
% theta411=theta0mat(:,1);

load('./Est411.mat');
clear Sofarbest1
Sofarbest1= Minimize_Apr2013(Est411,11);
aux11=Est411;

theta411=aux11;

funvalue=10;
funvalueend=0;
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
save Est411 Est411 Sofarbest1
%%

 OPTIONS=optimset('MaxIter',6000,'MaxFunEvals',1000000,'Display','iter');
  
%load initfun
%Minimize for 412
% aux=sortrows([funvalue12(funvalue12<min(funvalue12)+10).',theta0(:,funvalue12<min(funvalue12)+10).']);
% theta0mat=aux(:,2:end).';
% theta412=theta0mat(:,1);

load('./Est412.mat');
clear Sofarbest2
Sofarbest2= Minimize_Apr2013(Est412,12);
aux12=Est412;

theta412=aux12;

funvalue=10;
funvalueend=0;
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

Est412=aux12;
save Est412 Est412 Sofarbest2

 OPTIONS=optimset('MaxIter',6000,'MaxFunEvals',1000000,'Display','iter');
  
%load initfun
%Minimize for 413
% aux=sortrows([funvalue13(funvalue13<min(funvalue13)+10).',theta0(:,funvalue13<min(funvalue13)+10).']);
% theta0mat=aux(:,2:end).';
% theta413=theta0mat(:,1);

load('./Est413.mat');
clear Sofarbest3
Sofarbest3= Minimize_Apr2013(Est413,13);
aux13=Est413;

theta413=aux13;

funvalue=10;
funvalueend=0;
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

Est413=aux13;
save Est413 Est413 Sofarbest3



  %%%%%%%%%%
  %4-2: fundraising
  %%%%%%%%%%

%fundraising for young
OPTIONS=optimset('MaxIter',6000,'MaxFunEvals',1000000,'Display','iter');
%load initfun
%Minimize for 421
% aux=sortrows([funvalue21(funvalue21<min(funvalue21)+10).',theta0(:,funvalue21<min(funvalue21)+10).']);
% theta0mat=aux(:,2:end).';
% theta421=theta0mat(:,1);

load('./Est421.mat');
clear Sofarbest21
Sofarbest21= Minimize_Apr2013(Est421,21);
aux21=Est421;

theta421=aux21

funvalue=10;
funvalueend=0;
while abs(funvalue-funvalueend)>0.001
    funvalue=Minimize_Apr2013(theta421,21);
    [mintheta21,SRR]=fminsearch(@(x) Minimize_Apr2013(x,21),theta421,OPTIONS);
    funvalueend=Minimize_Apr2013(mintheta21,21);
%      if mod(sss,30)==1
%      theta0=mintheta+0.5*(rand(size(mintheta,1),1)-0.5).*mintheta;
%      theta0([9,17])=mintheta([9,17])+ones(2,1)*0.3;
%      else
        theta421=mintheta21;
%     end
    
    if SRR<Sofarbest21
        Sofarbest21=SRR;
        %    bestiter=iterate;
        %save Est4.txt mintheta -ASCII
        aux21=mintheta21;
    end
end

Est421=aux21;
save Est421 Est421 Sofarbest21




OPTIONS=optimset('MaxIter',6000,'MaxFunEvals',1000000,'Display','iter');
%load initfun
%Minimize for 422
% aux=sortrows([funvalue22(funvalue22<min(funvalue22)+10).',theta0(:,funvalue22<min(funvalue22)+10).']);
% theta0mat=aux(:,2:end).';
% theta422=theta0mat(:,1);

load('./Est422.mat');
clear Sofarbest22
Sofarbest22= Minimize_Apr2013(Est422,22);
aux22=Est422;

theta422=aux22;

funvalue=10;
funvalueend=0;
while abs(funvalue-funvalueend)>0.001
    funvalue=Minimize_Apr2013(theta422,22);
    [mintheta22,SRR]=fminsearch(@(x) Minimize_Apr2013(x,22),theta422,OPTIONS);
    funvalueend=Minimize_Apr2013(mintheta22,22);
%      if mod(sss,30)==1
%      theta0=mintheta+0.5*(rand(size(mintheta,1),1)-0.5).*mintheta;
%      theta0([9,17])=mintheta([9,17])+ones(2,1)*0.3;
%      else
        theta422=mintheta22;
%     end
    
    if SRR<Sofarbest22
        Sofarbest22=SRR;
        %    bestiter=iterate;
        %save Est4.txt mintheta -ASCII
        aux22=mintheta22;
    end
end

Est422=aux22;
save Est422 Est422 Sofarbest22

OPTIONS=optimset('MaxIter',6000,'MaxFunEvals',1000000,'Display','iter');
%load initfun
%Minimize for 423
% aux=sortrows([funvalue23(funvalue23<min(funvalue23)+10).',theta0(:,funvalue23<min(funvalue23)+10).']);
% theta0mat=aux(:,2:end).';
% theta423=theta0mat(:,1);

load('./Est423.mat');
clear Sofarbest23
Sofarbest23= Minimize_Apr2013(Est423,23);
aux23=Est423;

theta423=aux23;

funvalue=10;
funvalueend=0;
while abs(funvalue-funvalueend)>0.001
    funvalue=Minimize_Apr2013(theta423,23);
    [mintheta23,SRR]=fminsearch(@(x) Minimize_Apr2013(x,23),theta423,OPTIONS);
    funvalueend=Minimize_Apr2013(mintheta23,23);
%      if mod(sss,30)==1
%      theta0=mintheta+0.5*(rand(size(mintheta,1),1)-0.5).*mintheta;
%      theta0([9,17])=mintheta([9,17])+ones(2,1)*0.3;
%      else
        theta423=mintheta23;
%     end
    
    if SRR<Sofarbest23
        Sofarbest23=SRR;
        %    bestiter=iterate;
        %save Est4.txt mintheta -ASCII
        aux23=mintheta23;
    end
end

Est423=aux23;
save Est423 Est423 Sofarbest23


  


  %%%%%%%%%%
  %4-3: saving
  %DIRECT ESTIMATION
  %%%%%%%%%%

%Saving for young
% OPTIONS=optimset('MaxIter',6000,'MaxFunEvals',1000000,'Display','iter');
% load initfun
% %Minimize for 431
% aux=sortrows([funvalue31(funvalue31<min(funvalue31)+10).',theta0(:,funvalue31<min(funvalue31)+10).']);
% theta0mat=aux(:,2:end).';
% theta431=theta0mat(:,1);
% 
% load('./Est431.mat');
% clear Sofarbest31
% Sofarbest31= Minimize_Apr2013(Est431,31);
% aux31=Est431;
% 
% funvalue=10;
% funvalueend=0;
% while abs(funvalue-funvalueend)>0.001
%     funvalue=Minimize_Apr2013(theta431,31);
%     [mintheta31,SRR]=fminsearch(@(x) Minimize_Apr2013(x,31),theta431,OPTIONS);
%     funvalueend=Minimize_Apr2013(mintheta31,31);
% %      if mod(sss,30)==1
% %      theta0=mintheta+0.5*(rand(size(mintheta,1),1)-0.5).*mintheta;
% %      theta0([9,17])=mintheta([9,17])+ones(2,1)*0.3;
% %      else
%         theta431=mintheta31;
% %     end
%     
%     if SRR<Sofarbest31
%         Sofarbest31=SRR;
%         %    bestiter=iterate;
%         %save Est4.txt mintheta -ASCII
%         aux31=mintheta31;
%     end
% end
% 
% 
% Est431=aux31;
% save Est431 Est431 Sofarbest31
% %%
% 
% OPTIONS=optimset('MaxIter',6000,'MaxFunEvals',1000000,'Display','iter');
% load initfun
% %Minimize for 432
% aux=sortrows([funvalue32(funvalue32<min(funvalue32)+10).',theta0(:,funvalue32<min(funvalue32)+10).']);
% theta0mat=aux(:,2:end).';
% theta432=theta0mat(:,1);
% 
% load('./Est432.mat');
% clear Sofarbest32
% Sofarbest32= Minimize_Apr2013(Est432,32);
% aux32=Est432;
% 
% funvalue=10;
% funvalueend=0;
% while abs(funvalue-funvalueend)>0.001
%     funvalue=Minimize_Apr2013(theta432,32);
%     [mintheta32,SRR]=fminsearch(@(x) Minimize_Apr2013(x,32),theta432,OPTIONS);
%     funvalueend=Minimize_Apr2013(mintheta32,32);
% %      if mod(sss,30)==1
% %      theta0=mintheta+0.5*(rand(size(mintheta,1),1)-0.5).*mintheta;
% %      theta0([9,17])=mintheta([9,17])+ones(2,1)*0.3;
% %      else
%         theta432=mintheta32;
% %     end
%     
%     if SRR<Sofarbest32
%         Sofarbest32=SRR;
%         %    bestiter=iterate;
%         %save Est4.txt mintheta -ASCII
%         aux32=mintheta32;
%     end
% end
% 
% 
% Est432=aux32;
% save Est432 Est432 Sofarbest32
% %%
% 
% OPTIONS=optimset('MaxIter',6000,'MaxFunEvals',1000000,'Display','iter');
% load initfun
% %Minimize for 433
% aux=sortrows([funvalue33(funvalue33<min(funvalue33)+10).',theta0(:,funvalue33<min(funvalue33)+10).']);
% theta0mat=aux(:,2:end).';
% theta433=theta0mat(:,1);
% 
% load('./Est433.mat');
% clear Sofarbest33
% Sofarbest33= Minimize_Apr2013(Est433,33);
% aux33=Est433;
% 
% funvalue=10;
% funvalueend=0;
% while abs(funvalue-funvalueend)>0.001
%     funvalue=Minimize_Apr2013(theta433,33);
%     [mintheta33,SRR]=fminsearch(@(x) Minimize_Apr2013(x,33),theta433,OPTIONS);
%     funvalueend=Minimize_Apr2013(mintheta33,33);
% %      if mod(sss,30)==1
% %      theta0=mintheta+0.5*(rand(size(mintheta,1),1)-0.5).*mintheta;
% %      theta0([9,17])=mintheta([9,17])+ones(2,1)*0.3;
% %      else
%         theta433=mintheta33;
% %     end
%     
%     if SRR<Sofarbest33
%         Sofarbest33=SRR;
%         %    bestiter=iterate;
%         %save Est4.txt mintheta -ASCII
%         aux33=mintheta33;
%     end
% end
% 
% Est433=aux33;
% save Est433 Est433 Sofarbest33


%%%%%%%%%
%Alternative estimating strategy for saving:
%Estimate both spending and fundraising conditional on winnning, and use
%the budget constraint to recover saving.
%%%%%%%%%

OPTIONS=optimset('MaxIter',6000,'MaxFunEvals',1000000,'Display','iter');
%load initfun
%Minimize for 4311
% aux=sortrows([funvalue311(funvalue311<min(funvalue311)+10).',theta0(:,funvalue311<min(funvalue311)+10).']);
% theta0mat=aux(:,2:end).';
% theta4311=theta0mat(:,1);

load('./Est4311.mat');
clear Sofarbest311
Sofarbest311= Minimize_Apr2013(Est4311,311);
aux311=Est4311;

theta4311=aux311;

funvalue=10;
funvalueend=0;
while abs(funvalue-funvalueend)>0.001
    funvalue=Minimize_Apr2013(theta4311,311);
    [mintheta311,SRR]=fminsearch(@(x) Minimize_Apr2013(x,311),theta4311,OPTIONS);
    funvalueend=Minimize_Apr2013(mintheta311,311);
%      if mod(sss,30)==1
%      theta0=mintheta+0.5*(rand(size(mintheta,1),1)-0.5).*mintheta;
%      theta0([9,17])=mintheta([9,17])+ones(2,1)*0.3;
%      else
        theta4311=mintheta311;
%     end
    
    if SRR<Sofarbest311
        Sofarbest311=SRR;
        %    bestiter=iterate;
        %save Est4.txt mintheta -ASCII
        aux311=mintheta311;
    end
end


Est4311=aux311;
save Est4311 Est4311 Sofarbest311



OPTIONS=optimset('MaxIter',6000,'MaxFunEvals',1000000,'Display','iter');
%load initfun
%Minimize for 4312
% aux=sortrows([funvalue312(funvalue312<min(funvalue312)+10).',theta0(:,funvalue312<min(funvalue312)+10).']);
% theta0mat=aux(:,2:end).';
% theta4312=theta0mat(:,1);

load('./Est4312.mat');
clear Sofarbest312
Sofarbest312= Minimize_Apr2013(Est4312,312);
aux312=Est4312;

theta4312=aux312;

funvalue=10;
funvalueend=0;
while abs(funvalue-funvalueend)>0.001
    funvalue=Minimize_Apr2013(theta4312,312);
    [mintheta312,SRR]=fminsearch(@(x) Minimize_Apr2013(x,312),theta4312,OPTIONS);
    funvalueend=Minimize_Apr2013(mintheta312,312);
%      if mod(sss,30)==1
%      theta0=mintheta+0.5*(rand(size(mintheta,1),1)-0.5).*mintheta;
%      theta0([9,17])=mintheta([9,17])+ones(2,1)*0.3;
%      else
        theta4312=mintheta312;
%     end
    
    if SRR<Sofarbest312
        Sofarbest312=SRR;
        %    bestiter=iterate;
        %save Est4.txt mintheta -ASCII
        aux312=mintheta312;
    end
end


Est4312=aux312;
save Est4312 Est4312 Sofarbest312




OPTIONS=optimset('MaxIter',6000,'MaxFunEvals',1000000,'Display','iter');
%load initfun
%Minimize for 4313
% aux=sortrows([funvalue313(funvalue313<min(funvalue313)+10).',theta0(:,funvalue313<min(funvalue313)+10).']);
% theta0mat=aux(:,2:end).';
% theta4313=theta0mat(:,1);

load('./Est4313.mat');
clear Sofarbest313
Sofarbest313= Minimize_Apr2013(Est4313,313);
aux313=Est4313;

theta4313=aux313;

funvalue=10;
funvalueend=0;
while abs(funvalue-funvalueend)>0.001
    funvalue=Minimize_Apr2013(theta4313,313);
    [mintheta313,SRR]=fminsearch(@(x) Minimize_Apr2013(x,313),theta4313,OPTIONS);
    funvalueend=Minimize_Apr2013(mintheta313,313);
%      if mod(sss,30)==1
%      theta0=mintheta+0.5*(rand(size(mintheta,1),1)-0.5).*mintheta;
%      theta0([9,17])=mintheta([9,17])+ones(2,1)*0.3;
%      else
        theta4313=mintheta313;
%    end
    
    if SRR<Sofarbest313
        Sofarbest313=SRR;
%            bestiter=iterate;
%         save Est4.txt mintheta -ASCII
        aux313=mintheta313;
    end
end

Est4313=aux313;
save Est4313 Est4313 Sofarbest313


%ESTIMATE SPENDING CONDITIONAL ON WINNING
OPTIONS=optimset('MaxIter',6000,'MaxFunEvals',1000000,'Display','iter');
%load initfun
%Minimize for 4321
% aux=sortrows([funvalue321(funvalue321<min(funvalue321)+10).',theta0(:,funvalue321<min(funvalue321)+10).']);
% theta0mat=aux(:,2:end).';
% theta4321=theta0mat(:,1);

load('./Est4321.mat');
clear Sofarbest321
Sofarbest321= Minimize_Apr2013(Est4321,321);
aux321=Est4321;

theta4321=aux321;

funvalue=10;
funvalueend=0;
while abs(funvalue-funvalueend)>0.001
    funvalue=Minimize_Apr2013(theta4321,321);
    [mintheta321,SRR]=fminsearch(@(x) Minimize_Apr2013(x,321),theta4321,OPTIONS);
    funvalueend=Minimize_Apr2013(mintheta321,321);
%      if mod(sss,30)==1
%      theta0=mintheta+0.5*(rand(size(mintheta,1),1)-0.5).*mintheta;
%      theta0([9,17])=mintheta([9,17])+ones(2,1)*0.3;
%      else
        theta4321=mintheta321;
%     end
    
    if SRR<Sofarbest321
        Sofarbest321=SRR;
        %    bestiter=iterate;
        %save Est4.txt mintheta -ASCII
        aux321=mintheta321;
    end
end


Est4321=aux321;
save Est4321 Est4321 Sofarbest321


OPTIONS=optimset('MaxIter',6000,'MaxFunEvals',1000000,'Display','iter');
%load initfun
%Minimize for 4322
% aux=sortrows([funvalue322(funvalue322<min(funvalue322)+10).',theta0(:,funvalue322<min(funvalue322)+10).']);
% theta0mat=aux(:,2:end).';
% theta4322=theta0mat(:,1);

load('./Est4322.mat');
clear Sofarbest322
Sofarbest322= Minimize_Apr2013(Est4322,322);
aux322=Est4322;

theta4322=aux322;

funvalue=10;
funvalueend=0;
while abs(funvalue-funvalueend)>0.001
    funvalue=Minimize_Apr2013(theta4322,322);
    [mintheta322,SRR]=fminsearch(@(x) Minimize_Apr2013(x,322),theta4322,OPTIONS);
    funvalueend=Minimize_Apr2013(mintheta322,322);
%      if mod(sss,30)==1
%      theta0=mintheta+0.5*(rand(size(mintheta,1),1)-0.5).*mintheta;
%      theta0([9,17])=mintheta([9,17])+ones(2,1)*0.3;
%      else
        theta4322=mintheta322;
%     end
    
    if SRR<Sofarbest322
        Sofarbest322=SRR;
        %    bestiter=iterate;
        %save Est4.txt mintheta -ASCII
        aux322=mintheta322;
    end
end


Est4322=aux322;
save Est4322 Est4322 Sofarbest322

OPTIONS=optimset('MaxIter',6000,'MaxFunEvals',1000000,'Display','iter');
%load initfun
%Minimize for 4323
% aux=sortrows([funvalue323(funvalue323<min(funvalue323)+10).',theta0(:,funvalue323<min(funvalue323)+10).']);
% theta0mat=aux(:,2:end).';
% theta4323=theta0mat(:,1);

load('./Est4323.mat');
clear Sofarbest323
Sofarbest323= Minimize_Apr2013(Est4323,323);
aux323=Est4323;

theta4323=aux323;

funvalue=10;
funvalueend=0;
while abs(funvalue-funvalueend)>0.001
    funvalue=Minimize_Apr2013(theta4323,323);
    [mintheta323,SRR]=fminsearch(@(x) Minimize_Apr2013(x,323),theta4323,OPTIONS);
    funvalueend=Minimize_Apr2013(mintheta323,323);
%      if mod(sss,30)==1
%      theta0=mintheta+0.5*(rand(size(mintheta,1),1)-0.5).*mintheta;
%      theta0([9,17])=mintheta([9,17])+ones(2,1)*0.3;
%      else
        theta4323=mintheta323;
%     end
    
    if SRR<Sofarbest323
        Sofarbest323=SRR;
        %    bestiter=iterate;
        %save Est4.txt mintheta -ASCII
        aux323=mintheta323;
    end
end

Est4323=aux323;
save Est4323 Est4323 Sofarbest323


%%
%when not computing number=4
load('Est411.mat');
load('Est412.mat');
load('Est413.mat');
load('Est421.mat');
load('Est422.mat');
load('Est423.mat');


%Specification check: See if actions have right derivative
spendderiv=[Est411(1),Est411(3);Est412(1),Est412(3);Est413(1),Est413(3)];
fundderiv=[Est421(1),Est421(3);Est422(1),Est422(3);Est423(1),Est423(3)];
%savederiv=[Est431(1),Est431(3);Est432(1),Est432(3);Est433(1),Est433(3)];

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



save Est1 Est1
save Est2 Est2
save Est3 Est3
save Est411 Est411
save Est412 Est412
save Est413 Est413
save Est421 Est421
save Est422 Est422
save Est423 Est423
save Est4311 Est4311
save Est4312 Est4312
save Est4313 Est4313
save Est4321 Est4321
save Est4322 Est4322
save Est4323 Est4323
save Est51 Est51
save Est52 Est52
save Est53 Est53

%%
%Vote share equation
unemp=0.05:0.001:0.1;


