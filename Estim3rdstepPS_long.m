% clear all
clear
% **
% make sure you get the sign of q_e right. Also note that f_d_qe, f_w_qe have the PARTY of Challenger.
% **

global iterate
global results
global LOGW_I
global Contest
global Party
global Delt
global Entry
global Winrnd
global EntryC
global WinrndC
global EntryR
global WinrndR
global Betawh
global Betaump
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
global LOGTot_NCE_V
global LOGW_E_VCT
global LOGW_E_VNCT
global LOGW_E_VCTwnxt
global LOGD_E_VCT
global LOGD_E_VNCT
global N
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
global IND6NCT
global C
global DC
global Cc
global DCc
global RTotDE_V
global RTotDE_VCT
global XQEV
global XQEV2
global XQEVCT
global Pen
global RTotDE_VNCT
global VSEVCT;
global XXX
global PrE3step
global q_c3step
global BX13step
global LOGD_E_VC3step
global LOGLOGD_E_V
global LOGLOGTot_E_V
global LOGTot_E_VC3step
global rtotd3step
global max_N
global qEntry
global LOGW_E_VCT3step
global qc_LB
global qc_UB
global LOGW_NXT3step
global XS3step_
global Tenure3step
global Party3step
global q_i3step
global TenureE_VCT3step
global interest
global coef
load('./op_inc_july8.mat')

load('./op_inc_iv_july8.mat')

load ('./State_Transition.mat')
load ('./retire.txt')
load ('./E_V_july8.mat')
load ('./Fststage_kettei.txt')
load ('./Sndstage_kettei_long.txt')
load ('./Trdstage_initial.txt')
load ('./q_C_E_VCT.mat') %% quality of the challengers
load ('./coef.mat')
load ('./BX1.mat') 
load ('./DEL_Pen.mat')
rand('state',2000);
randn('state',20);

dele1=find(sum(isnan(OP_INC_IV_july8),2)>0);
OP_INC_IV_july8(dele1,:)=[];
deleA=find(OP_INC_IV_july8(:,16)==0);
OP_INC_IV_july8(deleA,:)=[];
dele2=find(sum(isnan(OP_INC_july8),2)>0);
OP_INC_july8(dele2,:)=[];
deleB=find(((OP_INC_july8(:,16)==0).*OP_INC_july8(:,8))==1);
OP_INC_july8(deleB,:)=[];
dele3=find(sum(isnan(E_V_july8),2)>0);
E_V_july8(dele3,:)=[];
deleC=find((E_V_july8(:,16)==0).*E_V_july8(:,8)==1);
E_V_july8(deleC,:)=[];

Sofarbest=10^8;
bestiter=0;
Delt=0.5;
interest=0.1;
dfSim=100;              %dfSim is the number of simulations used to evaluate the conditional expectation of the entrant's quality.
N=20;                      % N is the number of simulations.
NumSim=50;              %NumSim is the number of simulations in the 1st step to obtain f(.,q_e)
max_N=max(OP_INC_july8(:,36));               % max_N is the maximum number of possible entrants in the Primary.
T=10;                      % T is the number of periods that we move the simlation forward.
qc_LB=-0.4;
qc_UB=0.4;
Entry=rand(T,N,length(E_V_july8));             %Simulation draw for computing the continuation value E_V
Winrnd=rand(T,N,length(E_V_july8));
EntryC=rand(T,N,length(E_V_july8));             %Simulation draw for computing the continuation value E_V (for challenger FOC)
WinrndC=rand(T,N,length(E_V_july8));
EntryR=rand(T,N,length(E_V_july8));             %Simulation draw for computing the continuation value E_V (for Reserv)
WinrndR=rand(T,N,length(E_V_july8));

Retirernd=(rand(T-1,N,length(E_V_july8))<repmat(retire',[1 N, length(E_V_july8)]));
Retirernd(T,:,:)=1;
RetirerndC=(rand(T-1,N,length(E_V_july8))<repmat(retire',[1 N, length(E_V_july8)]));
RetirerndC(T,:,:)=1;
RetirerndR=(rand(T-1,N,length(E_V_july8))<repmat(retire',[1 N, length(E_V_july8)]));
RetirerndR(T,:,:)=1;




for i=1:N
    for j=1:length(E_V_july8)
        Ret(i,j)=min(find(Retirernd(:,i,j)==1));
    end
end
for i=1:N
    for j=1:length(E_V_july8)
        RetC(i,j)=min(find(RetirerndC(:,i,j)==1));
    end
end
for i=1:N
    for j=1:length(E_V_july8)
        RetR(i,j)=min(find(RetirerndR(:,i,j)==1));
    end
end


minmax=min(max(Ret,[],1))

minmax=min(max(RetC,[],1))

minmax=min(max(RetR,[],1))


dF_gamma_ct=rand(T,N,length(E_V_july8));
dF_total_ct=rand(T,N,length(E_V_july8));
dF_nxt_nxt_ct=rand(T,N,length(E_V_july8));


dF_gamma_ctC=rand(T,N,length(E_V_july8));
dF_total_ctC=rand(T,N,length(E_V_july8));
dF_nxt_nxt_ctC=rand(T,N,length(E_V_july8));


dF_gamma_ctR=rand(T,N,length(E_V_july8));
dF_total_ctR=rand(T,N,length(E_V_july8));
dF_nxt_nxt_ctR=rand(T,N,length(E_V_july8));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Preliminary Step: Transition Estimation  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Regress unemployment Jump and the pct of white population on a constant
%% and one period lagged value %%
%% This can be done separately from the other first step estimation as
%% there is not interdependence between the Transition Estimation and the
%% other first-step moment restrictions.

Ywh=[];
Xw=[];
Yump=[];
Xu=[];

%% Create Ywh, Xw (Yump, Xu) %%
%% For each district, stack the percentage of whites into ever longer
%% column vector. Do the same for unemployment >> Ywh, Yump
%% We omit the first observation, because for that value, we don't observe
%% the 1 period lagged value.

%% Also create Xw and Xu analogously, but including all periods but the
%% last.
%% useornot=0 if we do not use the particular CD for any of the rest of the
%% estimation: due to the fact that we cannot invert any of the policy
%% functions.

unuse=(useornot==0);
State_Transition(unuse,:)=[];
pctblk(unuse,:)=[];
pctwh(unuse,:)=[];
pctothr(unuse,:)=[];
unemployment(unuse,:)=[];

for i=1:size(State_Transition,1)
    Jwh=pctwh(i,:);
    Jump=unemployment(i,:);
    Kwh=(Jwh==-1);
    Kump=(Jump==-1);
    Jwh(Kwh)=[];
    Jump(Kump)=[];
    Jwh=Jwh';
    Jump=Jump';
    Xwh(:,1)=ones(size(Jwh,1)-1,1);   %% Constant
    Xwh(:,2)=Jwh(1:size(Jwh,1)-1,1);      %% lag value.
    Xump(:,1)=ones(size(Jump,1)-1,1);
    Xump(:,2)=Jump(1:size(Jump,1)-1,1);
    Ywh=[Ywh;Jwh(2:size(Jwh,1),1)];
    Xw=[Xw;Xwh];
    Yump=[Yump;Jump(2:size(Jump,1),1)];
    Xu=[Xu;Xump];
    Xwh=[];
    Xump=[];
    Jwh=[];
    Jump=[];
end

Betawh=inv(Xw'*Xw)*Xw'*Ywh;
Betaump=inv(Xu'*Xu)*Xu'*Yump;

epswh=(1/size(Ywh,1))*sum((Ywh-Xw*Betawh).^2);
epsump=(1/size(Yump,1))*sum((Yump-Xu*Betaump).^2);
Shockump=sqrt(epsump)*randn(T,N,size(E_V_july8,1));
Shockwh=sqrt(epswh)*randn(T,N,size(E_V_july8,1));
epswh2=(1/size(Ywh,1))*sum((Ywh-Xw*Betawh).^2);
epsump2=(1/size(Yump,1))*sum((Yump-Xu*Betaump).^2);
Shockump2=sqrt(epsump)*randn(T,N,size(E_V_july8,1));
Shockwh2=sqrt(epswh)*randn(T,N,size(E_V_july8,1));
epswh3=(1/size(Ywh,1))*sum((Ywh-Xw*Betawh).^2);
epsump3=(1/size(Yump,1))*sum((Yump-Xu*Betaump).^2);
Shockump3=sqrt(epsump)*randn(T,N,size(E_V_july8,1));
Shockwh3=sqrt(epswh)*randn(T,N,size(E_V_july8,1));

Estimation_OP=OP_INC_july8;
Samplesize=size(Estimation_OP,1);
Year=Estimation_OP(:,2); % year of election.

Contest=Estimation_OP(:,8); % contest
LOGD_I=log(max(ones(Samplesize,1),Estimation_OP(:,10))); % spend
LOGW_I=log(max(ones(Samplesize,1),Estimation_OP(:,11))); % w_I
LOGW_INXT=log(max(ones(Samplesize,1),Estimation_OP(:,12))); % W_I'
VS=Estimation_OP(:,13); %V_I
%IDN=Estimation_OP(:,21); %Candidate ID


LOGTot_NC=log(max(ones(Samplesize,1),Estimation_OP(:,31))); % Tot_INC
LOGD_NC=log(max(ones(Samplesize,1),Estimation_OP(:,32))); % D_INC
LOGW_NC=log(max(ones(Samplesize,1),Estimation_OP(:,33))); % W_INC
LOGWNXT_NC=log(max(ones(Samplesize,1),Estimation_OP(:,34))); %W'_INC
RTotD_NC=(max(0,LOGD_NC).^(-1/2))./LOGTot_NC;  %% ratio of spending over total. (modify)
RTotD_NC=RTotD_NC.*(exp(LOGTot_NC)./exp(LOGD_NC));
% LOGLOGD_NC=log(max(ones(Samplesize,1),LOGD_NC));
% LOGLOGTot_NC=log(max(ones(Samplesize,1),LOGTot_NC));
Tenure_NC=Estimation_OP(:,35); % Tenure_INC
Tenure=Estimation_OP(:,15);
Party=3*ones(Samplesize,1)-2*Estimation_OP(:,3);  %%Party=1  if candidate i is Democrat and Party=-1 if candidate i is Republican.

White=Estimation_OP(:,21);
Black=Estimation_OP(:,22);
Other=Estimation_OP(:,23);
White_NC=Estimation_OP(:,27);
Black_NC=Estimation_OP(:,28);
Other_NC=Estimation_OP(:,29);
Unemployment=Estimation_OP(:,24);
Unemployment_NC=Estimation_OP(:,30);
XS_=[Unemployment,White];
%%B-Spline for q_I(RTotD_NC), where RTotD_NC=(LOGD_NC./LOGTot_NC)  etc.%%
%% take knots to be between 0.88 to 1.025 with 8 knots (8 basis functions). (Almost all lie
%% within this range) Let X_Knot be the matrix with 8 columns that contain
%% the value of the B-Spline basis function evaluated at each of the 8
%% basis functions.
% mesh=quantile(LOGLOGD_NC,[.125;.25;.375;.5;.625;.75;.875]);
% mesh=[min(LOGLOGD_NC);mesh;max(LOGLOGD_NC)];
% X_Knot1=(LOGLOGD_NC<mesh(2,1)).*(1-(LOGLOGD_NC-mesh(1,1))/(mesh(2,1)-mesh(1,1)));
% for i=1:6
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
% for i=1:6
%     PLUS=(LOGLOGTot_NC>=mesh(i+1,1)).*(LOGLOGTot_NC<mesh(i+2,1)).*((LOGLOGTot_NC-mesh(i+1,1))/(mesh(i+2,1)-mesh(i+1,1)))...
%         +(LOGLOGTot_NC>mesh(i+2,1)).*(LOGLOGTot_NC<mesh(i+3,1)).*(1-(LOGLOGTot_NC-mesh(i+2,1))/(mesh(i+3,1)-mesh(i+2,1)));
%     X_Knot2=[X_Knot2,PLUS];
% end
% PLUS=(LOGLOGTot_NC>=mesh(8,1)).*(LOGLOGTot_NC-mesh(8,1))/(mesh(9,1)-mesh(8,1));
% X_Knot2=[X_Knot2,PLUS];


mesh=quantile(RTotD_NC,[.125;.25;.375;.5;.625;.75;.875]);
mesh=[min(RTotD_NC);mesh;max(RTotD_NC)];
X_Knot1=(RTotD_NC<mesh(2,1)).*(1-(RTotD_NC-mesh(1,1))/(mesh(2,1)-mesh(1,1)));
for i=0:6
    PLUS=(RTotD_NC>=mesh(i+1,1)).*(RTotD_NC<mesh(i+2,1)).*((RTotD_NC-mesh(i+1,1))/(mesh(i+2,1)-mesh(i+1,1)))...
        +(RTotD_NC>=mesh(i+2,1)).*(RTotD_NC<mesh(i+3,1)).*(1-(RTotD_NC-mesh(i+2,1))/(mesh(i+3,1)-mesh(i+2,1)));
    X_Knot1=[X_Knot1,PLUS];
end
PLUS=(RTotD_NC>=mesh(8,1)).*(RTotD_NC-mesh(8,1))/(mesh(9,1)-mesh(8,1));
X_Knot1=[X_Knot1,PLUS];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Estimation of Probability     %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%S=(c0_s+c1_s*pctwhite*c2_s*unemployment)
%q_I=(c1+c2*d_nc+c3*w_nc+c4*w'_nc+c5*S+c6*(d_nc)^2+
%c7*(w_nc)^2+c8*(w'_nc)^2+c9*PartyDem+c10*PartyRep+c11*S*PartyDem+c12*S*PartyRep

%[[+c9_q*(w_i_nc)(w'_i_nc)+ c10_q*(w_i_nc)(d_i_nc)+c11_q*(w_i_nc)S+c12_q*(w'_i_nc)(d_i_nc)+
%c13_q*(w'_i_nc)S+c14_q*(d_i_nc)S]]

thetaProb0=zeros(9,1);
% OPTIONS=optimset('TolX',1e-5,'TolFun',1e-4,'Diagnostics','off','MaxIter',500,'MaxFunEvals',1000);
% [thetaEnt,fval]=fminsearch(@(theta) Entryprob(theta,XSQ),thetaProb0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%      Ai & Chen     %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AC=OP_INC_IV_july8;
SamplesizeAC=size(OP_INC_IV_july8,1);
EPS=randn(SamplesizeAC,NumSim);
EPS2=randn(SamplesizeAC,NumSim);
EPS3=randn(SamplesizeAC,NumSim);
LOGD_IAC=log(max(ones(SamplesizeAC,1),AC(:,10))); % spend
LOGW_IAC=log(max(ones(SamplesizeAC,1),AC(:,11))); % w_I
LOGW_INXTAC=log(max(ones(SamplesizeAC,1),AC(:,12))); % W_I'
VSAC=AC(:,13); %V_I
LOGD_NCAC=log(max(ones(SamplesizeAC,1),AC(:,32))); % D_INC
LOGW_NCAC=log(max(ones(SamplesizeAC,1),AC(:,33))); % W_INC
LOGTot_NCAC=log(max(ones(SamplesizeAC,1),AC(:,31))); % Tot_INC
RTotD_NCAC=(max(0,LOGD_NCAC).^(-1/2))./LOGTot_NCAC; %% D_INC./Tot_INC
RTotD_NCAC=RTotD_NCAC.*(exp(LOGTot_NCAC)./exp(LOGD_NCAC));
% LOGLOGD_NCAC=log(max(ones(SamplesizeAC,1),LOGD_NCAC));
% LOGLOGTot_NCAC=log(max(ones(SamplesizeAC,1),LOGTot_NCAC));
LOGWNXT_NCAC=log(max(ones(SamplesizeAC,1),AC(:,34))); %W'_INC
Tenure_NCAC=AC(:,35); % Tenure_INC
NCAC=[LOGD_NCAC,LOGW_NCAC,LOGWNXT_NCAC];
LOGD_CAC=log(max(ones(SamplesizeAC,1),AC(:,16))); %Sepnding of Challenger
LOGW_CAC=log(max(ones(SamplesizeAC,1),AC(:,17))); %Warchest of Challenger
LOGW_CNXTAC=log(max(ones(SamplesizeAC,1),AC(:,18))); %Savings of Challenger
YearAC=AC(:,2); % Current year.
PartyAC=3*ones(SamplesizeAC,1)-2*AC(:,3);  %%PartyAC=1  if candidate i is Democrat and PartyAC=-1 if candidate i is Republican.

TenureAC=AC(:,15);

WhiteAC=AC(:,21);
BlackAC=AC(:,22);
OtherAC=AC(:,23);
White_NCAC=AC(:,27);
Black_NCAC=AC(:,28);
Other_NCAC=AC(:,29);
UnemploymentAC=AC(:,24);
Unemployment_NCAC=AC(:,30);

%%B-Spline for q_I(RTotD_NCAC), where RTotD_NCAC=(LOGD_NCAC./LOGTot_NCAC)  etc.%%
%% take knots to be between 0.88 to 1.025 with 8 knots (8 basis functions). (Almost all lie
%% within this range) Let X_Knot be the matrix with 8 columns that contain
%% the value of the B-Spline basis function evaluated at each of the 8
%% basis functions.
% 
% mesh=quantile(LOGLOGD_NCAC,[.125;.25;.375;.5;.625;.75;.875]);
% mesh=[min(LOGLOGD_NCAC);mesh;max(LOGLOGD_NCAC)];
% X_KnotAC1=(LOGLOGD_NCAC<mesh(2,1)).*(1-(LOGLOGD_NCAC-mesh(1,1))/(mesh(2,1)-mesh(1,1)));
% for i=1:6
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
% for i=1:6
%     PLUS=(LOGLOGTot_NCAC>=mesh(i+1,1)).*(LOGLOGTot_NCAC<mesh(i+2,1)).*((LOGLOGTot_NCAC-mesh(i+1,1))/(mesh(i+2,1)-mesh(i+1,1)))...
%         +(LOGLOGTot_NCAC>mesh(i+2,1)).*(LOGLOGTot_NCAC<mesh(i+3,1)).*(1-(LOGLOGTot_NCAC-mesh(i+2,1))/(mesh(i+3,1)-mesh(i+2,1)));
%     X_KnotAC2=[X_KnotAC2,PLUS];
% end
% PLUS=(LOGLOGTot_NCAC>=mesh(8,1)).*(LOGLOGTot_NCAC-mesh(8,1))/(mesh(9,1)-mesh(8,1));
% X_KnotAC2=[X_KnotAC2,PLUS];
mesh=quantile(RTotD_NCAC,[.125;.25;.375;.5;.625;.75;.875]);
mesh=[min(RTotD_NCAC);mesh;max(RTotD_NCAC)];
X_KnotAC1=(RTotD_NCAC<mesh(2,1)).*(1-(RTotD_NCAC-mesh(1,1))/(mesh(2,1)-mesh(1,1)));
for i=0:6
    PLUS=(RTotD_NCAC>=mesh(i+1,1)).*(RTotD_NCAC<mesh(i+2,1)).*((RTotD_NCAC-mesh(i+1,1))/(mesh(i+2,1)-mesh(i+1,1)))...
        +(RTotD_NCAC>=mesh(i+2,1)).*(RTotD_NCAC<mesh(i+3,1)).*(1-(RTotD_NCAC-mesh(i+2,1))/(mesh(i+3,1)-mesh(i+2,1)));
    X_KnotAC1=[X_KnotAC1,PLUS];
end
PLUS=(RTotD_NCAC>=mesh(8,1)).*(RTotD_NCAC-mesh(8,1))/(mesh(9,1)-mesh(8,1));
X_KnotAC1=[X_KnotAC1,PLUS];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%               E_V               %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


E_VContestFUL=find(E_V_july8(:,8)==0);           %% E_VContestFUL is the elements in which contest==0 :277 distinct districts
E_VNContestFUL=find(E_V_july8(:,8)==1);          %% E_VNContestFUL is the elements in which contest==1 :297 distinct districts

LOGTot_NCE_V=log(max(ones(size(E_V_july8,1),1),E_V_july8(:,31)));
LOGD_NCE_V=log(max(ones(size(E_V_july8,1),1),E_V_july8(:,32)));
LOGW_NCE_V=log(max(ones(size(E_V_july8,1),1),E_V_july8(:,33)));
LOGWNXT_NCE_V=log(max(ones(size(E_V_july8,1),1),E_V_july8(:,34)));
Tenure_NCE_V=E_V_july8(:,35); % Tenure_NC


NCE_V=[LOGD_NCE_V,LOGW_NCE_V,LOGWNXT_NCE_V];
NCE_VCT=NCE_V;
NCE_VCT(E_VContestFUL,:)=[];        %NCE_V limited to samples in which there is entry.
NCE_VNCT=NCE_V;
NCE_VNCT(E_VNContestFUL,:)=[];       %NCE_V limited to samples in which there is no entry.
NCE_VTenCT=Tenure_NCE_V;
NCE_VTenCT(E_VContestFUL,:)=[];
NCE_VTenNCT=Tenure_NCE_V;
NCE_VTenNCT(E_VNContestFUL,:)=[];
LOGTot_NCE_VCT=LOGTot_NCE_V;
LOGTot_NCE_VCT(E_VContestFUL,:)=[];
LOGTot_NCE_VNCT=LOGTot_NCE_V;
LOGTot_NCE_VNCT(E_VNContestFUL,:)=[];
Tenure_INCCT=E_V_july8(:,35);
Tenure_INCCT(E_VContestFUL,:)=[];
VSEV=E_V_july8(:,13);
VSEVCT=VSEV;
VSEVCT(E_VContestFUL,:)=[];
VSEVNCT=VSEV;
VSEVNCT(E_VNContestFUL,:)=[];

PartyE_V=3*ones(size(E_V_july8,1),1)-2*E_V_july8(:,3);  %%PartyE_V=1  if candidate i is Democrat and PartyE_V=-1 if candidate i is Republican.
PartyE_VCT=PartyE_V;
PartyE_VCT(E_VContestFUL,:)=[];     %PartyE_VCT is the party indicator when there is entry next period.
PartyE_VNCT=PartyE_V;
PartyE_VNCT(E_VNContestFUL,:)=[];    %PartyE_VCT is the party indicator when there is NO entry next period.

YEARE_V=E_V_july8(:,2);           % YEARE_V is the year of the election.
YEARE_VCT=YEARE_V;
YEARE_VCT(E_VContestFUL,:)=[];
YEARE_VNCT=YEARE_V;
YEARE_VNCT(E_VNContestFUL,:)=[];

XSNCEV_=[E_V_july8(:,30),E_V_july8(:,27)];       % E_V(:,30)=unemployment_NC (%), E_V(:,27)=pctwhite_NC
XSNCEVCT_=XSNCEV_;
XSNCEVCT_(E_VContestFUL,:)=[];
XSNCEVNCT_=XSNCEV_;
XSNCEVNCT_(E_VNContestFUL,:)=[];

XSEV_=[E_V_july8(:,24),E_V_july8(:,21)];    % E_V(:,24)=unemployment (%), E_V(:,21)=pctwhite (%)
XS_EVCT_=XSEV_;
XS_EVCT_(E_VContestFUL,:)=[];
XS_EVNCT_=XSEV_;
XS_EVNCT_(E_VNContestFUL,:)=[];

TenureE_V=E_V_july8(:,15);      %Tenure=# of terms from first period.
TenureE_VCT=TenureE_V;
TenureE_VCT(E_VContestFUL,:)=[];
TenureE_VNCT=TenureE_V;
TenureE_VNCT(E_VNContestFUL,:)=[];

Win_E_V=(E_V_july8(:,13)>=0.5);  % Dummy for whether incumbent won the election in (t).
Win_E_VCT=Win_E_V;
Win_E_VCT(E_VContestFUL,:)=[];
IND5=find(Win_E_VCT==0);  % index for elections with incumbent loss (among those with entry).
                          % Number of distinct districts within contest=0 is 27.
                          % Number of distinct districts within contest=1 is 295.

NCE_VCTwnxt=NCE_VCT;
NCE_VCTwnxt(IND5,:)=[];
YEARE_VCTwnxt=YEARE_VCT;
YEARE_VCTwnxt(IND5,:)=[];
TenureE_VCTwnxt=TenureE_VCT;
TenureE_VCTwnxt(IND5,:)=[];
XS_EVCTwnxt_=XS_EVCT_;
XS_EVCTwnxt_(IND5,:)=[];
PartyE_VCTwnxt=PartyE_VCT;
PartyE_VCTwnxt(IND5,:)=[];
XSNCEVCTwnxt_=XSNCEVCT_;
XSNCEVCTwnxt_(IND5,:)=[];
NCE_VTenCTwnxt=NCE_VTenCT;
NCE_VTenCTwnxt(IND5,:)=[];


LOGTotal_E_V=log(max(ones(size(E_V_july8,1),1),E_V_july8(:,9)));  %log(Total)
LOGTotal_E_VCT=LOGTotal_E_V;
LOGTotal_E_VCT(E_VContestFUL,:)=[];                             %log(Total) limited to contested.
LOGTotal_E_VNCT=LOGTotal_E_V;
LOGTotal_E_VNCT(E_VNContestFUL,:)=[];                            %log(Total) limited to uncontested.

% LOGLOGD_E_V=log(max(ones(size(E_V_july8,1),1),NCE_V(:,1)));
% LOGLOGD_E_VCT=LOGLOGD_E_V;
% LOGLOGD_E_VCT(E_VContestFUL,:)=[]; 
% LOGLOGD_E_VNCT=LOGLOGD_E_V;
% LOGLOGD_E_VNCT(E_VNContestFUL,:)=[];
% LOGLOGD_E_VCTwnxt=LOGLOGD_E_VCT;
% LOGLOGD_E_VCTwnxt(IND5,:)=[];
% LOGLOGTot_E_V=log(max(ones(size(E_V_july8,1),1),LOGTot_NCE_V));
% LOGLOGTot_E_VCT=LOGLOGTot_E_V;
% LOGLOGTot_E_VCT(E_VContestFUL,:)=[];
% LOGLOGTot_E_VNCT=LOGLOGTot_E_V;
% LOGLOGTot_E_VNCT(E_VNContestFUL,:)=[];
% LOGLOGTot_E_VCTwnxt=LOGLOGTot_E_VCT;
% LOGLOGTot_E_VCTwnxt(IND5,:)=[];

RTotDE_V=(max(0,NCE_V(:,1)).^(-1/2))./LOGTot_NCE_V;
RTotDE_V=RTotDE_V.*(exp(LOGTot_NCE_V)./exp(NCE_V(:,1)));
RTotDE_VCT=(max(0,NCE_VCT(:,1)).^(-1/2))./LOGTot_NCE_VCT;
RTotDE_VCT=RTotDE_VCT.*(exp(LOGTot_NCE_VCT)./exp(NCE_VCT(:,1)));
RTotDE_VNCT=(max(0,NCE_VNCT(:,1)).^(-1/2))./LOGTot_NCE_VNCT;
RTotDE_VNCT=RTotDE_VNCT.*(exp(LOGTot_NCE_VNCT)./exp(NCE_VNCT(:,1)));

LOGD_E_V=log(max(ones(size(E_V_july8,1),1),E_V_july8(:,10)));  %log(Spend)
LOGD_E_VCT=LOGD_E_V;
LOGD_E_VCT(E_VContestFUL,:)=[];                             %log(Spend) limited to contested.
LOGD_E_VNCT=LOGD_E_V;
LOGD_E_VNCT(E_VNContestFUL,:)=[];                            %log(Spend) limited to uncontested.



LOGW_E_V=log(max(ones(size(E_V_july8,1),1),E_V_july8(:,11)));  %log(begining cash)
LOGW_E_VCT=LOGW_E_V;
LOGW_E_VCT(E_VContestFUL,:)=[];
LOGW_E_VCTwnxt=LOGW_E_VCT;
LOGW_E_VCTwnxt(IND5,:)=[];
LOGW_E_VNCT=LOGW_E_V;
LOGW_E_VNCT(E_VNContestFUL,:)=[];



LOGW_NXT_E_V=log(max(ones(size(E_V_july8,1),1),E_V_july8(:,12)));  %log(realendcash)
LOGW_NXT_E_VCT=LOGW_NXT_E_V;
LOGW_NXT_E_VCT(E_VContestFUL,:)=[];
LOGW_NXT_E_VCTwnxt=LOGW_NXT_E_VCT;
LOGW_NXT_E_VCTwnxt(IND5,:)=[];
LOGW_NXT_E_VNCT=LOGW_NXT_E_V;
LOGW_NXT_E_VNCT(E_VNContestFUL,:)=[];
IND6CT=find(LOGW_NXT_E_VCT==0);                      % IND6CT is the index where realendcash_i ==0
IND6NCT=find(LOGW_NXT_E_VNCT==0);                    % IND6NCT is the index where realendcash_i ==0

LOGW_NXT_E_VC=log(max(ones(size(E_V_july8,1),1),E_V_july8(:,18)));  %log(realendcash_c)
LOGW_NXT_E_VC(E_VContestFUL,:)=[];
LOGD_E_VC=log(max(ones(size(E_V_july8,1),1),E_V_july8(:,16)));  %log(realdisburse_c)
LOGD_E_VC(E_VContestFUL,:)=[];
LOGTot_E_VC=log(max(ones(size(E_V_july8,1),1),E_V_july8(:,4)));  % log(realtotal_c)
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
mesh=quantile(RTotDE_V,[.125;.25;.375;.5;.625;.75;.875]);
mesh=[min(RTotDE_V);mesh;max(RTotDE_V)];
X_KnotEV1=(RTotDE_V<mesh(2,1)).*(1-(RTotDE_V-mesh(1,1))/(mesh(2,1)-mesh(1,1)));
for i=0:6
    PLUS=(RTotDE_V>=mesh(i+1,1)).*(RTotDE_V<mesh(i+2,1)).*((RTotDE_V-mesh(i+1,1))/(mesh(i+2,1)-mesh(i+1,1)))...
        +(RTotDE_V>=mesh(i+2,1)).*(RTotDE_V<mesh(i+3,1)).*(1-(RTotDE_V-mesh(i+2,1))/(mesh(i+3,1)-mesh(i+2,1)));
    X_KnotEV1=[X_KnotEV1,PLUS];
end
PLUS=(RTotDE_V>=mesh(8,1)).*(RTotDE_V-mesh(8,1))/(mesh(9,1)-mesh(8,1));
X_KnotEV1=[X_KnotEV1,PLUS];

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

% theta0=[thetaS;thetaQ;theta;B_I;B_C;Q_C1;Q_C2;B_T;E_VCTa;E_VCTt;gammaCT;...
%     thetadNC;thetatNC;thetanxtNC;thetawin;cost1;ben1;sig;cost2;ben2];
% theta0(1:23,1)=initialvalue(1:23,1);
% theta0(24:42,1)=initialvalue(24:42,1);  %ACDa
% theta0(43:61,1)=initialvalue(43:61,1);  %ACDt
% theta0(62:80,1)=initialvalue(62:80,1);  %gammaACD
% theta0(81:86,1)=initialvalue(81:86,1); %ACdqea
% theta0(87:92,1)=initialvalue(87:92,1); %ACdqet
% theta0(93:98,1)=initialvalue(93:98,1); % gammadqe
% theta0(99:113,1)=initialvalue(99:113,1); % E_VCTa
% theta0(114:128,1)=initialvalue(114:128,1); %E_VCTt
% theta0(129:143,1)=initialvalue(129:143,1); %gammaCT
% theta0(144:158,1)=initialvalue(144:158,1); % E_VNCTa
% theta0(159:173,1)=initialvalue(159:173,1); % E_VNCTt
% theta0(174:188,1)=initialvalue(174:188,1); % gammaNCT
% theta0(189:194,1)=initialvalue(189:194,1);  % ACwqea
% theta0(195:200,1)=initialvalue(195:200,1);  % ACwqet
% theta0(201:206,1)=initialvalue(201:206,1);  % gammawqe
% theta0(207:213,1)=initialvalue(207:213,1);  %ACvqea
% theta0(214:220,1)=initialvalue(214:220,1);   %ACvqet
% theta0(221:227,1)=initialvalue(221:227,1);  %gammavqe
% theta0(228:234,1)=initialvalue(228:234,1);  % thetawin
% theta0(235,1)=initialvalue(235,1);  % beta_1
% theta0(236,1)=initialvalue(236,1);  % beta_2
% theta0(237,1)=initialvalue(237,1);  % sig


thetain1step=Fststage_kettei;

thetaS=thetain1step(1:2,1);
thetaS(1,1)=1;
thetaS2=thetain1step(22:23,1);

B_T=thetain1step(11,1);
theta=thetain1step(12:16,1);
theta2=thetain1step(17:21,1);
menc=32;
B_I=thetain1step(menc+1,1);
B_C=thetain1step(menc+2,1);
Q_C1=thetain1step(menc+3,1);
Q_C2=thetain1step(menc+4,1);
Q_C3=thetain1step(menc+5,1);
E_VCTa(:,1)=thetain1step(menc+6:menc+20,1);
E_VCTa(:,2)=thetain1step(menc+21:menc+35,1);
E_VCTa(:,3)=thetain1step(menc+36:menc+50,1);
E_VCTt(:,1)=thetain1step(menc+51:menc+65,1);
E_VCTt(:,2)=thetain1step(menc+66:menc+80,1);
E_VCTt(:,3)=thetain1step(menc+81:menc+95,1);
gammaCT(:,1)=thetain1step(menc+96:menc+110,1);
gammaCT(:,2)=thetain1step(menc+111:menc+125,1);
gammaCT(:,3)=thetain1step(menc+126:menc+140,1);
mend=172;
gammaNCT(:,1)=thetain1step(mend+1:mend+11,1);
tNCT(:,1)=thetain1step(mend+12:mend+22,1);
wNCT(:,1)=thetain1step(mend+23:mend+33,1);
thetawin=thetain1step(mend+34:mend+40,1);



iterate=1;

thetaQ(1,1)=1;
%     for ii=2:size(X_Knot1,2)
%         thetaQ(ii,1)=1-sum(abs(Est1(3:(ii+1),1)));
%     end
thetaQ(2:9,1)=thetain1step(3:10,1);
thetaQ2(1,1)=thetain1step(24,1);
%     for ii=2:size(X_Knot1,2)
%         thetaQ2(ii,1)=thetaQ2(1,1)-sum(abs(thetain(4:(ii+2),1)));
%     end
thetaQ2(2:9,1)=thetain1step(25:32,1);


XS_EVCT2=XS_EVCT_*thetaS2;
XS_EVCT=XS_EVCT_*thetaS;
% XSEV=XSEV_*thetaS;
XQ2=X_Knot1*thetaQ2;
XQEV2=X_KnotEV1*thetaQ2;  % q_I
XQEVCT2=XQEV2;
XQEVCT2(E_VContestFUL,:)=[];
XQEV=X_KnotEV1*thetaQ;  % q_I2
XQEVCT=XQEV;
XQEVCT(E_VContestFUL,:)=[];



q_c3step=q_C_E_VCT;
q_c3step(IND6CT,:)=[];    %% Challenger quality.
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
XS3step_(IND6CT,:)=[];
XS3step_(DEL_Pen,:)=[];
Party3step=(-1)*PartyE_VCT;  %% Challenger party=(-1)*Incumbent party
Party3step(IND6CT,:)=[];
Party3step(DEL_Pen,:)=[];
Tenure3step=zeros(size(Party3step,1),1);               %% Challenger starts with Tenure==0;
LOGW_NXT3step=LOGW_NXT_E_VC;
LOGW_NXT3step(IND6CT,:)=[];                                    %% challenger warchest.
LOGW_NXT3step(DEL_Pen,:)=[];
LOGW_E_VCT3step=LOGW_E_VCT;
LOGW_E_VCT3step(IND6CT,:)=[];
LOGW_E_VCT3step(DEL_Pen,:)=[];
TenureE_VCT3step=TenureE_VCT;
TenureE_VCT3step(IND6CT,:)=[];
TenureE_VCT3step(DEL_Pen,:)=[];

LOGD_E_VC3step=LOGD_E_VC;
LOGD_E_VC3step(IND6CT,:)=[];            %challenger spending
LOGD_E_VC3step(DEL_Pen,:)=[];
LOGTot_E_VC3step=LOGTot_E_VC;
LOGTot_E_VC3step(IND6CT,:)=[];              %challenger tot
LOGTot_E_VC3step(DEL_Pen,:)=[];
VSEVCT3step=VSEVCT;
VSEVCT3step(IND6CT,:)=[];
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
%% use LOGW_EVCT, etc, variables for incumbents %%
PrE=normcdf([ones(size(LOGW_E_VCT,1),1),LOGW_E_VCT,XQEVCT,XS_EVCT.*PartyE_VCT,log(TenureE_VCT+1)]*theta,0,1);
PrE3step=PrE;
PrE3step(IND6CT,:)=[];
PrE3step(DEL_Pen,:)=[];
E_PrimaryN=[ones(size(LOGW_E_VCT,1),1),LOGW_E_VCT,XQEVCT,XS_EVCT.*PartyE_VCT,log(TenureE_VCT+1)]*theta2;
E_PrimaryN3step=E_PrimaryN;
E_PrimaryN3step(IND6CT,:)=[];
E_PrimaryN3step(DEL_Pen,:)=[];


%BX1=abs(1/Sndstage_kettei_long(5,1))*(B_I*LOGD_E_VCT+B_C*LOGD_E_VC+XS_EVCT2.*PartyE_VCT+B_T*log(TenureE_VCT+1)+XQEVCT2+q_C_E_VCT);
BX13step=BX1;
BX13step(IND6CT,:)=[];
BX13step(DEL_Pen,:)=[];
thetain2=Sndstage_kettei_long;
% Estim3=[   -1.9713322e+01
%   -1.1704758e-02
%    1.6123112e+03
%   -7.4047131e+033]

for i=1:5
thetain3=Trdstage_initial(1:4,i);
[mintheta,SRR]=fminsearch(@(theta3) Minimize3PS_long(thetain1step,thetain2,theta3),thetain3)
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
