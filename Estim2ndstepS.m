clear
% **
% make sure you get the sign of q_e right. Also note that f_d_qe, f_w_qe have the PARTY of Challenger.
% **

rng(958)

%Globals

%Calculating C matrix
%global N
%global interest



E_V_july8=csvread('E_Vjuly8partisan.csv',1,0);

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
load('./residvar.mat');
load ('./retire.txt');

load ('./stateevol.mat')

%Construct dataset
dele3=find(sum(isnan(E_V_july8),2)>0); 
E_V_july8(dele3,:)=[]; %Drop NaN
E_V_july8(E_V_july8(:,2)>2002,:)=[]; %Drop year 2004 and on
deleC=find((E_V_july8(:,16)==0).*E_V_july8(:,8)==1);
E_V_july8(deleC,:)=[]; %Drop if oppornent disburse=0 and contested
E_V_july8(174:177,:)=[]; %Drop an Rtotd outlier
E_V_july8(11:13,:)=[]; %Drop an Rtotd outlier

deleC2=find(E_V_july8(:,12)<5000|E_V_july8(:,12)>1200000);
%deleC2=find(E_V_july8(:,12)<5000);
E_V_july8(deleC2,:)=[];%Drop if nonpositive saving
deleC3=find(E_V_july8(:,11)<5000);
E_V_july8(deleC3,:)=[];%Drop if nonpositive begcash
%deleC4=find((E_V_july8(:,16)<5000)&E_V_july8(:,8)==1);
deleC4=find((E_V_july8(:,16)<5000|E_V_july8(:,16)>1200000)&E_V_july8(:,8)==1);
E_V_july8(deleC4,:)=[];%Drop if nonpositive challenger spending

%Drop RTotD outlier
LOGTot_NCE_V=log(max(ones(length(E_V_july8),1),E_V_july8(:,31)));
LOGD_NCE_V=log(max(ones(length(E_V_july8),1),E_V_july8(:,32)));
RTotDE_V=((max(0,LOGD_NCE_V).^(-1/2))./LOGTot_NCE_V).*(exp(LOGTot_NCE_V)./exp(LOGD_NCE_V));

dd=find(RTotDE_V>0.07);
E_V_july8(dd,:)=[]; %Drop an Rtotd outlier
%E_V_july8([25,214,662,678,679,716,717],:)=[]; %Drop an Rtotd outlier


E_VContestFUL=find(E_V_july8(:,8)==0);           %% E_VContestFUL is the elements in which contest==0 :277 distinct districts
E_VNContestFUL=find(E_V_july8(:,8)==1);          %% E_VNContestFUL is the elements in which contest==1 :297 distinct districts

Win_E_V=(E_V_july8(:,13)>=0.5);  % Dummy for whether incumbent won the election in (t).
Win_E_VCT=Win_E_V;
Win_E_VCT(E_VContestFUL,:)=[];
IND5=find(Win_E_VCT==0);  % index for elections with incumbent loss (among those with entry).
                          % Number of distinct districts within contest=0 is 27.
                          % Number of distinct districts within contest=1 is 295.

contest=(E_V_july8(:,8)==1);

PartyE_V=3*ones(length(E_V_july8),1)-2*E_V_july8(:,3);  %%PartyE_V=1  if candidate i is Democrat and PartyE_V=-1 if candidate i is Republican.
PartyE_VCT=PartyE_V;
PartyE_VCT(E_VContestFUL,:)=[];     %PartyE_VCT is the party indicator when there is entry next period.
% PartyE_VCTwnxt=PartyE_VCT;
% PartyE_VCTwnxt(IND5,:)=[];
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

YEARE_V=E_V_july8(:,2);           % YEARE_V is the year of the election.

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


LOGTot_NCE_V=log(max(ones(length(E_V_july8),1),E_V_july8(:,31)));
LOGD_NCE_V=log(max(ones(length(E_V_july8),1),E_V_july8(:,32)));
LOGW_NCE_V=log(max(ones(length(E_V_july8),1),E_V_july8(:,33)));
LOGWNXT_NCE_V=log(max(ones(length(E_V_july8),1),E_V_july8(:,34)));
% Tenure_NCE_V=E_V_july8(:,35); % Tenure_NC


NCE_V=[LOGD_NCE_V,LOGW_NCE_V,LOGWNXT_NCE_V];
NCE_VCT=NCE_V;
NCE_VCT(E_VContestFUL,:)=[];        %NCE_V limited to samples in which there is entry.
NCE_VNCT=NCE_V;
NCE_VNCT(E_VNContestFUL,:)=[];       %NCE_V limited to samples in which there is no entry.

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
% within this range) Let X_Knot be the matrix with 8 columns that contain
% the value of the B-Spline basis function evaluated at each of the 8
% basis functions.
% 
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



idi=E_V_july8(:,66);
idc=E_V_july8(:,67);

coefentry=Est1(1:25);
coefprimaryN=Est1(26:50);
Eentry=[ones(length(LOGW_E_V),1),LOGW_E_V,LOGW_E_V.^2,XSEV_(:,1).*SameE_V,XSEV_(:,2).*PartyE_V,log(TenureE_V+1),MidtermE_V,PresdumE_V,MidtermE_V.*PresdumE_V,X_KnotEV1,X_KnotEV1.*(repmat(LOGW_E_V,1,size(X_KnotEV1,2)))]*coefentry;
EprimaryN=[ones(length(LOGW_E_V),1),LOGW_E_V,LOGW_E_V.^2,XSEV_(:,1).*SameE_V,XSEV_(:,2).*PartyE_V,log(TenureE_V+1),MidtermE_V,PresdumE_V,MidtermE_V.*PresdumE_V,X_KnotEV1,X_KnotEV1.*(repmat(LOGW_E_V,1,size(X_KnotEV1,2)))]*coefprimaryN;
coefvoteshare=Est2;


%%%%%%%%%%%%%%%%%%%%%%%% 
% Defining variables: done
%%%%%%%%%%%%%%%%%%%%%%%%




%Estimated parameters in the first stage
coefentry=Est1(1:(length(Est1)/2));
%Saving and fundraising parameters when contested
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

thetaQ2=Est2(numel(Est2)-4:numel(Est2));
%    if numel(thetaQ2)==3
aux=[RTotDE_V,RTotDE_V.^2*100,RTotDE_V.^3*1000,RTotDE_V.^4*10000,RTotDE_V.^5*100000]*thetaQ2;
%Normalize estimated XQEV: Mean zero
XQEV2=aux-mean(aux); % q_I

XQEVCT2=XQEV2;
XQEVCT2(E_VContestFUL,:)=[];
XQEVNCT2=XQEV2;
XQEVNCT2(E_VNContestFUL,:)=[];
XQEVCTwnxt2=XQEVCT2;
XQEVCTwnxt2(IND5,:)=[];


   matchXQEV2=zeros(length(idi),1);
for i=1:(length(idi)-1)
    matchXQEV2(i)=(idc(i)==idi(i+1))*XQEV2(i+1);
end
matchXQEVCT2=matchXQEV2;
matchXQEVCT2(E_VContestFUL,:)=[];


%assuming constant term of vote share eq=0, constant term of first stage +
%normalization of XQEV = - constant term of E(qc|X).
const3=Est2(1)+mean(aux);
expectedqc=-[EprimaryN,Eentry,EprimaryN.*Eentry,EprimaryN.^2,Eentry.^2]*coefvoteshare(10:14)-const3;
expectedqcCT=expectedqc;
expectedqcCT(E_VContestFUL,:)=[];

%Result from the first stage till here.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Simulation draws preparation
%Sofarbest=10^8;
bestiter=0;
iterate=1;
NN=100;                      % NN is the number of simulations.
NumSim=50;              %NumSim is the number of simulations in the 1st step to obtain f(.,q_e)

T=10;                      % T is the number of periods that we move the simlation forward.
interest=0.1;               %interest on money saved.
Entry=rand(T,NN,length(E_V_july8));             %Simulation draw for computing the continuation value E_V
Winrnd=rand(T,NN,length(E_V_july8));            %Simulation draw for computing the continuation value E_V (for challenger FOC)

Retirernd=(rand(T-1,NN,length(E_V_july8))<repmat(retire',[1 NN, length(E_V_july8)]));
Retirernd(T,:,:)=1;
% RetirerndC=(rand(T-1,N,length(E_V_july8))<repmat(retire',[1 N, length(E_V_july8)]));
% RetirerndC(T,:,:)=1;
% RetirerndR=(rand(T-1,N,length(E_V_july8))<repmat(retire',[1 N, length(E_V_july8)]));
% RetirerndR(T,:,:)=1;


for i=1:NN
    for j=1:length(E_V_july8)
        Ret(i,j)=find(Retirernd(:,i,j)==1,1,'first');
    end
end
% for i=1:N
%     for j=1:length(E_V_july8)
%         RetC(i,j)=find(RetirerndC(:,i,j)==1,1,'first');
%     end
% end
% for i=1:N
%     for j=1:length(E_V_july8)
%         RetR(i,j)=find(RetirerndR(:,i,j)==1,1,'first');
%     end
% end



dF_gamma_ct=rand(T,NN,length(E_V_july8));
dF_total_ct=rand(T,NN,length(E_V_july8));
dF_gammasv_ct=rand(T,NN,length(E_V_july8));
dF_totalsv_ct=rand(T,NN,length(E_V_july8));
dF_nxt_nxt_ct=rand(T,NN,length(E_V_july8));



%Draw evolution of X
%Parameters governing evolution of partisanship already adjusted to
%two-year window.
Shockump=sqrt(epsump)*randn(T,NN,length(E_V_july8));
Shockpartisan=sqrt(epspartisan)*randn(T,NN,length(E_V_july8));

Shockump2=sqrt(epsump)*randn(T,NN,length(E_V_july8));
Shockpartisan2=sqrt(epspartisan)*randn(T,NN,length(E_V_july8));

Shockump3=sqrt(epsump)*randn(T,NN,length(E_V_july8));
Shockpartisan3=sqrt(epspartisan)*randn(T,NN,length(E_V_july8));

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
presseqmidone=rand(T,NN,length(midone));

%I need to keep track of whether the current president can run for
%reelection or not.
presyearmidone=ones(T,NN,length(midone));

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
    for n=1:NN
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
presseqmidtwo=rand(T,NN,length(midtwo));

%I need to keep track of whether the current president can run for
%reelection or not.
presyearmidtwo=ones(T,NN,length(midtwo));

%Midterm election: No change in incumbency
presseqmidtwo([2,4,6,8,10],:,:)=0;
%Today:incumbent cannot run for reelection. 1 if party changes.
presseqmidtwo(1,:,:)=(presseqmidtwo(1,:,:)>newele);
presyearmidtwo(1,:,:)=1;
%At election i, if incumbent year=2,  new election. If year=1,
%election with incumbency.
for k=1:length(midtwo)
    for n=1:NN
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
presseqbegone=rand(T,NN,length(begone));

%I need to keep track of whether the current president can run for
%reelection or not.
presyearbegone=ones(T,NN,length(begone));

%If previous election is midterm election: No change in incumbency
presseqbegone([1,3,5,7,9],:,:)=0;
%2nd election=incumbent runs for reelection at 1st election. 1 if incumbent loses, 0 otherwise.
presseqbegone(2,:,:)=(presseqbegone(2,:,:)>incwin);
presyearbegone(2,:,:)=1+(1-presseqbegone(2,:,:));
%At election i, if incumbent year=2,  new election. If year=1,
%election with incumbency.
for k=1:length(begone)
    for n=1:NN
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
presseqbegtwo=rand(T,NN,length(begtwo));

%I need to keep track of whether the current president can run for
%reelection or not.
presyearbegtwo=ones(T,NN,length(begtwo));

%begterm election: No change in incumbency
presseqbegtwo([1,3,5,7,9],:,:)=0;
%2nd election=incumbent cannot run for reelection. 1 if party changes.
presseqbegtwo(2,:,:)=(presseqbegtwo(2,:,:)>newele);
presyearbegtwo(2,:,:)=1;
%At election i, if incumbent year=2,  new election. If year=1,
%election with incumbency.
for k=1:length(begtwo)
    for n=1:NN
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
presseq=zeros(T,NN,length(E_V_july8));
presseq(:,:,midone)=presseqmidone;
presseq(:,:,midtwo)=presseqmidtwo;
presseq(:,:,begone)=presseqbegone;
presseq(:,:,begtwo)=presseqbegtwo;

%Save it
save presseq presseq






%%


%Create C matrix

% C=zeros(5,T,NN,length(NCE_V));
% for i=1:length(NCE_V)
%          C(:,:,:,i)=Actions(i,XSEV_(i,:)',SameE_V(i,1), XQEV2(i,1),X_KnotEV1(i,:),RTotDE_V(i,:), TenureE_V(i,1),LOGW_NXT_E_V(i,1),PartyE_V(i,1),PresdumE_V(i,1),MidtermE_V(i,1),...
%          Shockpartisan(:,:,i), Shockump(:,:,i), presseq(:,:,i),Entry(:,:,i), coefentry, E_VCTa, E_VCTt, gammaCT, coefspend, coeffund,coefsave, dF_gamma_ct(:,:,i), dF_total_ct(:,:,i),...
%        dF_gammasv_ct(:,:,i),dF_totalsv_ct(:,:,i),Winrnd(:,:,i),coefprobwin,Ret(:,i)',Betapartisan,Betaump,NN,interest);
%     squeeze(C(:,:,:,i));
%     i
% end
% 
%   save Cfull.mat C

%%

%Create data set to estimate second stage with.
 load('./C.mat');
% C(:,:,:,dd)=[];
%  C(:,:,:,[25,214,662,678,679,716,717])=[];
 C=C(:,:,1:100,:);

%Squeeze dimensions so that we can calculate cont payoff using matrix form
%later.
%Make sure C(2,:,:,:) is nonnegative
C(2,:,:,:)=max(0,C(2,:,:,:));
C1=squeeze(C(1,:,:,:));
C2=squeeze(C(2,:,:,:));
C3=squeeze(C(3,:,:,:));
C5=squeeze(C(5,:,:,:));

% load('./est2ndstage.mat');
% inittheta=mintheta2ndstage;

%Compute probability of winning from step 1-3
probwin=[ones(length(NCE_VCT),1),LOGW_E_VCT,XS_EVCT_(:,1).*SameE_VCT,XS_EVCT_(:,1).^2.*SameE_VCT,XS_EVCT_(:,2).*PartyE_VCT,LOGW_E_VCT.^2,X_KnotE_VCT,log(TenureE_VCT+1),log(TenureE_VCT+1).^2,...
    X_KnotE_VCT.*(XS_EVCT_(:,1).*SameE_VCT*ones(1,8)),X_KnotE_VCT.*(XS_EVCT_(:,2).*PartyE_VCT*ones(1,8)),X_KnotE_VCT.*(log(TenureE_VCT+1)*ones(1,8)),...;
  LOGW_E_VCT.^3]*Est3;


%Input for the function evaluated
datasetV=[NCE_V, LOGTot_NCE_V,XQEV2,LOGW_NXT_E_V, SameE_V,PartyE_V,TenureE_V, XSEV_,PresdumE_V,MidtermE_V];
datasetVCT=[E_VNContestFUL,XQEVCT2,LOGTot_NCE_VCT, NCE_VCT,LOGW_NXT_E_VCT,LOGTotal_E_VCT,VSEVCT,XS_EVCT_,TenureE_VCT, PartyE_VCT,SameE_VCT, LOGD_E_VCT,LOGD_E_VC,PresdumE_VCT,MidtermE_VCT,matchXQEVCT2];
datasetVNCT=[E_VContestFUL,LOGW_NXT_E_VNCT,LOGTotal_E_VNCT,NCE_VNCT,LOGTot_NCE_VNCT,XQEVNCT2,LOGD_E_VNCT,TenureE_VNCT,XS_EVNCT_, PartyE_VNCT,SameE_VNCT,PresdumE_VNCT,MidtermE_VNCT];


%Create delta matrix
vdelta=0.90;
deltamat=zeros(T,NN,length(NCE_V));
for j=1:T
    deltamat(j,:,:)=((vdelta)^(j-1));
end



datasetVC=[LOGTot_E_VC,LOGW_NXT_E_VC,LOGD_E_VC];
%%
%Start from random points and see the likelihood
rngset=125;
rng(rngset)
nsample=1000000;
%  load output32
%  paratest=0.8*repmat(mintheta2ndstage.',nsample,1)+0.4*repmat(mintheta2ndstage.',nsample,1).*rand(nsample,7);
% paratest(:,2)=paratest(:,2)+10*rand(nsample,1);
% paratest(:,4)=3+rand(nsample,1);
 paratest=zeros(nsample,7);
paratest(:,1)=0.3*rand(nsample,1);
paratest(:,2)=1*rand(nsample,1);
paratest(:,3)=0.01+0.06*rand(nsample,1);
paratest(:,4)=1/10*paratest(:,1).*rand(nsample,1);
paratest(:,5)=0.3*rand(nsample,1);
paratest(:,6)=1+2*rand(nsample,1);
paratest(:,7)=paratest(:,6).*rand(nsample,1);
paratest(:,8)=(0.5+rand(nsample,1)).*paratest(:,1);
paratest(:,9)=(0.5+rand(nsample,1)).*paratest(:,6);
%paratest(:,10)=(0.5+rand(nsample,1)).*paratest(:,4);
%paratest(:,9)=1+3*rand(nsample,1);

const0=0;

funvalue=zeros(nsample,1);
rsq=zeros(nsample,1);
for i=1:nsample
    i
[funvalue(i),rsq(i)]=Minimize2S_new(Est2,paratest(i,:).',probwin,datasetV,datasetVCT,datasetVNCT,datasetVC,deltamat,const0,NN,C1,C2,C3,C5,expectedqcCT,residvar,T);
end
funvalue=funvalue(funvalue>0);
rsq=rsq(rsq>0);
paratest=paratest(funvalue>0,:);

tt=[funvalue,rsq,paratest];
tt2=sortrows(tt);
funvalue=tt2(1:100,1);
rsq=tt2(1:100,2);
paratest=tt2(1:100,3:end);

save testfunvalue7 funvalue rsq paratest rngset
%%
%Combine multiple testfunvalue
likelimat0=zeros(1,11);
for i=7
     filename=sprintf('testfunvalue%d.mat', i);
     if exist(filename, 'file') == 2
     load(filename)
     aux=[funvalue,paratest,rsq];
     aux2=aux(funvalue>0,:);
     likelimat0=[likelimat0;aux2];
     end
end
aux=sortrows(likelimat0);
likelimatfull=aux(2:end,:);
likelimat=likelimatfull(:,1);
parasetmat=likelimatfull(:,2:10);
rsqmat=likelimatfull(:,11);

save testfunvalue likelimat parasetmat rsqmat


%%
%Minimizer

%         load testfunvalue
%       inittheta=parasetmat(13,1:9).';
% % 
          load output2_8para
          inittheta=mintheta2ndstage;
%      inittheta(4)=0.0001; 
%     inittheta(8)=inittheta(7); 
% inittheta(7)=0.5; 
% mintheta=[];
% SRR=[];
 options=optimset('MaxIter',600,'MaxFunEvals',1000000,'Display','iter');
  %options=optimoptions('patternsearch','MaxIter',100,'MaxFunEvals',1000000,'Display','iter');
% iter=1;  
%  Sofarbest=10^8;
%    minthetadist=zeros(8,iter);
%    initthetadist=zeros(8,iter);
%    srrdist=zeros(1,iter);
%    tic
 %for i=1:iter

 theta0=inittheta;
 const0=0;
%theta0(6)=theta0(1);
 funvalue=Minimize2S_new(Est2,theta0,probwin,datasetV,datasetVCT,datasetVNCT,datasetVC,deltamat,const0,NN,C1,C2,C3,C5,expectedqcCT,residvar,T)
 i=11;

funvalueend=0;
thetaup=theta0;
sss=1;
while abs(funvalue-funvalueend)>10^(-5)
   sss=sss+1;
   funvalue=Minimize2S_new(Est2,thetaup,probwin,datasetV,datasetVCT,datasetVNCT,datasetVC,deltamat,const0,NN,C1,C2,C3,C5,expectedqcCT,residvar,T);
 [mintheta,SRR]=fminsearch(@(thetaup) Minimize2S_new(Est2,thetaup,probwin,datasetV,datasetVCT,datasetVNCT,datasetVC,deltamat,const0,NN,C1,C2,C3,C5,expectedqcCT,residvar,T),thetaup,options);
 % [mintheta,SRR]=patternsearch(@(thetaup) Minimize2S_new(Est2,thetaup,probwin,datasetV,datasetVCT,datasetVNCT,datasetVC,deltamat,const0,NN,C1,C2,C3,C5,expectedqcCT,residvar,T),thetaup,[],[],[],[],[],[],[],options);
    funvalueend=Minimize2S_new(Est2,mintheta,probwin,datasetV,datasetVCT,datasetVNCT,datasetVC,deltamat,const0,NN,C1,C2,C3,C5,expectedqcCT,residvar,T);
%       if mod(sss,30)==1
%       thetaup=mintheta+0.5*(rand(size(mintheta,1),1)-0.5).*mintheta;
%       else
        thetaup=mintheta;
%       end
     
        m=matfile(sprintf('output%d.mat', i),'writable',true);
        m.mintheta2ndstage=thetaup;
        m.objsize=funvalueend;
        m.sss=sss;
end
%  minthetadist(:,i)=mintheta;
%  srrdist(i)=SRR;
%  initthetadist(:,i)=theta0;
% %theta2ndstep=aux;
% %save minimizedtheta.txt mintheta SRR -ASCII
% % 
% % end
% %  minthetadist=minthetadist(:,srrdist<mean(srrdist));
% %  initthetadist=initthetadist(:,srrdist<mean(srrdist));
% %   srrdist=srrdist(srrdist<mean(srrdist));
% mintheta2ndstage=mintheta;
%  toc
  
  %Target=1.52695
  

%%
  
%%%%%%  
% Specification check and saving the result
%%%%%%
%load est2ndstage
 %theta2=[mintheta2ndstage;0.5;2];
load output2_8para
theta2=mintheta2ndstage;
const0=0;







vdelta=0.90;
% thetaS=Est2(1:2,1);
cc=const0;
B_I=Est2(2,1);
B_C=Est2(3,1);
B_T=Est2(9,1);
thetaS2=Est2(4:8,1);


if length(theta2)==8
cost1=abs(theta2(1,1));  %% coefficient on the cost function of incumbent, contested
ben1=abs(theta2(2,1));   %% coefficient on the benefit function of incumbent, contested
ben2=ben1;
benc=ben1;


alpha=abs(theta2(7,1));
alphac=abs(theta2(8,1));
beta=abs(theta2(6,1));

cost2=abs(theta2(5,1));
sig=abs(theta2(3,1));
a=abs(theta2(4,1));

costc=cost1;
ac=a;
betac=beta;
end
if length(theta2)==9
cost1=abs(theta2(1,1));  %% coefficient on the cost function of incumbent, contested
ben1=abs(theta2(2,1));   %% coefficient on the benefit function of incumbent, contested
ben2=ben1;
benc=ben1;


alpha=abs(theta2(7,1));
alphac=alpha;
beta=abs(theta2(6,1));

cost2=abs(theta2(5,1));
sig=abs(theta2(3,1));
a=abs(theta2(4,1));


costc=abs(theta2(8,1));
ac=a;
betac=abs(theta2(9,1));
end


BB=1;
%Continuation(ST,q_I,Ten,w_I,epswh, epsump, E_VCTa, E_VCTt, gamma,dF_gamma_ct, dF_total_ct, dF_nxt_nxt,Winrnd,thetawin,Ret,Betawh,Betaump,cost,ben,thetaS,Party)
%%%%%%%%          E_V         %%%%%%%%%
% Given State, Tenure, warchest, compute incumbent continuation value %%
%Contcontest=payoff if contested

%Costpara with Rtotd
%costpara=permute(repmat((ones(N,1)*((exp(LOGTot_NCE_V(:,1))./exp(NCE_V(:,1))).*...
%                    ((max(0,NCE_V(:,1)).^(alpha-1))./(max(0,LOGTot_NCE_V(:,1)).^(beta-1)))).'),[1 1 10]),[3 1 2]);

%Costpara with flexible XQEV2 with constant
costpara=permute(repmat((ones(NN,numel(XQEV2))+a*ones(NN,1)*(1./exp(XQEV2).')),[1 1 10]),[3 1 2]);                

%Without constant
%costpara=permute(repmat((ones(N,1)*XQEV2.'+a*ones(N,1)*XQEV2.^2*10.'+b*ones(N,1)*XQEV2.^3*100.'),[1 1 10]),[3 1 2]);   

aux1=C2.^alpha;
aux2=(1/beta)*costpara.*C3.^beta;
Contcontest=(ben1*aux1-cost1*aux2).*C1+BB*C5; %Benefit
Contuncontest=(ben2*aux1-cost2*aux2).*(1-C1); %Benefit


%Discount

% for j=1:T
%     Contcontest(:,j,:,:)=((vdelta)^(j-1))*Contcontest(:,j,:,:);
%     Contuncontest(:,j,:,:)=((vdelta)^(j-1))*Contuncontest(:,j,:,:);
% end

%Sum over benefit+cost+winning and over t=1~T
Contcontestpath=squeeze(sum(deltamat.*Contcontest,1));
Contuncontestpath=squeeze(sum(deltamat.*Contuncontest,1));

%Sum over contest and uncontest
Continue1=(Contcontestpath+Contuncontestpath).';

%This gives us sample size* number of simulation matrix of continuation
%payoff



Continue=mean(Continue1,2);
Continuation1=Continue;
Continuation1(E_VContestFUL,:)=[];           %% Continuation1 is the Continuation value for periods in which incumebent is contested.
Continuation2=Continue;
Continuation2(E_VNContestFUL,:)=[];          %%  Continuation1 is the Continuation value for periods in which incumebent is uncontested.


regXS=[XSEV_(:,1),XSEV_(:,2),XSEV_(:,1).^2,XSEV_(:,2).^2,log(XSEV_(:,1)),XSEV_(:,1).*XSEV_(:,2),...
    XSEV_(:,1).^2.*SameE_V,XSEV_(:,2).^2.*PartyE_V,SameE_V,PartyE_V];
nregXS=size(regXS,2);
%Regress Continue on State variables to find the derivative. %
Regressand0=[ones(length(NCE_V),1),LOGW_NXT_E_V,LOGW_NXT_E_V.^2/10,LOGW_NXT_E_V.^3/100,...
    XQEV2,LOGW_NXT_E_V.*XQEV2,...    
    exp(XQEV2),...
    LOGW_NXT_E_V.*TenureE_V,XQEV2.*TenureE_V,...
    LOGW_NXT_E_V.*PresdumE_V,XQEV2.*PresdumE_V,...
    LOGW_NXT_E_V.*XSEV_(:,1).*SameE_V,LOGW_NXT_E_V.*XSEV_(:,2).*PartyE_V,...
    XQEV2.*XSEV_(:,1).*SameE_V,...
    TenureE_V,TenureE_V.^2,XSEV_(:,1).*SameE_V,XSEV_(:,1).^2.*SameE_V,XSEV_(:,2).*PartyE_V,XSEV_(:,2).^2.*PartyE_V,...
    PresdumE_V,MidtermE_V,PresdumE_V.*MidtermE_V,...
    SameE_V.*XSEV_(:,2),SameE_V.*PartyE_V,SameE_V.*TenureE_V,SameE_V.*PresdumE_V,SameE_V.*MidtermE_V,...
    PartyE_V.*XSEV_(:,1),PartyE_V.*TenureE_V,PartyE_V.*MidtermE_V,...
    PresdumE_V.*XSEV_(:,1),PresdumE_V.*TenureE_V,...
    MidtermE_V.*XSEV_(:,1),MidtermE_V.*XSEV_(:,2),...
    regXS.*LOGW_NXT_E_V,regXS];



Regressand=Regressand0(:,:);


coef=(Regressand.'*Regressand)\(Regressand.'*Continue);
predictedvalue=Regressand*coef;
error=Continue-predictedvalue;
    sstot=(sum((Continue.'-mean(Continue)).^2));
                     rsq=1-sum(error.^2)/sstot ;
 % scatter(error,Continue(Continue>0))
 

%Chop off negative cont payoff AFTER having fit the regression
Continuation1(Continuation1<0)=0;


regXS_VCT=[XS_EVCT_(:,1),XS_EVCT_(:,2),XS_EVCT_(:,1).^2,XS_EVCT_(:,2).^2,log(XS_EVCT_(:,1)),XS_EVCT_(:,1).*XS_EVCT_(:,2),...
    XS_EVCT_(:,1).^2.*SameE_VCT,XS_EVCT_(:,2).^2.*PartyE_VCT,SameE_VCT,PartyE_VCT];
%Use estimated derivative in computing FOC.
Deriv=[zeros(length(Continuation1),1),ones(length(Continuation1),1),2*LOGW_NXT_E_VCT/10,3*LOGW_NXT_E_VCT.^2/100,...
    zeros(length(Continuation1),1),XQEVCT2,zeros(length(Continuation1),1),TenureE_VCT,zeros(length(Continuation1),1),...
    PresdumE_VCT,zeros(length(Continuation1),1),...
    XS_EVCT_(:,1).*SameE_VCT,XS_EVCT_(:,2).*PartyE_VCT,zeros(length(Continuation1),22),regXS_VCT,zeros(length(Continuation1),nregXS)]*coef;
    
regXS_VNCT=[XS_EVNCT_(:,1),XS_EVNCT_(:,2),XS_EVNCT_(:,1).^2,XS_EVNCT_(:,2).^2,log(XS_EVNCT_(:,1)),XS_EVNCT_(:,1).*XS_EVNCT_(:,2),...
    XS_EVNCT_(:,1).^2.*SameE_VNCT,XS_EVNCT_(:,2).^2.*PartyE_VNCT,SameE_VNCT,PartyE_VNCT];

DerivNCT=[zeros(length(Continuation2),1),ones(length(Continuation2),1),2*LOGW_NXT_E_VNCT/10,3*LOGW_NXT_E_VNCT.^2/100,...
    zeros(length(Continuation2),1),XQEVNCT2,zeros(length(Continuation2),1),TenureE_VNCT,zeros(length(Continuation2),1),...
    PresdumE_VNCT,zeros(length(Continuation2),1),...
    XS_EVNCT_(:,1).*SameE_VNCT,XS_EVNCT_(:,2).*PartyE_VNCT,zeros(length(Continuation2),22),regXS_VNCT,zeros(length(Continuation2),nregXS)]*coef;

SRR17=mean(Deriv.^2.*(Deriv<0));
std17=max(10^(-8),std(Deriv.^2.*(Deriv<0)));
SRR18=mean(Continue.^2.*(Continue<0));
std18=max(10^(-8),std(Continue.^2.*(Continue<0)));

SRR19=(coef(6)<0).*coef(6).^2;
SRR20=(coef(2)<0).*coef(2).^2;

Deriv=max(Deriv,0.00001);
DerivNCT=max(DerivNCT,0.00001);

expratioOUT=exp(LOGW_NXT_E_VCT)./exp(LOGTotal_E_VCT(:,1));
OUT=((beta*cost1*(1/beta)*(ones(length(XQEVCT2),1)+a*(1./exp(XQEVCT2))).*LOGTotal_E_VCT.^(beta-1))).*expratioOUT./(vdelta*Deriv);

%DE=exp(LOGTotal_E_VCT(:,1)).*((vdelta*Deriv).*(1./exp(LOGW_NXT_E_VCT)));

SRR13=mean((OUT-(probwin.*(probwin<1)+0.99.*(probwin>=1)))).^2; %Averaged ex-post = interim (belief of the candidate)
std13=max(10^(-8),std((OUT-(probwin.*(probwin<1)+0.99.*(probwin>=1))).^2));
SRR14=mean(OUT>1);
std14=max(10^(-8),std(OUT>1));
SRR14B=mean((OUT>1).*(OUT-1).^2);
std14B=max(10^(-8),std((OUT>1).*(OUT-1).^2));

SRR15=mean(OUT<min(probwin));
std15=max(10^(-8),std(OUT<min(probwin)));

%Pen=find(OUT<1&OUT>0.0002);
Pen=1:length(OUT);
%Pen2(IND6CT,:)=[];
%DEL_Pen=[find(Pen2<quantile(Pen2,.025));find(Pen2>quantile(Pen2,.975))];    %Cannot invert some observations:
%DEL2=find(Pen2>0);
%Pen2(DEL_Pen,:)=[];

%ako=1/2+(1/pi)*atan((OUT-0.5)*h);
ako=min(max(OUT,0.000001),0.999999);
BX1=norminv(ako);
% BX1=norminv(max(0.01,min(0.99,(beta*(alpha/beta)*ben2*...
%     ((max(0,NCE_VCT(:,1)).^(alpha-1))./(max(0,LOGTot_NCE_VCT(:,1)).^(beta-1))).*LOGTotal_E_VCT.^(beta-1))./(vdelta*Deriv))));  %%BX1=(1/sig)*(-0.5+B_I*d_I-B_C*d_C+...)=(K) in the paper.
BX1sig=BX1;
%BX1sig(IND6CT,:)=[];


%Obtain q_C_E_VCT.
q_C_E_VCT=(-1)*(sig*BX1-cc-B_I*LOGD_E_VCT-B_C*LOGD_E_VC-[XS_EVCT_(:,1),XS_EVCT_(:,1).^2,XS_EVCT_(:,1).*SameE_VCT,XS_EVCT_(:,1).^2.*SameE_VCT,XS_EVCT_(:,2).*PartyE_VCT]*thetaS2-B_T*log(TenureE_VCT+1)-XQEVCT2);
const=const0;





VV=VSEVCT-0.5-const-B_I*LOGD_E_VCT-B_C*LOGD_E_VC-[XS_EVCT_(:,1),XS_EVCT_(:,1).^2,XS_EVCT_(:,1).*SameE_VCT,XS_EVCT_(:,1).^2.*SameE_VCT,XS_EVCT_(:,2).*PartyE_VCT]*thetaS2-XQEVCT2+q_C_E_VCT-B_T*log(TenureE_VCT+1);
%VV(IND6CT,:)=[];
%VV(DEL_Pen,:)=[];
LOGD_E_VCT_m=LOGD_E_VCT;
%LOGD_E_VCT_m(IND6CT,:)=[];
%LOGD_E_VCT_m(DEL_Pen,:)=[];
LOGD_E_VC_m=LOGD_E_VC;
%LOGD_E_VC_m(IND6CT,:)=[];
%LOGD_E_VC_m(DEL_Pen,:)=[];
XS_m=[XS_EVCT_(:,1),XS_EVCT_(:,1).^2,XS_EVCT_(:,1).*SameE_VCT,XS_EVCT_(:,1).^2.*SameE_VCT,XS_EVCT_(:,2).*PartyE_VCT];
%XS_m(IND6CT,:)=[];
%XS_m(DEL_Pen,:)=[];
XQEVCT_m=XQEVCT2;
%XQEVCT_m(IND6CT,:)=[];
%XQEVCT_m(DEL_Pen,:)=[];
q_C_E_VCT_m=q_C_E_VCT;
%q_C_E_VCT_m(IND6CT,:)=[];
%q_C_E_VCT_m(DEL_Pen,:)=[];
TenureE_VCT_m=TenureE_VCT;
%TenureE_VCT_m(IND6CT,:)=[];
%TenureE_VCT_m(DEL_Pen,:)=[];

SRR11_1=mean(VV(Pen))^2;     %To estimate sigma.
std11_1=std(VV(Pen));
SRR11_2=mean(VV(Pen).*LOGD_E_VCT_m(Pen))^2;
std11_2=std(VV(Pen).*LOGD_E_VCT_m(Pen));
SRR11_3=mean(VV(Pen).*LOGD_E_VC_m(Pen))^2;
std11_3=std(VV(Pen).*LOGD_E_VC_m(Pen));
SRR11_4=mean(VV(Pen).*XS_m((Pen),1))^2;
std11_4=std(VV(Pen).*XS_m((Pen),1));
SRR11_5=mean(VV(Pen).*XS_m((Pen),2))^2;
std11_5=std(VV(Pen).*XS_m((Pen),2));
SRR11_6=mean(VV(Pen).*XQEVCT_m(Pen))^2;
std11_6=std(VV(Pen).*XQEVCT_m(Pen));
SRR11_7=mean(VV(Pen).*q_C_E_VCT_m(Pen))^2;
std11_7=std(VV(Pen).*q_C_E_VCT_m(Pen));
SRR11_8=mean(VV(Pen).*log(TenureE_VCT_m(Pen)+1))^2;
std11_8=std(VV(Pen).*log(TenureE_VCT_m(Pen)+1));
SRR11_9=mean(VV(Pen).*XS_m((Pen),3))^2;
std11_9=std(VV(Pen).*XS_m((Pen),3));
SRR11_10=mean(VV(Pen).*XS_m((Pen),4))^2;
std11_10=std(VV(Pen).*XS_m((Pen),4));
SRR11_11=mean(VV(Pen).*XS_m((Pen),5))^2;
std11_11=std(VV(Pen).*XS_m((Pen),5));
SRR11=SRR11_1/std11_1+SRR11_2/std11_2+SRR11_3/std11_3+SRR11_4/std11_4...
   +SRR11_5/std11_5+SRR11_6/std11_6+SRR11_7/std11_7+SRR11_8/std11_8+SRR11_9/std11_9...
   +SRR11_10/std11_10+SRR11_11/std11_11;


%Extra set of moments: in expectation, q_C_E_VCT without constant term has
%to be on average equal to the average of first stage E(qC|X) without constant term.
qcvv=q_C_E_VCT-expectedqcCT; %corresponding to eta.


SRR31_1=mean(qcvv(Pen))^2;     %To estimate sigma.
std31_1=max(0.01,std(qcvv(Pen)));
% SRR31_2=mean(qcvv(Pen).*LOGD_E_VCT_m(Pen))^2;
% std31_2=max(0.01,std(qcvv(Pen).*LOGD_E_VCT_m(Pen)));
% SRR31_3=mean(qcvv(Pen).*LOGD_E_VC_m(Pen))^2;
% std31_3=max(0.01,std(qcvv(Pen).*LOGD_E_VC_m(Pen)));
SRR31_4=mean(qcvv(Pen).*XS_m((Pen),1))^2;
std31_4=max(0.01,std(qcvv(Pen).*XS_m((Pen),1)));
SRR31_5=mean(qcvv(Pen).*XS_m((Pen),2))^2;
std31_5=max(0.01,std(qcvv(Pen).*XS_m((Pen),2)));
SRR31_6=mean(qcvv(Pen).*XQEVCT_m(Pen))^2;
std31_6=max(0.01,std(qcvv(Pen).*XQEVCT_m(Pen)));
% SRR31_7=mean(qcvv(Pen).*q_C_E_VCT_m(Pen))^2;
% std31_7=max(0.01,std(qcvv(Pen).*q_C_E_VCT_m(Pen)));
SRR31_8=mean(qcvv(Pen).*log(TenureE_VCT_m(Pen)+1))^2;
std31_8=max(0.01,std(qcvv(Pen).*log(TenureE_VCT_m(Pen)+1)));
SRR31_9=mean(qcvv(Pen).*XS_m((Pen),3))^2;
std31_9=max(0.01,std(qcvv(Pen).*XS_m((Pen),3)));
SRR31_10=mean(qcvv(Pen).*XS_m((Pen),4))^2;
std31_10=max(0.01,std(qcvv(Pen).*XS_m((Pen),4)));
SRR31_11=mean(qcvv(Pen).*XS_m((Pen),5))^2;
std31_11=max(0.01,std(qcvv(Pen).*XS_m((Pen),5)));
SRR31=SRR31_1/std31_1+SRR31_4/std31_4...
    +SRR31_5/std31_5+SRR31_6/std31_6+SRR31_8/std31_8+SRR31_9/std31_9...
    +SRR31_10/std31_10+SRR31_11/std31_11;

% SRR31=SRR31_1+SRR31_2+SRR31_3+SRR31_4...
%     +SRR31_5+SRR31_6+SRR31_7+SRR31_8+SRR31_9...
%     +SRR31_10+SRR31_11;

SRR32=(residvar-(mean(qcvv.^2)+sig^2))^2;

SRR9=-log(sig)-(1/length(VV(Pen)))*sum((VV(Pen).^2)/(2*sig^2));
SRR9=-SRR9;
std9=sum((-log(sig)-(VV(Pen).^2)/(2*sig^2)).^2)/length(VV(Pen));

SRR9P=mean((VV(Pen)-mean(VV(Pen))).^2-sig^2).^2;
%SRR9P_1=mean(((VV-mean(VV)).^2-sig^2))^2;
std9P=std(((VV(Pen)-mean(VV(Pen))).^2-sig^2));
%SRR9P_2=mean(sqrt(VV.^2)-sig)^2;
%std9P_2=std(sqrt(VV.^2)-sig);
% foc of incumbent, contested periods %%
%FOC11=(beta*cost1*(alpha/beta)*ben2*(exp(LOGTot_NCE_VCT(:,1))./exp(NCE_VCT(:,1))).*...
%    ((max(0,NCE_VCT(:,1)).^(alpha-1))./(max(0,LOGTot_NCE_VCT(:,1)).^(beta-1))).*LOGTotal_E_VCT.^(beta-1)).*(1./exp(LOGTotal_E_VCT(:,1)))...      %% d/dI C_I(total)
%    -(B_I/(sig*exp(LOGD_E_VCT(:,1))))*normpdf(BX1).*(1+vdelta*Continuation1)-alpha*ben1*(1./exp(LOGD_E_VCT(:,1))).*(max(0,(LOGD_E_VCT))).^(alpha-1); %% d/d_I(P_2)(B+delta*EV+(d/dI)H_I())
expratio11=exp(LOGD_E_VCT(:,1))./exp(LOGTotal_E_VCT(:,1));
FOC11=(beta*cost1*(1/beta)*(ones(length(XQEVCT2),1)+a*(1./exp(XQEVCT2))).*LOGTotal_E_VCT.^(beta-1)).*expratio11...      %% d/dI C_I(total)
    -(B_I/sig).*normpdf(BX1).*(BB+vdelta*Continuation1)-alpha*ben1.*(max(0,(LOGD_E_VCT))).^(alpha-1); %% d/d_I(P_2)(B+delta*EV+(d/dI)H_I())



%FOC11(IND6CT,:)=[];
%FOC11(DEL_Pen,:)=[];
%FOC11(DEL2,:)=[];
SRR10P=mean(FOC11(Pen).^2,1);
std10P=std(FOC11(Pen).^2);
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
expratio21=exp(LOGTotal_E_VNCT(:,1))./exp(LOGW_NXT_E_VNCT);
FOC21=vdelta*DerivNCT.*expratio21-... 
(beta*cost2*(1/beta)*(ones(length(XQEVNCT2),1)+a*(1./exp(XQEVNCT2))).*LOGTotal_E_VNCT.^(beta-1));   %% delta*(d/dw_I)E_V-C_I'
%FOC21(IND6NCT,:)=[];

expratio22=exp(LOGTotal_E_VNCT(:,1))./exp(LOGD_E_VNCT(:,1));
FOC22=alpha*ben2.*(max(0,(LOGD_E_VNCT))).^(alpha-1).*expratio22-... 
(beta*cost2*(1/beta)*(ones(length(XQEVNCT2),1)+a*(1./exp(XQEVNCT2))).*LOGTotal_E_VNCT.^(beta-1));   %% delta*(d/dw_I)E_V-C_I'
%FOC21(IND6NCT,:)=[];

% FOC21delta=0.9*(DContinuation2b2-Continuation2b2)-(ben2+incr)*RTotDE_VNCT.*((LOGTotal_E_VNCT+Delt).^2-LOGTotal_E_VNCT.^2);   %% delta*(d/dw_I)E_V-C_I'
% FOC21delta(IND6NCT,:)=[];
% SRR12=((1/incr)*((1/290)*sum(FOC21delta.^2,1)-(1/290)*sum(FOC21.^2,1)))^2;
SRR12P=mean(FOC21.^2,1);
std12P=std(FOC21.^2);

SRR16P=mean(FOC22.^2,1);
std16P=std(FOC22.^2);
% SRR9
% SRR10P
% % SRR11P
% SRR12P


%ADD CHALLENGER SIDE
LOGTot_E_VC=datasetVC(:,1);
LOGW_NXT_E_VC=datasetVC(:,2);
LOGD_E_VC=datasetVC(:,3);


%Compute challenger's continuation payoff using the estimated coefficients

%Compute continue, then find the minimum, pick q_C_E_VCT<that minimum,
%equalize the value, and then take max(that value, 0).
regXS_VC=[XS_EVCT_(:,1),XS_EVCT_(:,2),XS_EVCT_(:,1).^2,XS_EVCT_(:,2).^2,log(XS_EVCT_(:,1)),XS_EVCT_(:,1).*XS_EVCT_(:,2),...
    XS_EVCT_(:,1).^2.*(-1*SameE_VCT),XS_EVCT_(:,2).^2.*(-1*PartyE_VCT),(-1*SameE_VCT),(-1*PartyE_VCT)];

Continuec=[ones(length(LOGW_NXT_E_VC),1),LOGW_NXT_E_VC,LOGW_NXT_E_VC.^2/10,LOGW_NXT_E_VC.^3/100,...
    q_C_E_VCT,LOGW_NXT_E_VC.*q_C_E_VCT,...
    exp(q_C_E_VCT),...
    zeros(length(LOGW_NXT_E_VC),2),...
    LOGW_NXT_E_VC.*PresdumE_VCT,q_C_E_VCT.*PresdumE_VCT,...
    LOGW_NXT_E_VC.*XS_EVCT_(:,1).*(-1*SameE_VCT),LOGW_NXT_E_VC.*XS_EVCT_(:,2).*(-1*PartyE_VCT),...   
    q_C_E_VCT.*XS_EVCT_(:,1).*(-1*SameE_VCT),...
    zeros(length(LOGW_NXT_E_VC),2),XS_EVCT_(:,1).*(-1*SameE_VCT),XS_EVCT_(:,1).^2.*(-1*SameE_VCT),XS_EVCT_(:,2).*(-1*PartyE_VCT),XS_EVCT_(:,2).^2.*(-1*PartyE_VCT),...
    PresdumE_VCT,MidtermE_VCT,PresdumE_VCT.*MidtermE_VCT,...
    (-1*SameE_VCT).*XS_EVCT_(:,2),(-1*SameE_VCT).*(-1*PartyE_VCT),(-1*SameE_VCT).*TenureE_VCT,(-1*SameE_VCT).*PresdumE_VCT,(-1*SameE_VCT).*MidtermE_VCT,...
    (-1*PartyE_VCT).*XS_EVCT_(:,1),(-1*PartyE_VCT).*TenureE_VCT,(-1*PartyE_VCT).*MidtermE_VCT,...
    PresdumE_VCT.*XS_EVCT_(:,1),PresdumE_VCT.*TenureE_VCT,...
    MidtermE_VCT.*XS_EVCT_(:,1),MidtermE_VCT.*XS_EVCT_(:,2),...
    regXS_VC.*LOGW_NXT_E_VC,regXS_VC]*coef;    

% %Benchmark-minimum utility type
% aux0=q_C_E_VCT(Continuec==min(Continuec));
% 
% %Find all people such that either (1)theta<benchmark in all dimensions
% %or alternatively, (2) mean(theta)<mean(benchmark).
% %dropind1=initbelief(1,:)<aux0(1)&initbelief(2,:)<aux0(2)&initbelief(3,:)<aux0(3)...
% %    &initbelief(4,:)<aux0(4)&initbelief(5,:)<aux0(5);
% dropind2=q_C_E_VCT<aux0;
% 
% Continuec(dropind2)=min(Continuec);
Continuec(Continuec<0)=0;


Derivc=[zeros(length(LOGW_NXT_E_VC),1),ones(length(LOGW_NXT_E_VC),1),2*LOGW_NXT_E_VC/10,3*LOGW_NXT_E_VC.^2/100,...
    zeros(length(LOGW_NXT_E_VC),1),q_C_E_VCT,zeros(length(LOGW_NXT_E_VC),1),zeros(length(LOGW_NXT_E_VC),1),zeros(length(LOGW_NXT_E_VC),1),...
    PresdumE_VCT,zeros(length(LOGW_NXT_E_VC),1),...
    XS_EVCT_(:,1).*(-1*SameE_VCT),XS_EVCT_(:,2).*(-1*PartyE_VCT),zeros(length(LOGW_NXT_E_VC),22),regXS_VC,zeros(length(LOGW_NXT_E_VC),nregXS)]*coef;
    
Derivc(Continuec==min(min(Continuec),0))=0;

%Alternative, assuming quadratic cont payoff
%If quality at the range where continuation payoff is decreasing in
%quality, assign zero.
% threshqual=log(-1/coef(7)*(coef(5)+LOGW_NXT_E_VC*coef(6)+PresdumE_VCT*coef(11)+XS_EVCT_(:,1).*(-1*SameE_VCT)*coef(14)));
% 
% Continuec=(threshqual<q_C_E_VCT).*aux;
% 
% aux2=[zeros(length(LOGW_NXT_E_VC),1),ones(length(LOGW_NXT_E_VC),1),2*LOGW_NXT_E_VCT/10,3*LOGW_NXT_E_VCT.^2/100,...
%     zeros(length(LOGW_NXT_E_VC),1),q_C_E_VCT,zeros(length(LOGW_NXT_E_VC),1),zeros(length(LOGW_NXT_E_VC),1),zeros(length(LOGW_NXT_E_VC),1),...
%     PresdumE_VCT,zeros(length(LOGW_NXT_E_VC),1),...
%     XS_EVCT_(:,1).*(-1*SameE_VCT),XS_EVCT_(:,2).*(-1*PartyE_VCT),zeros(length(LOGW_NXT_E_VC),10)]*coef;
%     
% Derivc=(threshqual<q_C_E_VCT).*aux2;


BX1c=abs(1/sig)*(-1)*([LOGD_E_VCT,LOGD_E_VC,[XS_EVCT_(:,1),XS_EVCT_(:,1).^2,XS_EVCT_(:,1).*SameE_VCT,XS_EVCT_(:,1).^2.*SameE_VCT,XS_EVCT_(:,2).*PartyE_VCT],log(TenureE_VCT+1)]*Est2(2:9)+XQEVCT2-q_C_E_VCT+const);

expratio31=exp(LOGD_E_VC)./exp(LOGTot_E_VC);
FOC31=costc*betac*(1/betac)...
   .*(ones(length(q_C_E_VCT),1)+ac*(1./exp(q_C_E_VCT))).*(LOGTot_E_VC).^(betac-1).*expratio31...      %% d/dI C_I(total)
   +(B_C/sig).*normpdf(BX1c).*(BB+vdelta*Continuec)-alphac*benc*(LOGD_E_VC.^(alphac-1)); %% d/d_I(P_2)(B+delta*EV+(d/dI)H_I())

expratio32=exp(LOGW_NXT_E_VC)./exp(LOGTot_E_VC);
FOC32=costc*betac*(1/betac)...
   .*(ones(length(q_C_E_VCT),1)+ac*(1./exp(q_C_E_VCT))).*(LOGTot_E_VC).^(betac-1).*expratio32...  
   -vdelta*normcdf(BX1c).*Derivc;

FOC321=FOC32.*(LOGW_NXT_E_VC>0.1);
FOC322=min(0,FOC32).*(LOGW_NXT_E_VC<0.1);

FOC32=FOC321+FOC322;

SRR21P=mean(FOC31(Pen).^2,1);
std21P=std(FOC31(Pen).^2);

SRR23=mean(FOC32(Pen).^2,1);
std23=std(FOC32(Pen).^2);

%sp1 - preferred specification
% SRR2step=SRR9P/std9P+SRR10P+SRR11+SRR12P+SRR16P+10*SRR14+SRR13/std13+SRR21P+SRR23+SRR18...%SRR17/std17+SRR18/std18+SRR19...
%    +SRR31+10000*SRR32+100*(a<0)*a^2+100*(ac<0)*(ac)^2;%+100*(beta<1)*(1-beta)^2;%+SRR20*100;%+10000*(sigma>residstd);%+;%;%

%sp2
SRR2step=SRR9P/std9P+SRR10P+... %Incumbent FOC contested
    SRR11+SRR31+... %Orthogonality between vote share errors and observables
    SRR12P+SRR16P+... %Incumbent FOC uncontested
    10*SRR14+SRR13/std13+... %OUT<1
    SRR21P+SRR23+... %Challenger FOC
    100*SRR18+... %Continuation payoff>0
   10000*SRR32+... %Variance of errors match with first stage
   100*(a<0)*a^2+100*(ac<0)*(ac)^2; %Cost decreasing in quality


if isnan(SRR2step)==1
    SRR2step=9999;
end
%%

% regX=[ones(length(LOGD_E_VCT),1),LOGD_E_VCT,LOGD_E_VC,[XS_EVCT_(:,1),XS_EVCT_(:,1).^2,XS_EVCT_(:,2).*PartyE_VCT],log(TenureE_VCT+1),q_C_E_VCT,XQEVCT2];
% [b,~,r]=regress(VSEVCT,regX);
% rr=sqrt(sum(r.^2/403));
rr=sig;
%%
  save q_C_E_VCT q_C_E_VCT
  save XQEV2 XQEV2
  save est2ndstage mintheta2ndstage
 % save DEL_Pen DEL_Pen
  save Continue Continue
  save OUT OUT
  coef2ndstep=coef;
  save coef2ndstep coef2ndstep
  save const const
  save rr rr

mintheta2ndstage(3)=rr;
  save est2ndstage mintheta2ndstage
  coef2ndstep=coef;
  save coef2ndstep coef2ndstep
  %%
%Outsheet the csv for creating quality estimation data
q_C_E_V=zeros(length(XQEV2),1);
q_C_E_V(E_VNContestFUL)=q_C_E_VCT;
computeincadv=[E_V_july8(:,[1:4,9,10,16,13,24,63]),q_C_E_V,XQEV2];
dlmwrite('incquality.csv',computeincadv,'delimiter', ',', 'precision', 16)
%appearence	year	party	opponent_total	realtotal	realdisburse	opponent_disburse	voteshare	unemploymentrate	partindex	qc	qi




%%
%%%%%%%%%%%
%Specification check
%%%%%%%%%%


%Cost function
A=@(test,test2)cost1*(1/beta)*(1+a*(1./exp(test))).*test2.^beta;
Ac=@(test,test2)costc*(1/betac)*(1+a*(1./exp(test))).*test2.^betac;
space1=-0.2:0.001:0.2;
space2=10:0.5:17;
mat1=A(space1,11.0);
mat1c=Ac(space1,11.0);
mat2=A(0.0382,space2);
mat2c=Ac(0.0382,space2);
figure(1)
scatter(space1,mat1)
%hold on
%scatter(space2,mat2c)
%scatter(space2,mat2)
%%
%Cost derivative
Aderiv=@(test,test2)beta*cost1*(1/beta)*(1+a*(1./exp(test))).*test2.^(beta-1)./exp(test2);
space1=0.01:0.001:0.0382;
space2=9:0.5:15;
mat1=Aderiv(space1,9);
mat2=Aderiv(0.03,space2);
scatter(space1,mat1)
%scatter(space2,mat2)
%%
%Derivative w.r.t LOGW
%Flex

B=@(test,test2)coef(2)+2*test2/10*coef(3)+3*test2.^2/100*coef(4)+test*coef(6)+mean(TenureE_V)*coef(8)...
   +mean(PresdumE_V)*coef(10)+mean(XSEV_(:,1).*SameE_V)*coef(12)+mean(XSEV_(:,2).*PartyE_V)*coef(13)...
   +mean(regXS)*coef(36:45);

space1=-0.01:0.001:0.01;
space2=9:0.5:15;
mat1=max(0,B(space1,12));
mat2=B(0.01,space2);
scatter(space1,mat1)
%%
%Cont payoff itself
% 
% contpayofffunc=@(test,test2) coef(1)+test2*coef(2)+test2.^2/10*coef(3)+test2.^3/100*coef(4)+test*coef(5)+test*test2*coef(6)+exp(test)*coef(7)...
%    +mean(TenureE_V)*test2*coef(8)+mean(TenureE_V)*test*coef(9)+mean(PresdumE_V)*test2*coef(10)+mean(PresdumE_V)*test*coef(11)...
%    +mean(XSEV_(:,1).*SameE_V)*test2*coef(12)+mean(XSEV_(:,2).*PartyE_V)*test2*coef(13)...
%    +mean(XSEV_(:,1).*SameE_V)*test*coef(14)+mean([TenureE_V,TenureE_V.^2,XSEV_(:,1).*SameE_V,XSEV_(:,1).^2.*SameE_V,XSEV_(:,2).*PartyE_V,XSEV_(:,2).^2.*PartyE_V,...
%     PresdumE_V,MidtermE_V,PresdumE_V.*MidtermE_V])*coef(15:23);

contpayofffunc=@(test,test2) coef(1)+test2*coef(2)+test2.^2/10*coef(3)+test2.^3/100*coef(4)+test*coef(5)+test*test2*coef(6)+exp(test)*coef(7)...
   +0*test2*coef(8)+0*test*coef(9)+mean(PresdumE_V)*test2*coef(10)+mean(PresdumE_V)*test*coef(11)...
   +mean(XSEV_(:,1).*SameE_V)*test2*coef(12)+mean(XSEV_(:,2).*PartyE_V)*test2*coef(13)...
   +mean(XSEV_(:,1).*SameE_V)*test*coef(14)+mean([XSEV_(:,1).*SameE_V,XSEV_(:,1).^2.*SameE_V,XSEV_(:,2).*PartyE_V,XSEV_(:,2).^2.*PartyE_V,...
    PresdumE_V,MidtermE_V,PresdumE_V.*MidtermE_V,...
    SameE_V.*XSEV_(:,2),SameE_V.*PartyE_V,SameE_V.*TenureE_V,SameE_V.*PresdumE_V,SameE_V.*MidtermE_V,...
    PartyE_V.*XSEV_(:,1),PartyE_V.*TenureE_V,PartyE_V.*MidtermE_V,...
    PresdumE_V.*XSEV_(:,1),PresdumE_V.*TenureE_V,...
    MidtermE_V.*XSEV_(:,1),MidtermE_V.*XSEV_(:,2)])*coef(17:35);




space1=0:0.001:0.05;
space2=0:0.5:15;
mat1=contpayofffunc(space1,13);
mat2=contpayofffunc(0,space2);
scatter(space2,mat2)

%%
%H function
C=@(test)ben1*(max(0,(test))).^(alpha); 
CC=@(test)ben1*(max(0,(test))).^(alphac); 
space1=0:0.5:15;
mat1=C(space1);
mat2=CC(space1);
scatter(space1,mat1)
hold on
scatter(space1,mat2)
%%
%Derivative of H function
CCC=@(test)alpha*ben1*(1./exp(test)).*(max(0,(test))).^(alpha-1); %% d/

space1=9:0.5:15;
mat1=CCC(space1);
scatter(space1,mat1)

%%

%h=1;
cost1=abs(theta2(1,1));  %% coefficient on the cost function of incumbent, contested
ben1=abs(theta2(2,1));   %% coefficient on the benefit function of incumbent, contested
ben2=ben1;
benc=ben1;

% alpha=abs(theta2(5,1));
% beta=abs(theta2(6,1));
alpha=0.5;
beta=2;

sig=abs(theta2(3,1));
a=abs(theta2(4,1));
costc=cost1;%abs(theta2(6,1));
ac=abs(theta2(4,1));



%Cost function
costchal=@(test,test2)costc*(1/beta)*(1+ac*(1./exp(test))).*test2.^beta;
costinc=@(test,test2)cost1*(1/beta)*(1+a*(1./exp(test))).*test2.^beta;

%Specification 1
% space1=0.001:0.01:0.201;
% matchal1=costchal(space1,10);
% matinc1=costinc(space1,10);
%Specification 2
space1=0:0.001:0.038;
space3=-0.2:0.01:0.2;
matchal1=costchal(space3,12);
matinc1=costinc(space1,12);

scatter(space3,matchal1)
hold on
scatter(space1,matinc1);

%%

hist(OUT,20)
%%
hist(OUT(OUT<1),20)
%%
hist(OUT(q_C_E_VCT>0.1),20)
%%
hist(DE(q_C_E_VCT>0.1),20)

%%
size(OUT(OUT>1))
%%
hist(XQEV2,20)
%%
hist(q_C_E_VCT,20)
%%
hist(q_C_E_VCT(OUT>1),20)
%%
hist(Continue,20)
%%
scatter(Continue,XQEV2)
%%
hist(Derivhist1(Derivhist1>-0.2&Derivhist1<0.2),50)
size(Derivhist1(Derivhist1<0),1)+size(Derivhist2(Derivhist2<0),1)
%%
hist(Derivhist2,50)
%%

bins=linspace(-0.5,0.5,100);
[hist1,scale1]=hist(q_C_E_VCT,bins);
[hist2,scale2]=hist(XQEV2,bins);
bar(hist2,'c')
hold on;
bar(hist1,'EdgeColor','r','FaceColor','none')

%%

figure(1)
hist(LOGW_NXT_E_VCT(Derivhist1<0),20)
figure(2)
hist(LOGW_NXT_E_VCT,20)
%%

figure(1)
hist(LOGW_E_VCT(Derivhist1<0),20)
figure(2)
hist(LOGW_E_VCT,20)
%%
figure(1)
hist(LOGD_E_VC,20)
figure(2)
hist(LOGD_E_VC(Derivhist1<0),20)
%%
figure(1)
hist(LOGD_E_VCT,20)
figure(2)
hist(LOGD_E_VCT(Derivhist1<0),20)
%%
figure(1)
hist(LOGTotal_E_VCT,20)
figure(2)
hist(LOGTotal_E_VCT(Derivhist1<0),20)