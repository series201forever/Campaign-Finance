function [SRR2step,rsq]=Minimize2S_new(Est2,theta2,probwin,datasetV,datasetVCT,datasetVNCT,datasetVC,deltamat,const0,NN,C1,C2,C3,C5,expectedqcCT,residvar,T)
%%
% Est2 : the estimated first step estimator.
% theta2 : the parameters that we estimate in the 2nd step.

%Calculation of continuation payoff
NCE_V=datasetV(:,1:3);
LOGTot_NCE_V=datasetV(:,4);
E_VContestFUL=datasetVNCT(:,1);
E_VNContestFUL=datasetVCT(:,1);

%Calculation of derivative of payoff
XQEV2=datasetV(:,5);
LOGW_NXT_E_V=datasetV(:,6);
SameE_V=datasetV(:,7);
PartyE_V=datasetV(:,8);
TenureE_V=datasetV(:,9);
XSEV_=datasetV(:,10:11);
PresdumE_V=datasetV(:,12);
MidtermE_V=datasetV(:,13);




%Derivation of challenger quality and incumbent FOC contested
XQEVCT2=datasetVCT(:,2);
LOGTot_NCE_VCT=datasetVCT(:,3);
NCE_VCT=datasetVCT(:,4:6);
LOGW_NXT_E_VCT=datasetVCT(:,7);
LOGTotal_E_VCT=datasetVCT(:,8);

VSEVCT=datasetVCT(:,9);
XS_EVCT_=datasetVCT(:,10:11);
TenureE_VCT=datasetVCT(:,12);
PartyE_VCT=datasetVCT(:,13);
SameE_VCT=datasetVCT(:,14);
LOGD_E_VCT=datasetVCT(:,15);
LOGD_E_VC=datasetVCT(:,16);
PresdumE_VCT=datasetVCT(:,17);
MidtermE_VCT=datasetVCT(:,18);
matchXQEVCT2=datasetVCT(:,19);

%Derication of incumbent FOC uncontested
LOGW_NXT_E_VNCT=datasetVNCT(:,2);
LOGTotal_E_VNCT=datasetVNCT(:,3);
NCE_VNCT=datasetVNCT(:,4:6);
LOGTot_NCE_VNCT=datasetVNCT(:,7);
XQEVNCT2=datasetVNCT(:,8);
LOGD_E_VNCT=datasetVNCT(:,9);
TenureE_VNCT=datasetVNCT(:,10);
XS_EVNCT_=datasetVNCT(:,11:12);
PartyE_VNCT=datasetVNCT(:,13);
SameE_VNCT=datasetVNCT(:,14);
PresdumE_VNCT=datasetVNCT(:,15);
MidtermE_VNCT=datasetVNCT(:,16);



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

%sp1
SRR2step=SRR9P/std9P+SRR10P+... %Incumbent FOC contested
    SRR11+SRR31+... %Orthogonality between vote share errors and observables
    SRR12P+SRR16P+... %Incumbent FOC uncontested
    10*SRR14+SRR13/std13+... %OUT<1
    SRR21P+SRR23+... %Challenger FOC
    100*SRR18+... %Continuation payoff>0
   10000*SRR32+... %Variance of errors match with first stage
   100*(a<0)*a^2+100*(ac<0)*(ac)^2; %Cost decreasing in quality

%sp1
% SRR2step=SRR9P/std9P+SRR10P/std10P+... %Incumbent FOC contested
%     SRR11+SRR31+... %Orthogonality between vote share errors and observables
%     SRR12P/std12P+SRR16P/std16P+... %Incumbent FOC uncontested
%     SRR21P/std21P+SRR23/std23+... %Challenger FOC
%     10*SRR14+SRR13/std13+... %OUT<1
%     100*SRR18+... %Continuation payoff>0
%    10000*SRR32+... %Variance of errors match with first stage
%    100*(a<0)*a^2+100*(ac<0)*(ac)^2; %Cost decreasing in quality
% 

if isnan(SRR2step)==1
    SRR2step=9999;
end
end

