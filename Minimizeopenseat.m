function SRR4step=Minimizeopenseat(mintheta2ndstage,theta4,Est2,datasetVCT,coef2ndstep)

%% thetain : the estimated first step estimator.
%% theta2 : the parameters that we estimate in the 2nd step.


%Calculation of derivative of payoff
XQEV=datasetVCT(:,1);
LOGW_NXT_E_V=datasetVCT(:,2);
LOGTotal_E_V=datasetVCT(:,3);
LOGD_E_V=datasetVCT(:,4);
XSEV_=datasetVCT(:,5:6);
SameE_V=datasetVCT(:,7);
PartyE_V=datasetVCT(:,8);
PresdumE_V=datasetVCT(:,9);
MidtermE_V=datasetVCT(:,10);
LOGW_NXT_E_VC=datasetVCT(:,11);
LOGTot_E_VC=datasetVCT(:,12);
LOGD_E_VC=datasetVCT(:,13);
VSEVCT=datasetVCT(:,14);

vdelta=0.90;
% thetaS=thetain(1:2,1);
% B_I=thetain(1,1);
% B_C=thetain(2,1);
% B_T=thetain(5,1);
% thetaS2=thetain(3:4,1);

BB=1;
alpha=abs(mintheta2ndstage(7));
beta=abs(mintheta2ndstage(9));
ben1=abs(mintheta2ndstage(2));
ben2=ben1;
sig=abs(mintheta2ndstage(3));

%costc=abs(theta4(2));
costc=abs(mintheta2ndstage(8));
a=abs(mintheta2ndstage(4));
B_I=abs(theta4(1));


thetaS2=Est2(4:8,1);


%Calculate cont payoff
regXS=[XSEV_(:,1),XSEV_(:,2),XSEV_(:,1).^2,XSEV_(:,2).^2,log(XSEV_(:,1)),XSEV_(:,1).*XSEV_(:,2),...
    XSEV_(:,1).^2.*SameE_V,XSEV_(:,2).^2.*PartyE_V,SameE_V,PartyE_V];
nregXS=size(regXS,2);

Continuei=max(0,[ones(length(LOGW_NXT_E_V),1),LOGW_NXT_E_V,LOGW_NXT_E_V.^2/10,LOGW_NXT_E_V.^3/100,...
    XQEV,LOGW_NXT_E_V.*XQEV,...
    exp(XQEV),...
    zeros(length(LOGW_NXT_E_V),2),...
    LOGW_NXT_E_V.*PresdumE_V,XQEV.*PresdumE_V,...
    LOGW_NXT_E_V.*XSEV_(:,1).*SameE_V,LOGW_NXT_E_V.*XSEV_(:,2).*PartyE_V,...
    XQEV.*XSEV_(:,1).*SameE_V,...
    zeros(length(LOGW_NXT_E_V),2),XSEV_(:,1).*SameE_V,XSEV_(:,1).^2.*SameE_V,XSEV_(:,2).*PartyE_V,XSEV_(:,2).^2.*PartyE_V,...
    PresdumE_V,MidtermE_V,PresdumE_V.*MidtermE_V,...
    SameE_V.*XSEV_(:,2),SameE_V.*PartyE_V,zeros(length(LOGW_NXT_E_V),1),SameE_V.*PresdumE_V,SameE_V.*MidtermE_V,...
    PartyE_V.*XSEV_(:,1),zeros(length(LOGW_NXT_E_V),1),PartyE_V.*MidtermE_V,...
    PresdumE_V.*XSEV_(:,1),zeros(length(LOGW_NXT_E_V),1),...
    MidtermE_V.*XSEV_(:,1),MidtermE_V.*XSEV_(:,2),regXS.*LOGW_NXT_E_V,regXS]*coef2ndstep);

Derivi=max(0.001,[zeros(length(Continuei),1),ones(length(Continuei),1),2*LOGW_NXT_E_V/10,3*LOGW_NXT_E_V.^2/100,...
    zeros(length(Continuei),1),XQEV,zeros(length(Continuei),1),zeros(length(Continuei),2),...
    PresdumE_V,zeros(length(Continuei),1),...
    XSEV_(:,1).*SameE_V,XSEV_(:,2).*PartyE_V,zeros(length(Continuei),22),regXS,zeros(length(Continuei),nregXS)]*coef2ndstep);
    
 


OUTI=((beta*costc*(1/beta)*(ones(length(XQEV),1)+a*(1./exp(XQEV))).*LOGTotal_E_V.^(beta-1)).*(1./exp(LOGTotal_E_V(:,1))))./((vdelta*Derivi).*(1./exp(LOGW_NXT_E_V)));

SRR14=(OUTI>1);
SRR15=(OUTI<0);

Pen=max(OUTI-1,0)-min(OUTI,0);
Pen2=Pen;
%Pen2(IND6CT,:)=[];
%DEL_Pen=[find(Pen2<quantile(Pen2,.025));find(Pen2>quantile(Pen2,.975))];    %Cannot invert some observations:
DEL_Pen=[find(OUTI>1);find(OUTI<0)];
%DEL2=find(Pen2>0);
%Pen2(DEL_Pen,:)=[];

%ako=1/2+(1/pi)*atan((OUT-0.5)*h);
akoi=min(max(OUTI,0.000001),0.999999);
BX1i=norminv(akoi);
% BX1=norminv(max(0.01,min(0.99,(beta*(alpha/beta)*ben2*...
%     ((max(0,NCE_VCT(:,1)).^(alpha-1))./(max(0,LOGTot_NCE_VCT(:,1)).^(beta-1))).*LOGTotal_E_VCT.^(beta-1))./(vdelta*Deriv))));  %%BX1=(1/sig)*(-0.5+B_I*d_I-B_C*d_C+...)=(K) in the paper.
%BX1sig=BX1;
%BX1sig(IND6CT,:)=[];

% foc
FOC11i=(beta*costc*(1/beta)*(ones(length(XQEV),1)+a*(1./exp(XQEV))).*LOGTotal_E_V.^(beta-1)).*(1./exp(LOGTotal_E_V(:,1)))...      %% d/dI C_I(total)
    -(B_I./(sig*exp(LOGD_E_V(:,1)))).*normpdf(BX1i).*(BB+vdelta*Continuei)-alpha*ben1*(1./exp(LOGD_E_V(:,1))).*(max(0,(LOGD_E_V))).^(alpha-1); %% d/d_I(P_2)(B+delta*EV+(d/dI)H_I())

%FOC11(IND6CT,:)=[];
FOC11i(DEL_Pen,:)=[];
%FOC11(DEL2,:)=[];
SRR10P=mean(FOC11i.^2,1);
std10P=std(FOC11i.^2);
% 
% %If directly estimate B_I only using incumbent FOC
% numel=-mean(((beta*costc*(1/beta)*(ones(length(XQEV),1)+a*(1./exp(XQEV))).*LOGTotal_E_V.^(beta-1)).*(1./exp(LOGTotal_E_V(:,1)))...
%     -alpha*ben1*(1./exp(LOGD_E_V(:,1))).*(max(0,(LOGD_E_V))).^(alpha-1))./(sig*exp(LOGD_E_V(:,1))).*(normpdf(BX1).*(1+vdelta*Continuei)));
% denom=mean(((1./(sig*exp(LOGD_E_V(:,1)))).*normpdf(BX1).*(1+vdelta*Continuei)).^2);
% BI=numel/denom;

%Obtain challenger quality
q_C_E_VCT=(-1)*(sig*BX1i-B_I*LOGD_E_V+B_I*LOGD_E_VC-[XSEV_(:,1),XSEV_(:,1).^2,XSEV_(:,1).*SameE_V,XSEV_(:,1).^2.*SameE_V,XSEV_(:,2).*PartyE_V]*thetaS2-XQEV);


%Ex-post voteshare
VV=VSEVCT-0.5-B_I*LOGD_E_V+B_I*LOGD_E_VC-[XSEV_(:,1),XSEV_(:,1).^2,XSEV_(:,1).*SameE_V,XSEV_(:,1).^2.*SameE_V,XSEV_(:,2).*PartyE_V]*thetaS2-XQEV+q_C_E_VCT;

XS=[XSEV_(:,1),XSEV_(:,1).^2,XSEV_(:,1).*SameE_V,XSEV_(:,1).^2.*SameE_V,XSEV_(:,2).*PartyE_V];

Pen=1:length(OUTI);
SRR11_1=mean(VV(Pen))^2;     %To estimate sigma.
std11_1=std(VV(Pen));
SRR11_2=mean(VV(Pen).*LOGD_E_V(Pen))^2;
std11_2=std(VV(Pen).*LOGD_E_V(Pen));
SRR11_3=mean(VV(Pen).*LOGD_E_VC(Pen))^2;
std11_3=std(VV(Pen).*LOGD_E_VC(Pen));
SRR11_4=mean(VV(Pen).*XS((Pen),1))^2;
std11_4=std(VV(Pen).*XS((Pen),1));
SRR11_5=mean(VV(Pen).*XS((Pen),2))^2;
std11_5=std(VV(Pen).*XS((Pen),2));
SRR11_6=mean(VV(Pen).*XQEV(Pen))^2;
std11_6=std(VV(Pen).*XQEV(Pen));
SRR11_7=mean(VV(Pen).*q_C_E_VCT(Pen))^2;
std11_7=std(VV(Pen).*q_C_E_VCT(Pen));
SRR11_9=mean(VV(Pen).*XS((Pen),3))^2;
std11_9=std(VV(Pen).*XS((Pen),3));
SRR11_10=mean(VV(Pen).*XS((Pen),4))^2;
std11_10=std(VV(Pen).*XS((Pen),4));
SRR11_11=mean(VV(Pen).*XS((Pen),5))^2;
std11_11=std(VV(Pen).*XS((Pen),5));

SRR11=SRR11_1/std11_1+SRR11_2/std11_2+SRR11_3/std11_3+SRR11_4/std11_4...
   +SRR11_5/std11_5+SRR11_6/std11_6+SRR11_7/std11_7+SRR11_9/std11_9...
   +SRR11_10/std11_10+SRR11_11/std11_11;




regXS_C=[XSEV_(:,1),XSEV_(:,2),XSEV_(:,1).^2,XSEV_(:,2).^2,log(XSEV_(:,1)),XSEV_(:,1).*XSEV_(:,2),...
    XSEV_(:,1).^2.*(-1*SameE_V),XSEV_(:,2).^2.*(-1*PartyE_V),(-1*SameE_V),(-1*PartyE_V)];

Continuec=max(0,[ones(length(LOGW_NXT_E_VC),1),LOGW_NXT_E_VC,LOGW_NXT_E_VC.^2/10,LOGW_NXT_E_VC.^3/100,...
    q_C_E_VCT,LOGW_NXT_E_VC.*q_C_E_VCT,...
    exp(q_C_E_VCT),...
    zeros(length(LOGW_NXT_E_VC),2),...
    LOGW_NXT_E_VC.*PresdumE_V,q_C_E_VCT.*PresdumE_V,...
    LOGW_NXT_E_VC.*XSEV_(:,1).*(-1*SameE_V),LOGW_NXT_E_VC.*XSEV_(:,2).*(-1*PartyE_V),...
    q_C_E_VCT.*XSEV_(:,1).*(-1*SameE_V),...
    zeros(length(LOGW_NXT_E_VC),2),XSEV_(:,1).*(-1*SameE_V),XSEV_(:,1).^2.*(-1*SameE_V),XSEV_(:,2).*(-1*PartyE_V),XSEV_(:,2).^2.*(-1*PartyE_V),...
    PresdumE_V,MidtermE_V,PresdumE_V.*MidtermE_V,...
    (-1*SameE_V).*XSEV_(:,2),(-1*SameE_V).*(-1*PartyE_V),zeros(length(LOGW_NXT_E_VC),1),(-1*SameE_V).*PresdumE_V,(-1*SameE_V).*MidtermE_V,...
    (-1*PartyE_V).*XSEV_(:,1),zeros(length(LOGW_NXT_E_VC),1),(-1*PartyE_V).*MidtermE_V,...
    PresdumE_V.*XSEV_(:,1),zeros(length(LOGW_NXT_E_VC),1),...
    MidtermE_V.*XSEV_(:,1),MidtermE_V.*XSEV_(:,2),regXS_C.*LOGW_NXT_E_VC,regXS_C]*coef2ndstep);

Derivc=max(0.001,[zeros(length(Continuec),1),ones(length(Continuec),1),2*LOGW_NXT_E_VC/10,3*LOGW_NXT_E_VC.^2/100,...
    zeros(length(Continuec),1),q_C_E_VCT,zeros(length(Continuec),1),zeros(length(Continuec),2),...
    PresdumE_V,zeros(length(Continuec),1),...
    XSEV_(:,1).*(-1*SameE_V),XSEV_(:,2).*(-1*PartyE_V),zeros(length(Continuec),22),regXS_C,zeros(length(Continuec),nregXS)]*coef2ndstep);


BX1c=abs(1/sig)*(-1)*(B_I*LOGD_E_V-B_I*LOGD_E_VC+[XSEV_(:,1),XSEV_(:,1).^2,XSEV_(:,1).*SameE_V,XSEV_(:,1).^2.*SameE_V,XSEV_(:,2).*PartyE_V]*Est2(4:8)+XQEV-q_C_E_VCT);
%Because q_C_E_VCT is determined by the same equation with BX1c replaced by
%BX1i, BX1c=BX1i.

expratio31=exp(LOGD_E_VC)./exp(LOGTot_E_VC);
FOC31=costc*beta*(1/beta)...
   .*(ones(length(q_C_E_VCT),1)+a*(1./exp(q_C_E_VCT))).*(LOGTot_E_VC).^(beta-1).*expratio31...      %% d/dI C_I(total)
   -(B_I/sig).*normpdf(BX1c).*(BB+vdelta*Continuec)-alpha*ben1*(LOGD_E_VC.^(alpha-1)); %% d/d_I(P_2)(B+delta*EV+(d/dI)H_I())

expratio32=exp(LOGW_NXT_E_VC)./exp(LOGTot_E_VC);
FOC32=costc*beta*(1/beta)...
   .*(ones(length(q_C_E_VCT),1)+a*(1./exp(q_C_E_VCT))).*(LOGTot_E_VC).^(beta-1).*expratio32...  
   -vdelta*normcdf(BX1c).*Derivc;

FOC321=FOC32.*(LOGW_NXT_E_VC>0.1);
FOC322=min(0,FOC32).*(LOGW_NXT_E_VC<0.1);

FOC32=FOC321+FOC322;

SRR21P=mean(FOC31(Pen).^2,1);
std21P=std(FOC31(Pen).^2);

SRR23=mean(FOC32(Pen).^2,1);
std23=std(FOC32(Pen).^2);


%SRR4step=10^10*SRR10P;
SRR4step=10^10*SRR10P+... %Incumbent FOC
    10^10*SRR21P;%+...%+SRR23; %Challenger FOC
    %SRR11;%+... %Orthogonality between vote share errors and observables

    

end

