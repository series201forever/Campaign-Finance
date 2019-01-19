function [SRR4step,out,OUTI,OUTC,BX1i,BX1c,SRR11i,SRR11c,OUT2i,Derivi0,Derivc0]=Minimizequality(Est2,mintheta2ndstage,coef2ndstep,estopen,q_C_E_VCT,datasetV,XQEV,const,id)
%%

%This is the function to be minimized to obtain quality estimates, given
%parameter values estimated in the previous steps
%if id=1, Incumbent vs Challenger. id=2 uncontested. id=3, openseat.

%Calculation of derivative of payoff
LOGW_NXT_E_V=datasetV(:,1);
LOGTotal_E_V=datasetV(:,2);
LOGD_E_V=datasetV(:,3);
PartyE_V=datasetV(:,4);
SameE_V=datasetV(:,5);
PresdumE_V=datasetV(:,6);
MidtermE_V=datasetV(:,7);
XSEV_=datasetV(:,8:9);
TenureE_V=datasetV(:,10);
LOGW_NXT_E_VC=datasetV(:,11);
LOGTotal_E_VC=datasetV(:,12);
LOGD_E_VC=datasetV(:,13);

vdelta=0.90;
% thetaS=thetain(1:2,1);
% B_I=thetain(1,1);
% B_C=thetain(2,1);
% B_T=thetain(5,1);
% thetaS2=thetain(3:4,1);

if length(mintheta2ndstage)==9
 %alpha=abs(mintheta2ndstage(5));
% betai=abs(mintheta2ndstage(6));
alpha=abs(mintheta2ndstage(7));
alphac=alpha;%abs(mintheta2ndstage(8));
betai=abs(mintheta2ndstage(6));
betac=abs(mintheta2ndstage(9));
ben1=abs(mintheta2ndstage(2));
ben2=ben1;
sig=abs(mintheta2ndstage(3));

costi=abs(mintheta2ndstage(1));
a=abs(mintheta2ndstage(4));

costc=abs(mintheta2ndstage(8));
ac=a;
end
if length(mintheta2ndstage)==8
  %alpha=abs(mintheta2ndstage(5));
% betai=abs(mintheta2ndstage(6));
alpha=abs(mintheta2ndstage(7));
alphac=abs(mintheta2ndstage(8));
betai=abs(mintheta2ndstage(6));
betac=betai;%abs(mintheta2ndstage(9));
ben1=abs(mintheta2ndstage(2));
ben2=ben1;
sig=abs(mintheta2ndstage(3));

costi=abs(mintheta2ndstage(1));
a=abs(mintheta2ndstage(4));

costc=costi;%abs(mintheta2ndstage(8));
ac=a;   
end

%CC=Est2(1)+const;
B_I=Est2(2);
B_C=Est2(3);
B_state=Est2(4:8);
B_T=Est2(9);

B_O=estopen(1);

qi=XQEV(1);
if id==1||id==3
qc=XQEV(2);
end

if id==1 %If contested, incumbent and challenger face different parameters
%Compute continuation payoff and its derivative for both people

regXS=[XSEV_(:,1),XSEV_(:,2),XSEV_(:,1).^2,XSEV_(:,2).^2,log(XSEV_(:,1)),XSEV_(:,1).*XSEV_(:,2),...
    XSEV_(:,1).^2.*SameE_V,XSEV_(:,2).^2.*PartyE_V,SameE_V,PartyE_V];
nregXS=size(regXS,2);
aux=max(0,[ones(length(LOGW_NXT_E_V),1),LOGW_NXT_E_V,LOGW_NXT_E_V.^2/10,LOGW_NXT_E_V.^3/100,...
    qi,LOGW_NXT_E_V.*qi,...
    exp(qi),...
    LOGW_NXT_E_V.*TenureE_V,qi.*TenureE_V,...
    LOGW_NXT_E_V.*PresdumE_V,qi.*PresdumE_V,...
    LOGW_NXT_E_V.*XSEV_(:,1).*SameE_V,LOGW_NXT_E_V.*XSEV_(:,2).*PartyE_V,...
    qi.*XSEV_(:,1).*SameE_V,...
    TenureE_V,TenureE_V.^2,XSEV_(:,1).*SameE_V,XSEV_(:,1).^2.*SameE_V,XSEV_(:,2).*PartyE_V,XSEV_(:,2).^2.*PartyE_V,...
    PresdumE_V,MidtermE_V,PresdumE_V.*MidtermE_V,...
    SameE_V.*XSEV_(:,2),SameE_V.*PartyE_V,SameE_V.*TenureE_V,SameE_V.*PresdumE_V,SameE_V.*MidtermE_V,...
    PartyE_V.*XSEV_(:,1),PartyE_V.*TenureE_V,PartyE_V.*MidtermE_V,...
    PresdumE_V.*XSEV_(:,1),PresdumE_V.*TenureE_V,...
    MidtermE_V.*XSEV_(:,1),MidtermE_V.*XSEV_(:,2),regXS.*LOGW_NXT_E_V,regXS]*coef2ndstep);

if coef2ndstep(7)>0
threshquali=log(-1/coef2ndstep(7)*(coef2ndstep(5)+LOGW_NXT_E_V*coef2ndstep(6)+TenureE_V*coef2ndstep(9)+PresdumE_V*coef2ndstep(11)+XSEV_(:,1).*SameE_V*coef2ndstep(14)));
Continuei=(threshquali<qi).*aux;
else
    threshquali=-99999;
    Continuei=aux;
end


regXS_VC=[XSEV_(:,1),XSEV_(:,2),XSEV_(:,1).^2,XSEV_(:,2).^2,log(XSEV_(:,1)),XSEV_(:,1).*XSEV_(:,2),...
    XSEV_(:,1).^2.*(-1*SameE_V),XSEV_(:,2).^2.*(-1*PartyE_V),(-1*SameE_V),(-1*PartyE_V)];

auxc=max(0,[ones(length(LOGW_NXT_E_VC),1),LOGW_NXT_E_VC,LOGW_NXT_E_VC.^2/10,LOGW_NXT_E_VC.^3/100,...
    qc,LOGW_NXT_E_VC.*qc,...
    exp(qc),...
    zeros(length(LOGW_NXT_E_VC),2),...
    LOGW_NXT_E_VC.*PresdumE_V,qc.*PresdumE_V,...
    LOGW_NXT_E_VC.*XSEV_(:,1).*(-1*SameE_V),LOGW_NXT_E_VC.*XSEV_(:,2).*(-1*PartyE_V),...   
    qc.*XSEV_(:,1).*(-1*SameE_V),...
    zeros(length(LOGW_NXT_E_VC),2),XSEV_(:,1).*(-1*SameE_V),XSEV_(:,1).^2.*(-1*SameE_V),XSEV_(:,2).*(-1*PartyE_V),XSEV_(:,2).^2.*(-1*PartyE_V),...
    PresdumE_V,MidtermE_V,PresdumE_V.*MidtermE_V,...
    (-1*SameE_V).*XSEV_(:,2),(-1*SameE_V).*(-1*PartyE_V),(-1*SameE_V).*TenureE_V,(-1*SameE_V).*PresdumE_V,(-1*SameE_V).*MidtermE_V,...
    (-1*PartyE_V).*XSEV_(:,1),(-1*PartyE_V).*TenureE_V,(-1*PartyE_V).*MidtermE_V,...
    PresdumE_V.*XSEV_(:,1),PresdumE_V.*TenureE_V,...
    MidtermE_V.*XSEV_(:,1),MidtermE_V.*XSEV_(:,2),regXS_VC.*LOGW_NXT_E_VC,regXS_VC]*coef2ndstep);

if coef2ndstep(7)>0
threshqualc=log(-1/coef2ndstep(7)*(coef2ndstep(5)+LOGW_NXT_E_VC*coef2ndstep(6)+PresdumE_V*coef2ndstep(11)+XSEV_(:,1).*(-1*SameE_V)*coef2ndstep(14)));
Continuec=(threshqualc<qc).*auxc;
else
    threshqualc=-99999;
    Continuec=auxc;
end

Derivi=max(0.0001,[zeros(length(LOGW_NXT_E_V),1),ones(length(LOGW_NXT_E_V),1),2*LOGW_NXT_E_V/10,3*LOGW_NXT_E_V.^2/100,...
    zeros(length(LOGW_NXT_E_V),1),qi,zeros(length(LOGW_NXT_E_V),1),TenureE_V,zeros(length(LOGW_NXT_E_V),1),...
    PresdumE_V,zeros(length(LOGW_NXT_E_V),1),...
    XSEV_(:,1).*SameE_V,XSEV_(:,2).*PartyE_V,zeros(length(LOGW_NXT_E_V),22),regXS,zeros(length(LOGW_NXT_E_V),nregXS)]*coef2ndstep)+(threshquali>qi)*0.0001;

Derivi0=[zeros(length(LOGW_NXT_E_V),1),ones(length(LOGW_NXT_E_V),1),2*LOGW_NXT_E_V/10,3*LOGW_NXT_E_V.^2/100,...
    zeros(length(LOGW_NXT_E_V),1),qi,zeros(length(LOGW_NXT_E_V),1),TenureE_V,zeros(length(LOGW_NXT_E_V),1),...
    PresdumE_V,zeros(length(LOGW_NXT_E_V),1),...
    XSEV_(:,1).*SameE_V,XSEV_(:,2).*PartyE_V,zeros(length(LOGW_NXT_E_V),22),regXS,zeros(length(LOGW_NXT_E_V),nregXS)]*coef2ndstep;

Derivc=max(0.0001,[zeros(length(LOGW_NXT_E_VC),1),ones(length(LOGW_NXT_E_VC),1),2*LOGW_NXT_E_VC/10,3*LOGW_NXT_E_VC.^2/100,...
    zeros(length(LOGW_NXT_E_VC),1),qc,...
    zeros(length(LOGW_NXT_E_VC),1),...
    zeros(length(LOGW_NXT_E_VC),1),zeros(length(LOGW_NXT_E_VC),1),...
    PresdumE_V,zeros(length(LOGW_NXT_E_VC),1),...
    XSEV_(:,1).*(-1*SameE_V),XSEV_(:,2).*(-1*PartyE_V),zeros(length(LOGW_NXT_E_VC),22),regXS_VC,zeros(length(LOGW_NXT_E_VC),nregXS)]*coef2ndstep)+(threshqualc>qc)*0.0001;

Derivc0=[zeros(length(LOGW_NXT_E_VC),1),ones(length(LOGW_NXT_E_VC),1),2*LOGW_NXT_E_VC/10,3*LOGW_NXT_E_VC.^2/100,...
    zeros(length(LOGW_NXT_E_VC),1),qc,...
    zeros(length(LOGW_NXT_E_VC),1),...
    zeros(length(LOGW_NXT_E_VC),1),zeros(length(LOGW_NXT_E_VC),1),...
    PresdumE_V,zeros(length(LOGW_NXT_E_VC),1),...
    XSEV_(:,1).*(-1*SameE_V),XSEV_(:,2).*(-1*PartyE_V),zeros(length(LOGW_NXT_E_VC),22),regXS_VC,zeros(length(LOGW_NXT_E_VC),nregXS)]*coef2ndstep;


%Calculate K
OUTI=((betai*costi*(1/betai)*(max(10^(-8),ones(length(qi),1)+a*(1./exp(qi)))).*LOGTotal_E_V.^(betai-1)).*(1./exp(LOGTotal_E_V(:,1))))./((vdelta*Derivi).*(1./exp(LOGW_NXT_E_V)));
OUTC=((betac*costc*(1/betac)*(max(10^(-8),ones(length(qc),1)+ac*(1./exp(qc)))).*LOGTotal_E_VC.^(betac-1)).*(1./exp(LOGTotal_E_VC(:,1))))./((vdelta*Derivc).*(1./exp(LOGW_NXT_E_VC)));

 
SRR14=(OUTI>1);
SRR15=(OUTI<0);
SRR16=(OUTC>1);
SRR17=(OUTC<0);

Pen=max(OUTI-1,0)-min(OUTI,0);
%  Pen2=Pen;
%  DEL_Peni=[find(Pen2<quantile(Pen2,.025));find(Pen2>quantile(Pen2,.975))];    %Cannot invert some observations:
ako=min(max(OUTI,0.000001),0.999999);
BX1i=norminv(ako);

 Penc=max(OUTC-1,0)-min(OUTC,0);
%  Pen2c=Penc;
%  DEL_Penc=[find(Pen2c<quantile(Pen2c,.025));find(Pen2c>quantile(Pen2c,.975))];    %Cannot invert some observations:
akoc=min(max(OUTC,0.000001),0.999999);
BX1c=norminv(akoc);


 BX2i=1./sig*(const+B_I*LOGD_E_V+B_C*LOGD_E_VC+qi-qc+[XSEV_(:,1),XSEV_(:,1).^2,XSEV_(:,1).*SameE_V,XSEV_(:,1).^2.*SameE_V,XSEV_(:,2).*PartyE_V]*B_state+B_T*log(1+TenureE_V));
 BX2c=1./sig*(-const-B_C*LOGD_E_VC-B_I*LOGD_E_V+qc-qi-[XSEV_(:,1),XSEV_(:,1).^2,XSEV_(:,1).*SameE_V,XSEV_(:,1).^2.*SameE_V,XSEV_(:,2).*PartyE_V]*B_state-B_T*log(1+TenureE_V));

% 
  OUT2i=normcdf(BX2i);
  OUT2c=1-OUT2i;
%Two FOCs used as moments
expratio11i=exp(LOGD_E_V(:,1))./exp(LOGTotal_E_V(:,1));
FOC11i=(betai*costi*(1/betai)*(ones(length(qi),1)+a*(1./exp(qi))).*LOGTotal_E_V.^(betai-1)).*expratio11i...      %% d/dI C_I(total)
    -(B_I./sig).*normpdf(BX1i).*(1+vdelta*Continuei)-alpha*ben1.*(max(0,(LOGD_E_V))).^(alpha-1); %% d/d_I(P_2)(B+delta*EV+(d/dI)H_I())
expratio11c=exp(LOGD_E_VC(:,1))./exp(LOGTotal_E_VC(:,1));
FOC11c=(betac*costc*(1/betac)*(ones(length(qc),1)+ac*(1./exp(qc))).*LOGTotal_E_VC.^(betac-1)).*expratio11c...      %% d/dI C_I(total)
    +(B_C./sig).*normpdf(BX1c).*(1+vdelta*Continuec)-alphac*ben1.*(max(0,(LOGD_E_VC))).^(alphac-1); %% d/d_I(P_2)(B+delta*EV+(d/dI)H_I())


%  FOC11i(DEL_Peni2,:)=[];
%  FOC11c(DEL_Penc2,:)=[];

SRR10Pi=mean(FOC11i.^2,1);
%std10Pi=std(FOC11i.^2);
SRR10Pc=mean(FOC11c.^2,1);
%std10Pc=std(FOC11c.^2);


%Use K to calculate two ex-ante vote shares
 %SRR11i=mean((BX1i-BX2i).^2,1);
 %SRR11c=mean((BX1c-BX2c).^2,1);
 SRR11i=mean((OUTI-OUT2i).^2,1);
 SRR11c=mean((OUTC-OUT2c).^2,1);
%SRR12=mean((OUTI+OUTC-1).^2,1);

SRR4step=10^6*SRR10Pi+10^6*SRR10Pc+SRR11i+SRR11c;%+SRR12;%+sum(SRR14)+sum(SRR15)+sum(SRR16)+sum(SRR17);%+0.0000*SRR13/std13;%+10000*(sigma>residstd);%+;%;%

end


if id==3 %If openseat, both people face the same parameters
%Compute continuation payoff and its derivative for both people

regXS=[XSEV_(:,1),XSEV_(:,2),XSEV_(:,1).^2,XSEV_(:,2).^2,log(XSEV_(:,1)),XSEV_(:,1).*XSEV_(:,2),...
    XSEV_(:,1).^2.*SameE_V,XSEV_(:,2).^2.*PartyE_V,SameE_V,PartyE_V];
nregXS=size(regXS,2);

aux=max(0,[ones(length(LOGW_NXT_E_V),1),LOGW_NXT_E_V,LOGW_NXT_E_V.^2/10,LOGW_NXT_E_V.^3/100,...
    qi,LOGW_NXT_E_V.*qi,...
    exp(qi),...
    zeros(length(LOGW_NXT_E_V),2),...
    LOGW_NXT_E_V.*PresdumE_V,qi.*PresdumE_V,...
    LOGW_NXT_E_V.*XSEV_(:,1).*SameE_V,LOGW_NXT_E_V.*XSEV_(:,2).*PartyE_V,...
    qi.*XSEV_(:,1).*SameE_V,...
    zeros(length(LOGW_NXT_E_V),2),XSEV_(:,1).*SameE_V,XSEV_(:,1).^2.*SameE_V,XSEV_(:,2).*PartyE_V,XSEV_(:,2).^2.*PartyE_V,...
    PresdumE_V,MidtermE_V,PresdumE_V.*MidtermE_V,...
    SameE_V.*XSEV_(:,2),SameE_V.*PartyE_V,zeros(length(LOGW_NXT_E_V),1),SameE_V.*PresdumE_V,SameE_V.*MidtermE_V,...
    PartyE_V.*XSEV_(:,1),zeros(length(LOGW_NXT_E_V),1),PartyE_V.*MidtermE_V,...
    PresdumE_V.*XSEV_(:,1),zeros(length(LOGW_NXT_E_V),1),...
    MidtermE_V.*XSEV_(:,1),MidtermE_V.*XSEV_(:,2),regXS.*LOGW_NXT_E_V,regXS]*coef2ndstep);

if coef2ndstep(7)>0
threshquali=log(-1/coef2ndstep(7)*(coef2ndstep(5)+LOGW_NXT_E_V*coef2ndstep(6)+PresdumE_V*coef2ndstep(11)+XSEV_(:,1).*SameE_V*coef2ndstep(14)));
Continuei=(threshquali<qi).*aux;
else
threshquali=-99999;
Continuei=aux;
end

regXS_VC=[XSEV_(:,1),XSEV_(:,2),XSEV_(:,1).^2,XSEV_(:,2).^2,log(XSEV_(:,1)),XSEV_(:,1).*XSEV_(:,2),...
    XSEV_(:,1).^2.*(-1*SameE_V),XSEV_(:,2).^2.*(-1*PartyE_V),(-1*SameE_V),(-1*PartyE_V)];

auxc=max(0,[ones(length(LOGW_NXT_E_VC),1),LOGW_NXT_E_VC,LOGW_NXT_E_VC.^2/10,LOGW_NXT_E_VC.^3/100,...
    qc,LOGW_NXT_E_VC.*qc,...
    exp(qc),...
    zeros(length(LOGW_NXT_E_VC),2),...
    LOGW_NXT_E_VC.*PresdumE_V,qc.*PresdumE_V,...
    LOGW_NXT_E_VC.*XSEV_(:,1).*(-1*SameE_V),LOGW_NXT_E_VC.*XSEV_(:,2).*(-1*PartyE_V),...   
    qc.*XSEV_(:,1).*(-1*SameE_V),...
    zeros(length(LOGW_NXT_E_VC),2),XSEV_(:,1).*(-1*SameE_V),XSEV_(:,1).^2.*(-1*SameE_V),XSEV_(:,2).*(-1*PartyE_V),XSEV_(:,2).^2.*(-1*PartyE_V),...
    PresdumE_V,MidtermE_V,PresdumE_V.*MidtermE_V,...
    (-1*SameE_V).*XSEV_(:,2),(-1*SameE_V).*(-1*PartyE_V),zeros(length(LOGW_NXT_E_VC),1),(-1*SameE_V).*PresdumE_V,(-1*SameE_V).*MidtermE_V,...
    (-1*PartyE_V).*XSEV_(:,1),zeros(length(LOGW_NXT_E_VC),1),(-1*PartyE_V).*MidtermE_V,...
    PresdumE_V.*XSEV_(:,1),zeros(length(LOGW_NXT_E_VC),1),...
    MidtermE_V.*XSEV_(:,1),MidtermE_V.*XSEV_(:,2),regXS_VC.*LOGW_NXT_E_VC,regXS_VC]*coef2ndstep);

if coef2ndstep(7)>0
threshqualc=log(-1/coef2ndstep(7)*(coef2ndstep(5)+LOGW_NXT_E_VC*coef2ndstep(6)+PresdumE_V*coef2ndstep(11)+XSEV_(:,1).*(-1*SameE_V)*coef2ndstep(14)));
Continuec=(threshqualc<qc).*auxc;
else
threshqualc=-99999;
Continuec=auxc;
end

Derivi=max(0.0001,[zeros(length(LOGW_NXT_E_V),1),ones(length(LOGW_NXT_E_V),1),2*LOGW_NXT_E_V/10,3*LOGW_NXT_E_V.^2/100,...
    zeros(length(LOGW_NXT_E_V),1),qi,zeros(length(LOGW_NXT_E_V),3),...
    PresdumE_V,zeros(length(LOGW_NXT_E_V),1),...
    XSEV_(:,1).*SameE_V,XSEV_(:,2).*PartyE_V,zeros(length(LOGW_NXT_E_V),22),regXS,zeros(length(LOGW_NXT_E_V),nregXS)]*coef2ndstep)+(threshquali>qi)*0.0001;
%(threshquali<qi).*
Derivc=max(0.0001,[zeros(length(LOGW_NXT_E_VC),1),ones(length(LOGW_NXT_E_VC),1),2*LOGW_NXT_E_VC/10,3*LOGW_NXT_E_VC.^2/100,...
    zeros(length(LOGW_NXT_E_VC),1),qc,zeros(length(LOGW_NXT_E_VC),3),...
    PresdumE_V,zeros(length(LOGW_NXT_E_VC),1),...
    XSEV_(:,1).*(-1*SameE_V),XSEV_(:,2).*(-1*PartyE_V),zeros(length(LOGW_NXT_E_VC),22),regXS_VC,zeros(length(LOGW_NXT_E_VC),nregXS)]*coef2ndstep)+(threshqualc>qc)*0.0001;
%(threshqualc<qc).*




%Calculate K
OUTI=(betac*costc*(1/betac)*(max(10^(-8),ones(length(qi),1)+ac*(1./exp(qi))).*LOGTotal_E_V.^(betac-1)).*(1./exp(LOGTotal_E_V(:,1))))./((vdelta*Derivi).*(1./exp(LOGW_NXT_E_V)));
OUTC=(betac*costc*(1/betac)*(max(10^(-8),ones(length(qc),1)+ac*(1./exp(qc))).*LOGTotal_E_VC.^(betac-1)).*(1./exp(LOGTotal_E_VC(:,1))))./((vdelta*Derivc).*(1./exp(LOGW_NXT_E_VC)));
%OUTI=((beta*costc*(alpha/beta)*ben2*(ones(length(qi),1)+ac*(qi+abs(min(q_C_E_VCT))+0.001)+bc*(qi+abs(min(q_C_E_VCT))+0.001).^2+cc*(qi+abs(min(q_C_E_VCT))+0.001).^3).*LOGTotal_E_V.^(beta-1)).*(1./exp(LOGTotal_E_V(:,1))))./((vdelta*Derivi).*(1./exp(LOGW_NXT_E_V)));
%OUTC=((beta*costc*(alpha/beta)*ben2*(ones(length(qc),1)+ac*(qc+abs(min(q_C_E_VCT))+0.001)+bc*(qc+abs(min(q_C_E_VCT))+0.001).^2+cc*(qc+abs(min(q_C_E_VCT))+0.001).^3).*LOGTotal_E_VC.^(beta-1)).*(1./exp(LOGTotal_E_VC(:,1))))./((vdelta*Derivc).*(1./exp(LOGW_NXT_E_VC)));

 
SRR14=(OUTI>1);
SRR15=(OUTI<0);
SRR16=(OUTC>1);
SRR17=(OUTC<0);

 Pen=max(OUTI-1,0)-min(OUTI,0);
%  Pen2=Pen;
%  DEL_Peni=[find(Pen2<quantile(Pen2,.025));find(Pen2>quantile(Pen2,.975))];    %Cannot invert some observations:
ako=min(max(OUTI,0.000001),0.999999);
BX1i=norminv(ako);

 Penc=max(OUTC-1,0)-min(OUTC,0);
%  Pen2c=Penc;
%  DEL_Penc=[find(Pen2c<quantile(Pen2c,.025));find(Pen2c>quantile(Pen2c,.975))];    %Cannot invert some observations:
akoc=min(max(OUTC,0.000001),0.999999);
BX1c=norminv(akoc);

BX2i=1./sig*(B_O*LOGD_E_V-B_O*LOGD_E_VC+qi-qc+[XSEV_(:,1),XSEV_(:,1).^2,XSEV_(:,1).*SameE_V,XSEV_(:,1).^2.*SameE_V,XSEV_(:,2).*PartyE_V]*B_state);
BX2c=1./sig*(B_O*LOGD_E_VC-B_O*LOGD_E_V+qc-qi-[XSEV_(:,1),XSEV_(:,1).^2,XSEV_(:,1).*SameE_V,XSEV_(:,1).^2.*SameE_V,XSEV_(:,2).*PartyE_V]*B_state);
% OUTI2=normcdf(BX2i);
% OUTC2=normcdf(BX2c);
%Two FOCs used as moments 
FOC11i=(betac*costc*(1/betac)*(ones(length(qi),1)+ac*(1./exp(qi))).*LOGTotal_E_V.^(betac-1)).*(1./exp(LOGTotal_E_V(:,1)))...      %% d/dI C_I(total)
    -(B_O./(sig*exp(LOGD_E_V(:,1)))).*normpdf(BX1i).*(1+vdelta*Continuei)-alphac*ben1*(1./exp(LOGD_E_V(:,1))).*(max(0,(LOGD_E_V))).^(alphac-1); %% d/d_I(P_2)(B+delta*EV+(d/dI)H_I())
FOC11c=(betac*costc*(1/betac)*(ones(length(qc),1)+ac*(1./exp(qc))).*LOGTotal_E_VC.^(betac-1)).*(1./exp(LOGTotal_E_VC(:,1)))...      %% d/dI C_I(total)
    -(B_O./(sig*exp(LOGD_E_VC(:,1)))).*normpdf(BX1c).*(1+vdelta*Continuec)-alphac*ben1*(1./exp(LOGD_E_VC(:,1))).*(max(0,(LOGD_E_VC))).^(alphac-1); %% d/d_I(P_2)(B+delta*EV+(d/dI)H_I())


% FOC11i(DEL_Peni,:)=[];
% FOC11c(DEL_Penc,:)=[];

SRR10Pi=mean(FOC11i.^2,1);
%std10Pi=std(FOC11i);
SRR10Pc=mean(FOC11c.^2,1);
%std10Pc=std(FOC11c);
% 
% 
% %Use K to calculate two ex-ante vote shares
SRR11i=mean((BX1i-BX2i).^2,1);
SRR11c=mean((BX1c-BX2c).^2,1);

%SRR12=mean((OUTI+OUTC-1).^2,1);

SRR4step=10^6*SRR10Pi+10^6*SRR10Pc+SRR11i+SRR11c;%+sum(SRR14)+sum(SRR15)+sum(SRR16)+sum(SRR17);%+SRR12;

end


if id==2 %If uncontested, only incumbent needed
    
    
regXS=[XSEV_(:,1),XSEV_(:,2),XSEV_(:,1).^2,XSEV_(:,2).^2,log(XSEV_(:,1)),XSEV_(:,1).*XSEV_(:,2),...
    XSEV_(:,1).^2.*SameE_V,XSEV_(:,2).^2.*PartyE_V,SameE_V,PartyE_V];
nregXS=size(regXS,2);

Derivi=max(0.0001,[zeros(length(LOGW_NXT_E_V),1),ones(length(LOGW_NXT_E_V),1),2*LOGW_NXT_E_V/10,3*LOGW_NXT_E_V.^2/100,...
    zeros(length(LOGW_NXT_E_V),1),qi,zeros(length(LOGW_NXT_E_V),1),TenureE_V,zeros(length(LOGW_NXT_E_V),1),...
    PresdumE_V,zeros(length(LOGW_NXT_E_V),1),...
    XSEV_(:,1).*SameE_V,XSEV_(:,2).*PartyE_V,zeros(length(LOGW_NXT_E_V),22),regXS,zeros(length(LOGW_NXT_E_V),nregXS)]*coef2ndstep);


%Two FOCs used as moments 
expratio11i=exp(LOGD_E_V(:,1))./exp(LOGTotal_E_V(:,1));
FOC11i=(betai*costi*(1/betai)*(ones(length(qi),1)+a*(1./exp(qi))).*LOGTotal_E_V.^(betai-1)).*expratio11i...      %% d/dI C_I(total)
    -alpha*ben1.*(max(0,(LOGD_E_V))).^(alpha-1); 

expratio12i=exp(LOGW_NXT_E_V)./exp(LOGTotal_E_V(:,1));
FOC12i=(betai*costi*(1/betai)*(ones(length(qi),1)+a*(1./exp(qi))).*LOGTotal_E_V.^(betai-1)).*expratio12i...
    -(vdelta*Derivi);

SRR10Pi=mean(FOC11i.^2,1);
SRR11Pi=mean(FOC12i.^2,1);
%std10Pi=std(FOC11i);


%Use K to calculate two ex-ante vote shares
% SRR11i=mean((BX1i-1./sig*(B_I*LOGD_E_V+B_C*LOGD_E_VC+qi-qc+XSEV_*B_state+B_T*log(1+TenureE_V))).^2,1);
% SRR11c=mean((BX1c-1./sig*(-B_C*LOGD_E_VC-B_I*LOGD_E_V+qc-qi-XSEV_*B_state)).^2,1);
% 
% SRR12=mean((OUTI+OUTC-1).^2,1);

SRR4step=10^6*SRR10Pi+10^6*SRR11Pi;%+;%;%

end

if id==1||id==3
out=max((SRR15+SRR14),(SRR16+SRR17));
end

%out=min((OUTI>=0.99)+(OUTI<=0)+(OUTC>=0.99)+(OUTC<=0),1);
end



