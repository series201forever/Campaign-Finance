function [SRR4step,out,OUTI,OUTC,BX1i,BX1c,SRR11i,SRR11c,test]=Minimizequality(Est2,mintheta2ndstage,est3rdstage,coef3rdstep,estopen,q_C_E_VCT,datasetV,XQEV,id)

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


alpha=1/2;
beta=2;
ben1=abs(mintheta2ndstage(2));
ben2=ben1;
sig=abs(mintheta2ndstage(3));

costi=abs(mintheta2ndstage(1));
a=mintheta2ndstage(4);
b=mintheta2ndstage(5);
c=mintheta2ndstage(6);

costc=abs(est3rdstage(1));
ac=est3rdstage(2);
bc=est3rdstage(3);
cc=est3rdstage(4);


B_I=Est2(1);
B_C=Est2(2);
B_O=abs(estopen(1));
B_state=Est2(3:4);
B_T=Est2(5);

qi=XQEV(1);
if id==1||id==3
qc=XQEV(2);
end

if id==1 %If contested, incumbent and challenger face different parameters
%Compute continuation payoff and its derivative for both people

Continuei=max(0,[ones(length(LOGW_NXT_E_V),1),LOGW_NXT_E_V,LOGW_NXT_E_V.^2/10,LOGW_NXT_E_V.^3/100,LOGW_NXT_E_V.^4/1000,...
    qi,LOGW_NXT_E_V.*qi,LOGW_NXT_E_V.^2/10.*qi,LOGW_NXT_E_V.^3/100.*qi,LOGW_NXT_E_V.^4/1000.*qi,...
    TenureE_V,XSEV_(:,1).*SameE_V,XSEV_(:,2).*PartyE_V,...
    PresdumE_V,MidtermE_V]*coef3rdstep);
Continuec=max(0,[ones(length(LOGW_NXT_E_VC),1),LOGW_NXT_E_VC,LOGW_NXT_E_VC.^2/10,LOGW_NXT_E_VC.^3/100,LOGW_NXT_E_VC.^4/1000,...
    qc,LOGW_NXT_E_VC.*qc,LOGW_NXT_E_VC.^2/10.*qc,LOGW_NXT_E_VC.^3/100.*qc,LOGW_NXT_E_VC.^4/1000.*qc,...
    zeros(length(LOGW_NXT_E_VC),1),XSEV_(:,1).*(-1).*SameE_V,XSEV_(:,2).*(-1).*PartyE_V,...
    PresdumE_V,MidtermE_V]*coef3rdstep);


Derivi=max(0.0001,[zeros(length(LOGW_NXT_E_V),1),ones(length(LOGW_NXT_E_V),1),2*LOGW_NXT_E_V/10,3*LOGW_NXT_E_V.^2/100,4*LOGW_NXT_E_V.^3/1000,...
    zeros(length(LOGW_NXT_E_V),1),qi,2*LOGW_NXT_E_V/10.*qi,3*LOGW_NXT_E_V.^2/100.*qi,4*LOGW_NXT_E_V.^3/1000.*qi,zeros(length(LOGW_NXT_E_V),5)]*coef3rdstep);
Derivc=max(0.0001,[zeros(length(LOGW_NXT_E_VC),1),ones(length(LOGW_NXT_E_VC),1),2*LOGW_NXT_E_VC/10,3*LOGW_NXT_E_VC.^2/100,4*LOGW_NXT_E_VC.^3/1000,...
    zeros(length(LOGW_NXT_E_VC),1),qc,2*LOGW_NXT_E_VC/10.*qc,3*LOGW_NXT_E_VC.^2/100.*qc,4*LOGW_NXT_E_VC.^3/1000.*qc,zeros(length(LOGW_NXT_E_VC),5)]*coef3rdstep);
% Derivi=[zeros(length(LOGW_NXT_E_V),1),ones(length(LOGW_NXT_E_V),1),2*LOGW_NXT_E_V/10,3*LOGW_NXT_E_V.^2/100,4*LOGW_NXT_E_V.^3/1000,...
%      zeros(length(LOGW_NXT_E_V),1),qi,2*LOGW_NXT_E_V/10.*qi,3*LOGW_NXT_E_V.^2/100.*qi,4*LOGW_NXT_E_V.^3/1000.*qi,zeros(length(LOGW_NXT_E_V),5)]*coef3rdstep;
% Derivc=[zeros(length(LOGW_NXT_E_VC),1),ones(length(LOGW_NXT_E_VC),1),2*LOGW_NXT_E_VC/10,3*LOGW_NXT_E_VC.^2/100,4*LOGW_NXT_E_VC.^3/1000,...
%      zeros(length(LOGW_NXT_E_VC),1),qc,2*LOGW_NXT_E_VC/10.*qc,3*LOGW_NXT_E_VC.^2/100.*qc,4*LOGW_NXT_E_VC.^3/1000.*qc,zeros(length(LOGW_NXT_E_VC),5)]*coef3rdstep;




%Calculate K
OUTI=((beta*costi*(alpha/beta)*ben2*(max(10^(-8),ones(length(qi),1)+a*qi+b*qi.^2+c*qi.^3)).*LOGTotal_E_V.^(beta-1)).*(1./exp(LOGTotal_E_V(:,1))))./((vdelta*Derivi).*(1./exp(LOGW_NXT_E_V)));
OUTC=((beta*costc*(alpha/beta)*ben2*(max(10^(-8),ones(length(qc),1)+ac*(qc+abs(min(q_C_E_VCT))+0.001)+bc*(qc+abs(min(q_C_E_VCT))+0.001).^2+cc*(qc+abs(min(q_C_E_VCT))+0.001).^3)).*LOGTotal_E_VC.^(beta-1)).*(1./exp(LOGTotal_E_VC(:,1))))./((vdelta*Derivc).*(1./exp(LOGW_NXT_E_VC)));

 
SRR14=(OUTI>1);
SRR15=(OUTI<0);
SRR16=(OUTC>1);
SRR17=(OUTC<0);

 Pen=max(OUTI-1,0)-min(OUTI,0);
 Pen2=Pen;
 DEL_Peni=[find(Pen2<quantile(Pen2,.025));find(Pen2>quantile(Pen2,.975))];    %Cannot invert some observations:
ako=min(max(OUTI,0.000001),0.999999);
BX1i=norminv(ako);

 Penc=max(OUTC-1,0)-min(OUTC,0);
 Pen2c=Penc;
 DEL_Penc=[find(Pen2c<quantile(Pen2c,.025));find(Pen2c>quantile(Pen2c,.975))];    %Cannot invert some observations:
akoc=min(max(OUTC,0.000001),0.999999);
BX1c=norminv(akoc);

%Two FOCs used as moments 
FOC11i=(beta*costi*(alpha/beta)*ben2*(ones(length(qi),1)+a*qi+b*qi.^2+c*qi.^3).*LOGTotal_E_V.^(beta-1)).*(1./exp(LOGTotal_E_V(:,1)))...      %% d/dI C_I(total)
    -(B_I./(sig*exp(LOGD_E_V(:,1)))).*normpdf(BX1i).*(1+vdelta*Continuei)-alpha*ben1*(1./exp(LOGD_E_V(:,1))).*(max(0,(LOGD_E_V))).^(alpha-1); %% d/d_I(P_2)(B+delta*EV+(d/dI)H_I())
FOC11c=(beta*costc*(alpha/beta)*ben2*(ones(length(qc),1)+ac*(qc+abs(min(q_C_E_VCT))+0.001)+bc*(qc+abs(min(q_C_E_VCT))+0.001).^2+cc*(qc+abs(min(q_C_E_VCT))+0.001).^3).*LOGTotal_E_VC.^(beta-1)).*(1./exp(LOGTotal_E_VC(:,1)))...      %% d/dI C_I(total)
    +(B_C./(sig*exp(LOGD_E_VC(:,1)))).*normpdf(BX1c).*(1+vdelta*Continuec)-alpha*ben1*(1./exp(LOGD_E_VC(:,1))).*(max(0,(LOGD_E_VC))).^(alpha-1); %% d/d_I(P_2)(B+delta*EV+(d/dI)H_I())


FOC11i(DEL_Peni,:)=[];
FOC11c(DEL_Penc,:)=[];

SRR10Pi=mean(FOC11i.^2,1);
%std10Pi=std(FOC11i);
SRR10Pc=mean(FOC11c.^2,1);
%std10Pc=std(FOC11c);


%Use K to calculate two ex-ante vote shares
SRR11i=mean((BX1i-1./sig*(B_I*LOGD_E_V+B_C*LOGD_E_VC+qi-qc+XSEV_*B_state+B_T*log(1+TenureE_V))).^2,1);
SRR11c=mean((BX1c-1./sig*(-B_C*LOGD_E_VC-B_I*LOGD_E_V+qc-qi-XSEV_*B_state)).^2,1);

SRR12=mean((OUTI+OUTC-1).^2,1);

SRR4step=SRR10Pi+SRR10Pc+SRR11i+SRR11c+SRR12+sum(SRR14)+sum(SRR15)+sum(SRR16)+sum(SRR17);%+0.0000*SRR13/std13;%+10000*(sigma>residstd);%+;%;%

test=1./sig*(B_I*LOGD_E_V+B_C*LOGD_E_VC+qi-qc+XSEV_*B_state+B_T*log(1+TenureE_V));
end



if id==3 %If openseat, both people face the same parameters
%Compute continuation payoff and its derivative for both people

Continuei=max(0,[ones(length(LOGW_NXT_E_V),1),LOGW_NXT_E_V,LOGW_NXT_E_V.^2/10,LOGW_NXT_E_V.^3/100,LOGW_NXT_E_V.^4/1000,...
    qi,LOGW_NXT_E_V.*qi,LOGW_NXT_E_V.^2/10.*qi,LOGW_NXT_E_V.^3/100.*qi,LOGW_NXT_E_V.^4/1000.*qi,...
    zeros(length(LOGW_NXT_E_VC),1),XSEV_(:,1).*SameE_V,XSEV_(:,2).*PartyE_V,...
    PresdumE_V,MidtermE_V]*coef3rdstep);
Continuec=max(0,[ones(length(LOGW_NXT_E_VC),1),LOGW_NXT_E_VC,LOGW_NXT_E_VC.^2/10,LOGW_NXT_E_VC.^3/100,LOGW_NXT_E_VC.^4/1000,...
    qc,LOGW_NXT_E_VC.*qc,LOGW_NXT_E_VC.^2/10.*qc,LOGW_NXT_E_VC.^3/100.*qc,LOGW_NXT_E_VC.^4/1000.*qc,...
    zeros(length(LOGW_NXT_E_VC),1),XSEV_(:,1).*(-1).*SameE_V,XSEV_(:,2).*(-1).*PartyE_V,...
    PresdumE_V,MidtermE_V]*coef3rdstep);


Derivi=max(0.0001,[zeros(length(LOGW_NXT_E_V),1),ones(length(LOGW_NXT_E_V),1),2*LOGW_NXT_E_V/10,3*LOGW_NXT_E_V.^2/100,4*LOGW_NXT_E_V.^3/1000,...
    zeros(length(LOGW_NXT_E_V),1),qi,2*LOGW_NXT_E_V/10.*qi,3*LOGW_NXT_E_V.^2/100.*qi,4*LOGW_NXT_E_V.^3/1000.*qi,zeros(length(LOGW_NXT_E_V),5)]*coef3rdstep);
Derivc=max(0.0001,[zeros(length(LOGW_NXT_E_VC),1),ones(length(LOGW_NXT_E_VC),1),2*LOGW_NXT_E_VC/10,3*LOGW_NXT_E_VC.^2/100,4*LOGW_NXT_E_VC.^3/1000,...
    zeros(length(LOGW_NXT_E_VC),1),qc,2*LOGW_NXT_E_VC/10.*qc,3*LOGW_NXT_E_VC.^2/100.*qc,4*LOGW_NXT_E_VC.^3/1000.*qc,zeros(length(LOGW_NXT_E_VC),5)]*coef3rdstep);
% Derivi=[zeros(length(LOGW_NXT_E_V),1),ones(length(LOGW_NXT_E_V),1),2*LOGW_NXT_E_V/10,3*LOGW_NXT_E_V.^2/100,4*LOGW_NXT_E_V.^3/1000,...
%      zeros(length(LOGW_NXT_E_V),1),qi,2*LOGW_NXT_E_V/10.*qi,3*LOGW_NXT_E_V.^2/100.*qi,4*LOGW_NXT_E_V.^3/1000.*qi,zeros(length(LOGW_NXT_E_V),5)]*coef3rdstep;
% Derivc=[zeros(length(LOGW_NXT_E_VC),1),ones(length(LOGW_NXT_E_VC),1),2*LOGW_NXT_E_VC/10,3*LOGW_NXT_E_VC.^2/100,4*LOGW_NXT_E_VC.^3/1000,...
%      zeros(length(LOGW_NXT_E_VC),1),qc,2*LOGW_NXT_E_VC/10.*qc,3*LOGW_NXT_E_VC.^2/100.*qc,4*LOGW_NXT_E_VC.^3/1000.*qc,zeros(length(LOGW_NXT_E_VC),5)]*coef3rdstep;



%Calculate K
OUTI=((beta*costc*(alpha/beta)*ben2*(max(10^(-8),ones(length(qi),1)+ac*(qi+abs(min(q_C_E_VCT))+0.001)+bc*(qi+abs(min(q_C_E_VCT))+0.001).^2+cc*(qi+abs(min(q_C_E_VCT))+0.001).^3)).*LOGTotal_E_V.^(beta-1)).*(1./exp(LOGTotal_E_V(:,1))))./((vdelta*Derivi).*(1./exp(LOGW_NXT_E_V)));
OUTC=((beta*costc*(alpha/beta)*ben2*(max(10^(-8),ones(length(qc),1)+ac*(qc+abs(min(q_C_E_VCT))+0.001)+bc*(qc+abs(min(q_C_E_VCT))+0.001).^2+cc*(qc+abs(min(q_C_E_VCT))+0.001).^3)).*LOGTotal_E_VC.^(beta-1)).*(1./exp(LOGTotal_E_VC(:,1))))./((vdelta*Derivc).*(1./exp(LOGW_NXT_E_VC)));
%OUTI=((beta*costc*(alpha/beta)*ben2*(ones(length(qi),1)+ac*(qi+abs(min(q_C_E_VCT))+0.001)+bc*(qi+abs(min(q_C_E_VCT))+0.001).^2+cc*(qi+abs(min(q_C_E_VCT))+0.001).^3).*LOGTotal_E_V.^(beta-1)).*(1./exp(LOGTotal_E_V(:,1))))./((vdelta*Derivi).*(1./exp(LOGW_NXT_E_V)));
%OUTC=((beta*costc*(alpha/beta)*ben2*(ones(length(qc),1)+ac*(qc+abs(min(q_C_E_VCT))+0.001)+bc*(qc+abs(min(q_C_E_VCT))+0.001).^2+cc*(qc+abs(min(q_C_E_VCT))+0.001).^3).*LOGTotal_E_VC.^(beta-1)).*(1./exp(LOGTotal_E_VC(:,1))))./((vdelta*Derivc).*(1./exp(LOGW_NXT_E_VC)));

 
SRR14=(OUTI>1);
SRR15=(OUTI<0);
SRR16=(OUTC>1);
SRR17=(OUTC<0);

 Pen=max(OUTI-1,0)-min(OUTI,0);
 Pen2=Pen;
 DEL_Peni=[find(Pen2<quantile(Pen2,.025));find(Pen2>quantile(Pen2,.975))];    %Cannot invert some observations:
ako=min(max(OUTI,0.000001),0.999999);
BX1i=norminv(ako);

 Penc=max(OUTC-1,0)-min(OUTC,0);
 Pen2c=Penc;
 DEL_Penc=[find(Pen2c<quantile(Pen2c,.025));find(Pen2c>quantile(Pen2c,.975))];    %Cannot invert some observations:
akoc=min(max(OUTC,0.000001),0.999999);
BX1c=norminv(akoc);

%Two FOCs used as moments 
FOC11i=(beta*costc*(alpha/beta)*ben2*(ones(length(qi),1)+ac*(qi+abs(min(q_C_E_VCT))+0.001)+bc*(qi+abs(min(q_C_E_VCT))+0.001).^2+cc*(qi+abs(min(q_C_E_VCT))+0.001).^3).*LOGTotal_E_V.^(beta-1)).*(1./exp(LOGTotal_E_V(:,1)))...      %% d/dI C_I(total)
    -(B_O./(sig*exp(LOGD_E_V(:,1)))).*normpdf(BX1i).*(1+vdelta*Continuei)-alpha*ben1*(1./exp(LOGD_E_V(:,1))).*(max(0,(LOGD_E_V))).^(alpha-1); %% d/d_I(P_2)(B+delta*EV+(d/dI)H_I())
FOC11c=(beta*costc*(alpha/beta)*ben2*(ones(length(qc),1)+ac*(qc+abs(min(q_C_E_VCT))+0.001)+bc*(qc+abs(min(q_C_E_VCT))+0.001).^2+cc*(qc+abs(min(q_C_E_VCT))+0.001).^3).*LOGTotal_E_VC.^(beta-1)).*(1./exp(LOGTotal_E_VC(:,1)))...      %% d/dI C_I(total)
    -(B_O./(sig*exp(LOGD_E_VC(:,1)))).*normpdf(BX1c).*(1+vdelta*Continuec)-alpha*ben1*(1./exp(LOGD_E_VC(:,1))).*(max(0,(LOGD_E_VC))).^(alpha-1); %% d/d_I(P_2)(B+delta*EV+(d/dI)H_I())


FOC11i(DEL_Peni,:)=[];
FOC11c(DEL_Penc,:)=[];

SRR10Pi=mean(FOC11i.^2,1);
%std10Pi=std(FOC11i);
SRR10Pc=mean(FOC11c.^2,1);
%std10Pc=std(FOC11c);


%Use K to calculate two ex-ante vote shares
SRR11i=mean((BX1i-1./sig*(B_O*LOGD_E_V-B_O*LOGD_E_VC+qi-qc+XSEV_*B_state)).^2,1);
SRR11c=mean((BX1c-1./sig*(B_O*LOGD_E_VC-B_O*LOGD_E_V+qc-qi-XSEV_*B_state)).^2,1);

%SRR12=mean((OUTI+OUTC-1).^2,1);

SRR4step=SRR10Pi+SRR10Pc+SRR11i+SRR11c+sum(SRR14)+sum(SRR15)+sum(SRR16)+sum(SRR17);%+SRR12;

test=1./sig*(B_O*LOGD_E_V-B_O*LOGD_E_VC+qi-qc+XSEV_*B_state);
end


if id==2 %If uncontested, only incumbent needed
%Compute continuation payoff and its derivative
% 
% Continuei=max(0,[ones(length(LOGW_NXT_E_V),1),LOGW_NXT_E_V,LOGW_NXT_E_V.^2/10,LOGW_NXT_E_V.^3/100,LOGW_NXT_E_V.^4/1000,...
%     qi,LOGW_NXT_E_V.*qi,LOGW_NXT_E_V.^2/10.*qi,LOGW_NXT_E_V.^3/100.*qi,LOGW_NXT_E_V.^4/1000.*qi,...
%     TenureE_V,XSEV_(:,1).*SameE_V,XSEV_(:,2).*PartyE_V,...
%     PresdumE_V,MidtermE_V]*coef3rdstep);
% 
% 
% Derivi=max(0.0001,[zeros(length(LOGW_NXT_E_V),1),ones(length(LOGW_NXT_E_V),1),2*LOGW_NXT_E_V/10,3*LOGW_NXT_E_V.^2/100,4*LOGW_NXT_E_V.^3/1000,...
%     zeros(length(LOGW_NXT_E_V),1),qi,2*LOGW_NXT_E_V/10.*qi,3*LOGW_NXT_E_V.^2/100.*qi,4*LOGW_NXT_E_V.^3/1000.*qi,zeros(length(LOGW_NXT_E_V),5)]*coef3rdstep);
% 
% 
% %Calculate K
% OUTI=((beta*costi*(alpha/beta)*ben2*(max(10^(-8),ones(length(qi),1)+a*qi+b*qi.^2+c*qi.^3)).*LOGTotal_E_V.^(beta-1)).*(1./exp(LOGTotal_E_V(:,1))))./((vdelta*Derivi).*(1./exp(LOGW_NXT_E_V)));
% 
%  
% SRR14=(OUTI>1);
% SRR15=(OUTI<0);


%  Pen=max(OUTI-1,0)-min(OUTI,0);
%  Pen2=Pen;
%  DEL_Peni=[find(Pen2<quantile(Pen2,.025));find(Pen2>quantile(Pen2,.975))];    %Cannot invert some observations:
% ako=min(max(OUTI,0.000001),0.999999);
% BX1i=norminv(ako);


%Two FOCs used as moments 
FOC11i=(beta*costi*(alpha/beta)*ben2*(ones(length(qi),1)+a*qi+b*qi.^2+c*qi.^3).*LOGTotal_E_V.^(beta-1)).*(1./exp(LOGTotal_E_V(:,1)))...      %% d/dI C_I(total)
    -alpha*ben1*(1./exp(LOGD_E_V(:,1))).*(max(0,(LOGD_E_V))).^(alpha-1); %% d/d_I(P_2)(B+delta*EV+(d/dI)H_I())


% FOC11i(DEL_Peni,:)=[];

SRR10Pi=mean(10^15*FOC11i.^2,1);
%std10Pi=std(FOC11i);


%Use K to calculate two ex-ante vote shares
% SRR11i=mean((BX1i-1./sig*(B_I*LOGD_E_V+B_C*LOGD_E_VC+qi-qc+XSEV_*B_state+B_T*log(1+TenureE_V))).^2,1);
% SRR11c=mean((BX1c-1./sig*(-B_C*LOGD_E_VC-B_I*LOGD_E_V+qc-qi-XSEV_*B_state)).^2,1);
% 
% SRR12=mean((OUTI+OUTC-1).^2,1);

SRR4step=SRR10Pi;%+;%;%

end

if id==1||id==3
out=max((SRR14+SRR15),(SRR16+SRR17));
end

%out=min((OUTI>=0.99)+(OUTI<=0)+(OUTC>=0.99)+(OUTC<=0),1);
end



