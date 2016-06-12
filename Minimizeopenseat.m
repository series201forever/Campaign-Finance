function SRR4step=Minimizeopenseat(thetain,theta3,theta4,q_C_E_VCT,datasetV)

%% thetain : the estimated first step estimator.
%% theta2 : the parameters that we estimate in the 2nd step.


%Calculation of derivative of payoff
XQEV=datasetV(:,1);
LOGW_NXT_E_V=datasetV(:,2);
LOGTotal_E_V=datasetV(:,3);
LOGD_E_V=datasetV(:,4);
Continuei=datasetV(:,5);
Deriv=datasetV(:,6);


vdelta=0.90;
% thetaS=thetain(1:2,1);
% B_I=thetain(1,1);
% B_C=thetain(2,1);
% B_T=thetain(5,1);
% thetaS2=thetain(3:4,1);


alpha=1/2;
beta=2;
ben1=abs(thetain(2));
ben2=ben1;
sig=abs(thetain(3));


costc=abs(theta3(1));
a=theta3(2);
b=theta3(3);
c=theta3(4);
B_I=abs(theta4(1));


OUT=((beta*costc*(alpha/beta)*ben2*(ones(length(XQEV),1)+a*(XQEV+abs(min(q_C_E_VCT))+0.001)+b*(XQEV+abs(min(q_C_E_VCT))+0.001).^2+c*(XQEV+abs(min(q_C_E_VCT))+0.001).^3).*LOGTotal_E_V.^(beta-1)).*(1./exp(LOGTotal_E_V(:,1))))./((vdelta*Deriv).*(1./exp(LOGW_NXT_E_V)));

 

%DE=exp(LOGTotal_E_VCT(:,1)).*((vdelta*Deriv).*(1./exp(LOGW_NXT_E_VCT)));

SRR14=(OUT>1);
SRR15=(OUT<0);

Pen=max(OUT-1,0)-min(OUT,0);
Pen2=Pen;
%Pen2(IND6CT,:)=[];
DEL_Pen=[find(Pen2<quantile(Pen2,.025));find(Pen2>quantile(Pen2,.975))];    %Cannot invert some observations:
%DEL2=find(Pen2>0);
%Pen2(DEL_Pen,:)=[];

%ako=1/2+(1/pi)*atan((OUT-0.5)*h);
ako=min(max(OUT,0.000001),0.999999);
BX1=norminv(ako);
% BX1=norminv(max(0.01,min(0.99,(beta*(alpha/beta)*ben2*...
%     ((max(0,NCE_VCT(:,1)).^(alpha-1))./(max(0,LOGTot_NCE_VCT(:,1)).^(beta-1))).*LOGTotal_E_VCT.^(beta-1))./(vdelta*Deriv))));  %%BX1=(1/sig)*(-0.5+B_I*d_I-B_C*d_C+...)=(K) in the paper.
%BX1sig=BX1;
%BX1sig(IND6CT,:)=[];


%% foc
%FOC11=(beta*cost1*(alpha/beta)*ben2*(exp(LOGTot_NCE_VCT(:,1))./exp(NCE_VCT(:,1))).*...
%    ((max(0,NCE_VCT(:,1)).^(alpha-1))./(max(0,LOGTot_NCE_VCT(:,1)).^(beta-1))).*LOGTotal_E_VCT.^(beta-1)).*(1./exp(LOGTotal_E_VCT(:,1)))...      %% d/dI C_I(total)
%    -(B_I/(sig*exp(LOGD_E_VCT(:,1))))*normpdf(BX1).*(1+vdelta*Continuation1)-alpha*ben1*(1./exp(LOGD_E_VCT(:,1))).*(max(0,(LOGD_E_VCT))).^(alpha-1); %% d/d_I(P_2)(B+delta*EV+(d/dI)H_I())
FOC11=(beta*costc*(alpha/beta)*ben2*(ones(length(XQEV),1)+a*(XQEV+abs(min(q_C_E_VCT))+0.001)+b*(XQEV+abs(min(q_C_E_VCT))+0.001).^2+c*(XQEV+abs(min(q_C_E_VCT))+0.001).^3).*LOGTotal_E_V.^(beta-1)).*(1./exp(LOGTotal_E_V(:,1)))...      %% d/dI C_I(total)
    -(B_I./(sig*exp(LOGD_E_V(:,1)))).*normpdf(BX1).*(1+vdelta*Continuei)-alpha*ben1*(1./exp(LOGD_E_V(:,1))).*(max(0,(LOGD_E_V))).^(alpha-1); %% d/d_I(P_2)(B+delta*EV+(d/dI)H_I())

%%%%%%%%%%%%%%%%%%%%%%
%Moment conditions concerning "truncated mean = realized epsilon" needed to identify B_I.
%%%%%%%%%%%%%%%%%%%%%%

%FOC11(IND6CT,:)=[];
FOC11(DEL_Pen,:)=[];
%FOC11(DEL2,:)=[];
SRR10P=mean(FOC11.^2,1);
std10P=std(FOC11);

SRR4step=10^12*SRR10P/std10P+10^3*sum(SRR15)+10^3*sum(SRR14);%+0.0000*SRR13/std13;%+10000*(sigma>residstd);%+;%;%

end

