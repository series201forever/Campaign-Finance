function SRR2step=Minimize2S_new(thetain,theta2,probwin,datasetV,datasetVCT,datasetVNCT,N,C,T,residstd)

%% thetain : the estimated first step estimator.
%% theta2 : the parameters that we estimate in the 2nd step.


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
XQEVCT2=datasetVCT(:,17);

%Derication of incumbent FOC uncontested
LOGW_NXT_E_VNCT=datasetVNCT(:,2);
LOGTotal_E_VNCT=datasetVNCT(:,3);
NCE_VNCT=datasetVNCT(:,4:6);
LOGTot_NCE_VNCT=datasetVNCT(:,7);
XQEVNCT2=datasetVNCT(:,8);


%LOGD_E_VC(LOGD_E_VC<9.21)=9.21;



vdelta=0.90;
% thetaS=thetain(1:2,1);
B_I=thetain(1,1);
B_C=thetain(2,1);
B_T=thetain(5,1);
thetaS2=thetain(3:4,1);
% for k=1:13
%     ACDa(:,k)=thetain(24+12*(k-1):24+12*k-1,1);
% end
% for k=1:13
%     ACDt(:,k)=thetain(180+12*(k-1):180+12*k-1,1);
% end
% for k=1:13
%     gammaACD(:,k)=thetain(336+12*(k-1):336+12*k-1,1);
% end


%h=1;
cost1=abs(theta2(1,1));  %% coefficient on the cost function of incumbent, contested
ben1=abs(theta2(2,1));   %% coefficient on the benefit function of incumbent, contested
ben2=ben1;
%ben2=abs(theta2(3,1));   %% coefficient on the benefit function of incumbent, uncontested
% cost_c=theta2(5,1); %% coefficient on the cost function of challenger
% ben_c=theta2(6,1);  %% coefficient on the benefit funciton of challenger.
% RES=theta2(7,1); %% reservation value

alpha=abs(theta2(4,1));
% alpha=cdf('norm',theta2(3,1),0,1);
% beta=1+abs(theta2(4,1));
beta=abs(theta2(3,1));
cost2=1;%abs(theta2(6,1));  %% coefficient on the cost function of the incumbent, uncontested >normalized to 1. Can be normalized
          % because we make C_I(x)=(x^alpha)*exp(f(q)) and we treat f(q)
          % non parametrically. If we set f(q)=C+f(q) then
          % C_I(x)=C*(x^alpha)*exp(f(q))
sig=abs(theta2(5,1));




%Continuation(ST,q_I,Ten,w_I,epswh, epsump, E_VCTa, E_VCTt, gamma,dF_gamma_ct, dF_total_ct, dF_nxt_nxt,Winrnd,thetawin,Ret,Betawh,Betaump,cost,ben,thetaS,Party)
%%%%%%%%          E_V         %%%%%%%%%
%% Given State, Tenure, warchest, compute incumbent continuation value %%
%Contcontest=payoff if contested
costpara=permute(repmat((ones(N,1)*((exp(LOGTot_NCE_V(:,1))./exp(NCE_V(:,1))).*...
                    ((max(0,NCE_V(:,1)).^(alpha-1))./(max(0,LOGTot_NCE_V(:,1)).^(beta-1)))).'),[1 1 10]),[3 1 2]);
                
Contcontest=zeros(3,T,N,length(NCE_V));
Contcontest(1,:,:,:)=ben1*(max(0,C(2,:,:,:))).^alpha.*C(1,:,:,:); %Benefit
Contcontest(2,:,:,:)=-cost1*(alpha/beta)*ben2*costpara.*squeeze(C(3,:,:,:)).^beta.*squeeze(C(1,:,:,:)); %Cost
Contcontest(3,:,:,:)=C(5,:,:,:).*C(1,:,:,:); %Return from winning

%Discount
for j=1:T
    Contcontest(:,j,:,:)=((vdelta)^(j-1))*Contcontest(:,j,:,:);
end
%Sum over benefit+cost+winning and over t=1~T
Contcontestpath=squeeze(sum(sum(Contcontest,1),2));

%Contuncontest=payoff if uncontested
Contuncontest=zeros(3,T,N,length(NCE_V));
Contuncontest(1,:,:,:)=ben2*(max(0,C(2,:,:,:))).^alpha.*(1-C(1,:,:,:)); %Benefit
Contuncontest(2,:,:,:)=-cost2*(alpha/beta)*ben2*costpara.*squeeze(C(3,:,:,:)).^beta.*squeeze(1-C(1,:,:,:)); %Cost
Contuncontest(3,:,:,:)=C(5,:,:,:).*(1-C(1,:,:,:)); %Return from winning

%Discount
for j=1:T
    Contuncontest(:,j,:,:)=((vdelta)^(j-1))*Contuncontest(:,j,:,:);
end

%Sum over benefit+cost+winning and over t=1~T
Contuncontestpath=squeeze(sum(sum(Contuncontest,1),2));

%Sum over contest and uncontest
Continue1=(Contcontestpath+Contuncontestpath).';

%This gives us sample size* number of simulation matrix of continuation
%payoff


% Continue1=zeros(length(NCE_V),N);
% %DContinue1=zeros(length(NCE_V),N);
% for i=1:length(NCE_V)  %% (all elements of E_V)
%     for k=1:N % simulation
%         for j=1:(min(find(C(5,:,k,i)==0))) % number of periods
%             if C(1,j,k,i)==1 %Contest
%                 Continue1(i,k)=Continue1(i,k)+((vdelta)^(j-1))*(ben1*max(0,C(2,j,k,i))^alpha-cost1*(alpha/beta)*ben2*(exp(LOGTot_NCE_V(i,1))/exp(NCE_V(i,1)))*...
%                     ((max(0,NCE_V(i,1))^(alpha-1))/(max(0,LOGTot_NCE_V(i,1))^(beta-1)))*C(3,j,k,i)^beta+C(5,j,k,i));
%             else
%                 Continue1(i,k)=Continue1(i,k)+((vdelta)^(j-1))*(ben2*max(0,C(2,j,k,i))^alpha-cost2*(alpha/beta)*ben2*(exp(LOGTot_NCE_V(i,1))/exp(NCE_V(i,1)))*...
%                     ((max(0,NCE_V(i,1))^(alpha-1))/(max(0,LOGTot_NCE_V(i,1))^(beta-1)))*C(3,j,k,i)^beta+C(5,j,k,i));
%             end
% %             if DC(1,j,k,i)==1 %Contest
% %                 DContinue1(i,k)=DContinue1(i,k)+((0.9)^(j-1))*(ben1*max(0,DC(2,j,k,i))^(1/2)-cost1*(1/4)*(ben2/cost2)*RTotDE_V(i,1)*DC(3,j,k,i)^2+DC(5,j,k,i));
% %             else
% %                 DContinue1(i,k)=DContinue1(i,k)+((0.9)^(j-1))*(ben2*max(0,DC(2,j,k,i))^(1/2)-ben2*(1/4)*RTotDE_V(i,1)*DC(3,j,k,i)^2+DC(5,j,k,i));
% %             end
%         end
%     end
% end
Continue=mean(Continue1,2);
Continuation1=Continue;
Continuation1(E_VContestFUL,:)=[];           %% Continuation1 is the Continuation value for periods in which incumebent is contested.
Continuation2=Continue;
Continuation2(E_VNContestFUL,:)=[];          %%  Continuation1 is the Continuation value for periods in which incumebent is uncontested.

%Regress Continue on State variables to find the derivative. %
Regressand=[ones(length(NCE_V),1),LOGW_NXT_E_V,LOGW_NXT_E_V.^2/10,LOGW_NXT_E_V.^3/100,LOGW_NXT_E_V.^4/1000,XQEV2,LOGW_NXT_E_V.*XQEV2,LOGW_NXT_E_V.^2/10.*XQEV2,LOGW_NXT_E_V.^3/100.*XQEV2,LOGW_NXT_E_V.^4/1000.*XQEV2,TenureE_V,XSEV_(:,1).*SameE_V,XSEV_(:,2).*PartyE_V,PresdumE_V,MidtermE_V];
%Outlier=(LOGW_NXT_E_V<quantile(LOGW_NXT_E_V,.05));
%Outlier=[];
%Continue_san_OL=Continue;
%Continue_san_OL(Outlier,:)=[];
%Regressand_san_OL=Regressand;
%Regressand_san_OL(Outlier,:)=[];
coef=Regressand\Continue;
%coef=(inv(Regressand_san_OL'*Regressand_san_OL))*Regressand_san_OL'*Continue_san_OL;

%Use estimated derivative in computing FOC.
Deriv=[zeros(length(Continuation1),1),ones(length(Continuation1),1),2*LOGW_NXT_E_VCT/10,3*LOGW_NXT_E_VCT.^2/100,4*LOGW_NXT_E_VCT.^3/1000,zeros(length(Continuation1),1),XQEVCT2,2*LOGW_NXT_E_VCT/10.*XQEVCT2,3*LOGW_NXT_E_VCT.^2/100.*XQEVCT2,4*LOGW_NXT_E_VCT.^3/1000.*XQEVCT2,zeros(length(Continuation1),5)]*coef;
DerivNCT=[zeros(length(Continuation2),1),ones(length(Continuation2),1),2*LOGW_NXT_E_VNCT/10,3*LOGW_NXT_E_VNCT.^2/100,4*LOGW_NXT_E_VNCT.^3/1000,zeros(length(Continuation2),1),XQEVNCT2,2*LOGW_NXT_E_VNCT/10.*XQEVNCT2,3*LOGW_NXT_E_VNCT.^2/100.*XQEVNCT2,4*LOGW_NXT_E_VNCT.^3/1000.*XQEVNCT2,zeros(length(Continuation2),5)]*coef;
Deriv=max(Deriv,0.00001);
DerivNCT=max(DerivNCT,0.00001);
% Regressand=[ones(length(NCE_V),1),LOGW_NXT_E_V,XQEV,LOGW_NXT_E_V.^2,LOGW_NXT_E_V.*XQEV,XSEV.*PartyE_V,log(TenureE_V+1)];
% coef=(inv(Regressand'*Regressand))*Regressand'*Continue;
% Deriv=[zeros(length(Continuation1),1),ones(length(Continuation1),1),zeros(length(Continuation1),1),2*LOGW_NXT_E_VCT,XQEVCT,zeros(length(Continuation1),2)]*coef;
% DerivNCT=[zeros(length(Continuation2),1),ones(length(Continuation2),1),zeros(length(Continuation2),1),2*LOGW_NXT_E_VNCT,XQEVNCT,zeros(length(Continuation2),2)]*coef;

OUT=((beta*cost1*(alpha/beta)*ben2*(exp(LOGTot_NCE_VCT(:,1))./exp(NCE_VCT(:,1))).*...
    ((max(0,NCE_VCT(:,1)).^(alpha-1))./(max(0,LOGTot_NCE_VCT(:,1)).^(beta-1))).*LOGTotal_E_VCT.^(beta-1)).*(1./exp(LOGTotal_E_VCT(:,1))))./((vdelta*Deriv).*(1./exp(LOGW_NXT_E_VCT)));

SRR13=mean(abs(OUT-probwin),1);
std13=std(OUT-probwin);
SRR14=(OUT>1);
SRR15=(OUT<min(probwin));

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
BX1sig=BX1;
%BX1sig(IND6CT,:)=[];


q_C_E_VCT=(-1)*(sig*BX1-B_I*LOGD_E_VCT-B_C*LOGD_E_VC-[XS_EVCT_(:,1).*SameE_VCT,XS_EVCT_(:,2).*PartyE_VCT]*thetaS2-B_T*log(TenureE_VCT+1)-XQEVCT2);


% XXX=Pen;


VV=VSEVCT-0.5-B_I*LOGD_E_VCT-B_C*LOGD_E_VC-[XS_EVCT_(:,1).*SameE_VCT,XS_EVCT_(:,2).*PartyE_VCT]*thetaS2-XQEVCT2+q_C_E_VCT-B_T*log(TenureE_VCT+1);
%VV(IND6CT,:)=[];
VV(DEL_Pen,:)=[];
LOGD_E_VCT_m=LOGD_E_VCT;
%LOGD_E_VCT_m(IND6CT,:)=[];
LOGD_E_VCT_m(DEL_Pen,:)=[];
LOGD_E_VC_m=LOGD_E_VC;
%LOGD_E_VC_m(IND6CT,:)=[];
LOGD_E_VC_m(DEL_Pen,:)=[];
XS_m=[XS_EVCT_(:,1).*SameE_VCT,XS_EVCT_(:,2).*PartyE_VCT];
%XS_m(IND6CT,:)=[];
XS_m(DEL_Pen,:)=[];
XQEVCT_m=XQEVCT2;
%XQEVCT_m(IND6CT,:)=[];
XQEVCT_m(DEL_Pen,:)=[];
q_C_E_VCT_m=q_C_E_VCT;
%q_C_E_VCT_m(IND6CT,:)=[];
q_C_E_VCT_m(DEL_Pen,:)=[];
TenureE_VCT_m=TenureE_VCT;
%TenureE_VCT_m(IND6CT,:)=[];
TenureE_VCT_m(DEL_Pen,:)=[];
SRR11_1=mean(VV.^2);     %To estimate sigma.
std11_1=std(VV);
SRR11_2=mean(VV.*LOGD_E_VCT_m)^2;
std11_2=std(VV.*LOGD_E_VCT_m);
SRR11_3=mean(VV.*LOGD_E_VC_m)^2;
std11_3=std(VV.*LOGD_E_VC_m);
SRR11_4=mean(VV.*XS_m(:,1))^2;
std11_4=std(VV.*XS_m(:,1));
SRR11_5=mean(VV.*XS_m(:,2))^2;
std11_5=std(VV.*XS_m(:,2));
SRR11_6=mean(VV.*XQEVCT_m)^2;
std11_6=std(VV.*XQEVCT_m);
SRR11_7=mean(VV.*q_C_E_VCT_m)^2;
std11_7=std(VV.*q_C_E_VCT_m);
SRR11_8=mean(VV.*log(TenureE_VCT_m+1))^2;
std11_8=std(VV.*log(TenureE_VCT_m+1));
SRR11=SRR11_1/std11_1+SRR11_2/std11_2+SRR11_3/std11_3+SRR11_4/std11_4...
    +SRR11_5/std11_5+SRR11_6/std11_6+SRR11_7/std11_7+SRR11_8/std11_8;
%VV(DEL2,:)=[];
%SRR11=(1/290)*(sum(VV))^2;
%SRR11P=(1/290)*sum(VV.^2);

SRR9=-log(sig)-(1/length(VV))*sum((VV.^2)/(2*sig^2));
SRR9=-SRR9;
%SRR9P=(1/290)*sum(VV.^2-sig^2);
%SRR9P_1=mean(((VV-mean(VV)).^2-sig^2))^2;
%std9P_1=std(((VV-mean(VV)).^2-sig^2));
%SRR9P_2=mean(sqrt(VV.^2)-sig)^2;
%std9P_2=std(sqrt(VV.^2)-sig);
%% foc of incumbent, contested periods %%
FOC11=(beta*cost1*(alpha/beta)*ben2*(exp(LOGTot_NCE_VCT(:,1))./exp(NCE_VCT(:,1))).*...
    ((max(0,NCE_VCT(:,1)).^(alpha-1))./(max(0,LOGTot_NCE_VCT(:,1)).^(beta-1))).*LOGTotal_E_VCT.^(beta-1)).*(1./exp(LOGTotal_E_VCT(:,1)))...      %% d/dI C_I(total)
    -(B_I/(sig*exp(LOGD_E_VCT(:,1))))*normpdf(BX1).*(1+vdelta*Continuation1)-alpha*ben1*(1./exp(LOGD_E_VCT(:,1))).*(max(0,(LOGD_E_VCT))).^(alpha-1); %% d/d_I(P_2)(B+delta*EV+(d/dI)H_I())

%FOC11(IND6CT,:)=[];
FOC11(DEL_Pen,:)=[];
%FOC11(DEL2,:)=[];
SRR10P=mean(FOC11.^2,1);
std10P=std(FOC11);
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
FOC21=vdelta*DerivNCT.*(1./exp(LOGW_NXT_E_VNCT))-... 
(beta*cost2*(alpha/beta)*ben2*(exp(LOGTot_NCE_VNCT(:,1))./exp(NCE_VNCT(:,1))).*...
    ((max(0,NCE_VNCT(:,1)).^(alpha-1))./(max(0,LOGTot_NCE_VNCT(:,1)).^(beta-1))).*LOGTotal_E_VNCT.^(beta-1)).*(1./exp(LOGTotal_E_VNCT(:,1)));   %% delta*(d/dw_I)E_V-C_I'
%FOC21(IND6NCT,:)=[];

% FOC21delta=0.9*(DContinuation2b2-Continuation2b2)-(ben2+incr)*RTotDE_VNCT.*((LOGTotal_E_VNCT+Delt).^2-LOGTotal_E_VNCT.^2);   %% delta*(d/dw_I)E_V-C_I'
% FOC21delta(IND6NCT,:)=[];
% SRR12=((1/incr)*((1/290)*sum(FOC21delta.^2,1)-(1/290)*sum(FOC21.^2,1)))^2;
SRR12P=mean(FOC21.^2,1);
std12P=std(FOC21);

% SRR9
% SRR10P
% % SRR11P
% SRR12P

SRR2step=SRR9+10^7*SRR10P/std10P+SRR11+10^4*SRR12P/std12P+10*sum(SRR15)+10*sum(SRR14)+0.0000*SRR13/std13;%+10000*(sigma>residstd);%+;%;%

%if bestiter==iterate-1
 %   save Best2step.txt theta2 SRR2step bestiter -ASCII
  %  save q_C_E_VCT.mat q_C_E_VCT
   % save BX1.mat BX1
   % save coef.mat coef
   % save DEL_Pen.mat DEL_Pen
%     subplot(2,1,1)
%     hist(XQEVCT2,30)
%     subplot(2,1,2)
%     hist(q_C_E_VCT,30)
%     saveas(figure(1), 'qual' , 'jpg')
end

