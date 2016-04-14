function SRR3step=Minimize3PS(thetain1,thetain2,theta3)

%% thetain : the estimated first step estimator.
%% theta2 : the parameters that we estimate in the 2nd step.

global iterate


global N
global Sofarbest
global bestiter

global Cc
global XQEV2

global q_c3step
global LOGD_E_VC3step
global rtotd3step
global LOGTot_E_VC3step
global BX13step

global LOGW_NXT3step
global XS3step_

global Party3step

global coef

B_I=thetain(1,1);
B_C=thetain(2,1);
B_T=thetain(5,1);
thetaS2=thetain(3:4,1);




cost1=abs(thetain2(1,1));  %% coefficient on the cost function of incumbent, contested
ben1=abs(thetain2(2,1));   %% coefficient on the benefit function of incumbent, contested
ben2=ben1;   %% coefficient on the benefit function of incumbent, uncontested
cost2=1;  %% coefficient on the cost function of the incumbent, uncontested >normalized to 1. Can be normalized
sig=abs(thetain2(3,1));
alpha=1/2;
%alpha=cdf('norm',thetain2(3,1),0,1);
%beta=1+abs(thetain2(4,1));
beta=2;


costc=abs(theta3(1,1)); %% coefficient on the cost function of challenger
benc=ben1;  %% coefficient on the benefit funciton of challenger
% pracrazy=abs(theta3(4,1));  %% probability of entry by crazy type.
% betacr1=abs(theta3(5,1));  %% crazy type: beta1
% betacr2=abs(theta3(6,1));  %% crazy tpe: beta2
%betara1=abs(theta3(3,1));  %% rational type: beta1
%betara2=abs(theta3(4,1));  %% rational type : beta2





% %%%%%%%%          E_V         %%%%%%%%%
% %% Given State, Tenure, warchest, compute challenger's continuation value %%
% %% use ben1 and cost1 to evaluate E_V as challenger becomes incumbemt  %%
% Continue13step=zeros(length(q_c3step),N);
% % % DContinue13step=zeros(length(q_c3step),N);
% for i=1:length(q_c3step)  %% (all elements of E_V(IND6CT,:)=[])
%     for k=1:N
%         for j=1:10
%             if Cc(1,j,k,i)==1 %Contest
%                 Continue13step(i,k)=Continue13step(i,k)+((0.9)^(j-1))*(ben1*max(0,Cc(2,j,k,i))^alpha-(alpha/beta)*cost1*ben2*...
%                     rtotd3step(i,1)*Cc(3,j,k,i)^beta+Cc(5,j,k,i)); %Note we use ben1, cost1 etc and NOT benc, costc.
%             else
%                 Continue13step(i,k)=Continue13step(i,k)+((0.9)^(j-1))*(ben2*max(0,Cc(2,j,k,i))^alpha-(alpha/beta)*ben2*...
%                     rtotd3step(i,1)*Cc(3,j,k,i)^beta+Cc(5,j,k,i));
%             end
%         end
%     end
% end
% Continue3step=mean(Continue13step,2);
% Regressand=[ones(length(q_c3step),1),LOGW_NXT3step,q_c3step,LOGW_NXT3step.*q_c3step,XS3step_.*repmat(Party3step,1,2)];
% Outlier=(LOGW_NXT3step<quantile(LOGW_NXT3step,.05));
% Continue3step_san_OL=Continue3step;
% Continue3step_san_OL(Outlier,:)=[];
% %egressand_san_OL=Regressand;
% Regressand_san_OL(Outlier,:)=[];
% %coef=(inv(Regressand_san_OL'*Regressand_san_OL))*Regressand_san_OL'*Continue3step_san_OL;
% %Regress Continue on State variables to find the derivative. %
% % Regressand=[ones(length(q_c3step),1),LOGW_NXT3step,q_c3step,LOGW_NXT3step.^2,LOGW_NXT3step.*q_c3step,XS3step.*Party3step];
% % coef=(inv(Regressand'*Regressand))*Regressand'*Continue3step;
% Continue3step=[ones(length(q_c3step),1),LOGW_NXT3step,q_c3step,ones(length(q_c3step),1),XS3step_(:,1).*Party3step,XS3step_(:,2).*Party3step]*coef;
% % Deriv=[zeros(length(q_c3step),1),ones(length(q_c3step),1),zeros(length(q_c3step),1),zeros(length(q_c3step),3)]*coef;
% % Deriv=max(Deriv,0.00001);
% 
% NOFOC=min([LOGTot_E_VC3step,LOGD_E_VC3step]')';
% IND8=find(NOFOC==0);


Continue3step=[ones(length(q_c3step),1),LOGW_NXT3step,q_c3step,ones(length(q_c3step),1),XS3step_(:,1).*Party3step,XS3step_(:,2).*Party3step]*coef;


[ones(length(NCE_V),1),LOGW_NXT_E_V,LOGW_NXT_E_V.^2/10,LOGW_NXT_E_V.^3/100,LOGW_NXT_E_V.^4/1000,...
    XQEV2,LOGW_NXT_E_V.*XQEV2,LOGW_NXT_E_V.^2/10.*XQEV2,LOGW_NXT_E_V.^3/100.*XQEV2,LOGW_NXT_E_V.^4/1000.*XQEV2,...
    TenureE_V,TenureE_V.^2,XSEV_(:,1).*SameE_V,XSEV_(:,1).^2.*SameE_V,XSEV_(:,2).*PartyE_V,XSEV_(:,2).^2.*PartyE_V,...
    TenureE_V.*XSEV_(:,1).*SameE_V,TenureE_V.*XSEV_(:,2).*PartyE_V,PresdumE_V,MidtermE_V]

% %% foc of challenger %%
FOC31=costc*beta*(alpha/beta)*ben2*(1./exp(LOGTot_E_VC3step)).*...
    rtotd3step.*(LOGTot_E_VC3step).^(beta-1)...      %% d/dI C_I(total)
    +(B_C./(sig*exp(LOGD_E_VC3step))).*normpdf(BX13step).*(1+0.9*Continue3step)-alpha*benc*(LOGD_E_VC3step.^(alpha-1)).*(1./exp(LOGD_E_VC3step)); %% d/d_I(P_2)(B+delta*EV+(d/dI)H_I())
% benefitspwin=mean(-1*(B_C/sig)*normpdf(BX13step).*(1+0.9*Continue3step))
% benefitspsp=mean(benc*((LOGD_E_VC3step+Delt).^(1/2)-(LOGD_E_VC3step).^(1/2)))
FOC31(IND8,:)=[];
HQ_Entrant=(q_c3step>=min(XQEV2));
HQ_Entrant(IND8,:)=[];
FOC31=FOC31(HQ_Entrant,:);
SRR13P=sum(FOC31,1)^2;
%XXX=(B_C/sig)*normpdf(BX13step).*(1+0.9*Continue3step)-benc*((LOGD_E_VC3step+Delt).^(1/2)-(LOGD_E_VC3step).^(1/2));

% %% foc of challenger %%
% FOC32=0.9*(1-normcdf(BX13step)).*Deriv-costc*beta*(alpha/beta)*ben2*...
%     rtotd3step.*(LOGTot_E_VC3step).^(beta-1); %% delta*(d/dw_I)E_V-C_I'
% FOC32(IND6CTto7,:)=[];  %only if warchest of challenger>0
% HQ_Entrant=(q_c3step>=min(XQEV2));
% HQ_Entrant(IND6CTto7,:)=[];
% FOC32=FOC32(HQ_Entrant,:);
% SRR14P=sum(FOC32.^2,1);
% 
% benefitsv=0.9*normcdf(BX13step).*(DContinue3step-Continue3step);
% benefitsv(IND6CTto7,:)=[];
% bensv=mean(benefitsv)
% 
% costsv=costc*(1/4)*ben2*rtotd3step.*((LOGTot_E_VC3step+Delt).^2-LOGTot_E_VC3step.^2);
% costsv(IND6CTto7,:)=[];
% costsvshow=mean(costsv)
% XXX=-(B_C/sig)*normcdf(BX13step).*(1+0.9*Continue3step)+benc*((LOGD_E_VC3step).^(1/2))...
%  -costc*(1/4)*(ben2/cost2)*rtotd3step.*LOGTot_E_VC3step.^2;
% 

%COMMENT OUT 7/25/2012
% q_c3step_=(q_c3step-qc_LB)/(qc_UB-qc_LB);
% span_ind=repmat((1:max(N_E)),length(q_c3step_),1);
% dift=max(cdf('beta',q_c3step_,betara1,betara2)-F_cutoff,10^(-5));
% SRR15PS=binopdf(span_ind,repmat(N_E,1,max(N_E)),repmat(F_cutoff,1,max(N_E))).*span_ind.*repmat(dift,1,max(N_E)).^(span_ind-1).*...
%     repmat(pdf('beta',q_c3step_,betara1,betara2),1,max(N_E))./repmat((1-F_cutoff),1,max(N_E)).^(span_ind);
% SRR15PS=log(sum(SRR15PS,2));
% SRR15PS=(1/290)*(-1)*sum(SRR15PS);


% (qc_UB-qc_LB)*icdf('beta',qEntry, betara1, betara2)+qc_LB;



% ER=(1/290)*sum(repmat((qe-q_c3step),1,4).*[q_i3step,LOGW_E_VCT3step,TenureE_VCT3step,(-1)*XS3step.*Party3step]);
% SRR15PS=ER*ER';

% 
% PrEra=(PrE3step-pracrazy)/(1-pracrazy);
% % SRR15P=(1/290)*sum(log((1-PrEra)*pracrazy.*pdf('beta',max(10^(-7),q_c3step-min(q_c3step)),betacr1,betacr2)...   % no ratinoal type entry.
% %     +PrEra.*(cdf('beta',max(10^(-7),q_c3step-min(q_c3step)),betara1,betara2)>(1-PrEra)).*pdf('beta',max(10^(-7),q_c3step-min(q_c3step)),betara1,betara2)...
% %     ./PrEra));
% % %cdf('beta',max(10^(-7),q_c3step-min(q_c3step)),betara1,betara2)
% % SRR15P=(-1)*SRR15P; %minimize (-1)*loglikelihood.
% %XXX=max(10^(-7),q_c3step-min(q_c3step))
% % SRR10
% % SRR12
% Maxq_c3step=max(q_c3step)+0.01;
% Minq_c3step=min(q_c3step)-0.01;
% Span=Maxq_c3step-Minq_c3step;
% MPrEra=mean(PrEra);
% threshold=betainv(1-MPrEra,betara1,betara2)
% TH=betainv(1-PrEra,betara1,betara2);
% SRR15PS=0;
% for i=2:10
%     cutoff=Minq_c3step+((i-1)/10)*Span;
%     COq_c3step=(q_c3step<cutoff);
% %     Drop=find(COq_c3step==1);
% %     MODq_c3step=q_c3step;
% %     MODq_c3step(Drop,:)=[];
% %     Fmomcr=Minq_c3step+Span*(1/(1-cdf('beta',((i-1)/10),betacr1,betacr2)))...
% %         *(beta(betacr1+1,betacr2)/beta(betacr1,betacr2))*(1-betainc(((i-1)/10),(betacr1+1),betacr2))
% %     Fmomra=Minq_c3step+Span*(1/(1-cdf('beta',max(((i-1)/10),threshold),betara1,betara2)))...
% %         *(beta(betara1+1,betara2)/beta(betara1,betara2))*(1-betainc(max(((i-1)/10),threshold),(betara1+1),betara2))
% %     predMean=(1/(pracrazy*(1-MPrEra)*(1-cdf('beta',((i-1)/10),betacr1,betacr2))+MPrEra*(1-cdf('beta',max(((i-1)/10),threshold),betara1,betara2))))...
% %     *(pracrazy*(1-cdf('beta',((i-1)/10),betacr1,betacr2))*(1-MPrEra)*Fmomcr+MPrEra*(1-cdf('beta',max(((i-1)/10),threshold),betara1,betara2))*Fmomra)
% %     SRR15PS=SRR15PS+sum((predMean-mean(MODq_c3step)).^2);
% %     mean(MODq_c3step)
%     
%     SampleFrac=sum(COq_c3step)/length(q_c3step)
% 
%     PredFrac=mean((1./(PrEra+pracrazy*(1-PrEra)))...
%         .*(pracrazy*(1-PrEra)*cdf('beta',((i-1)/10),betacr1,betacr2)...
%         +(1./(1-cdf('beta',TH,betara1,betara2))).*PrEra.*max(0,cdf('beta',((i-1)/10),betara1,betara2)-cdf('beta',TH,betara1,betara2))))
%     SRR15PS=SRR15PS+(SampleFrac-PredFrac).^2;
%     
% end

    
% 
% Fmom=(beta(betara1+1,betara2)/beta(betara1,betara2))*(1-betainc(min(0.9999999,max(10^(-7),threshold)),(betara1+1),betara2))./PrE3step;
% Smom=(beta(betara1+2,betara2)/beta(betara1,betara2))*(1-betainc(min(0.9999999,max(10^(-7),threshold)),(betara1+2),betara2))./PrE3step;
% 
% SRR15PS=(1/290)*sum((max(10^(-7),q_c3step-min(q_c3step))...
%     -(PrE3step.*(1./PrE3step).*Fmom)).^2); %% betara1+1 
%                                                                                     %%because E[x] w.r.t. beta is same as integrate with beta with first parameter++1
%                                                                                     
% SRR16PS=(1/290)*sum((((max(10^(-7),q_c3step-min(q_c3step))).^2)...
%     -(PrE3step.*(1./PrE3step).*Smom)).^2);
    
% U=sum(((Smom-Fmom.^2)<0));
% if U>0
%     XXX=(Smom-Fmom)
% end

% U=(1-PrEra).*(1./PrE3step)*pracrazy*(betacr1/(betacr1+betacr2))
% V=betainc(min(0.9999999,max(10^(-7),threshold)),(betara1+1),betara2)
% XXX=threshold


% %SRR13P
% SRR14P
% % SRR15P
% SRR15PS
% % SRR16PS
% Penalty=sum((PrEra<0).*(PrEra).^2)+sum((PrEra>1).*(PrEra-1).^2)

%SRR3step=SRR13P+SRR14P+SRR15PS+SRR16PS%+10000*Penalty
SRR3step=SRR13P;%+SRR14P;%+SRR15PS%+10000*Penalty%+SRR16PS
if SRR3step<Sofarbest
    Sofarbest=SRR3step;
    bestiter=iterate;
end

iterate=iterate+1;

if bestiter==iterate-1
    save Best3step.txt theta3 SRR3step bestiter -ASCII
end

