function SRR=Minimize_Apr2013(thetain,num)

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

global LOGW_NXT_E_VCTwnxt
global LOGW_E_VCTwnxt
global SameE_VCTwnxt
global XS_EVCTwnxt_
global PartyE_VCTwnxt

global RTotDE_VCT
global X_KnotEV1


if num==4

  %Define variables needed
    
    %thetaS2=Est2(1:2,1); % do not normalize.
    %thetaQ2(1,1)=Est2(3,1);
    %     for ii=2:size(X_Knot1,2)
    %         thetaQ2(ii,1)=thetaQ2(1,1)-sum(abs(thetain(4:(ii+2),1)));
    %     end
    %thetaQ2(2:9,1)=Est2(4:11,1);
    thetaQ2=Est2(9:numel(Est2));
%     menc=-5;
%     E_VCTa(:,1)=thetain(menc+6:menc+20,1);
%     E_VCTa(:,2)=thetain(menc+21:menc+35,1);
%     E_VCTa(:,3)=thetain(menc+36:menc+50,1);
%     E_VCTt(:,1)=thetain(menc+51:menc+65,1);
%     E_VCTt(:,2)=thetain(menc+66:menc+80,1);
%     E_VCTt(:,3)=thetain(menc+81:menc+95,1);
%     gammaCT(:,1)=thetain(menc+96:menc+110,1);
%     gammaCT(:,2)=thetain(menc+111:menc+125,1);
%     gammaCT(:,3)=thetain(menc+126:menc+140,1);

    E_VCTa(:,1)=thetain(1:18,1);
    E_VCTa(:,2)=thetain(19:36,1);
    E_VCTa(:,3)=thetain(37:54,1);
    E_VCTt(:,1)=thetain(55:72,1);
    E_VCTt(:,2)=thetain(73:90,1);
    E_VCTt(:,3)=thetain(91:108,1);
    gammaCT(:,1)=thetain(109:126,1);
    gammaCT(:,2)=thetain(127:144,1);
    gammaCT(:,3)=thetain(145:162,1);
    
    %XS_EVCT2=XS_EVCT_*thetaS2;
    %XS_EVNCT2=XS_EVNCT_*thetaS2;
    if numel(thetaQ2)==3
        XQEV2=[RTotDE_VCT,RTotDE_VCT.^2,RTotDE_VCT.^3]*thetaQ2; % q_I
        XQEVCT2=XQEV2;
        %XQEVCT2(E_VContestFUL,:)=[];
        XQEVCTwnxt2=XQEVCT2;
        XQEVCTwnxt2(IND5,:)=[];
    else
        XQEV2=X_KnotEV1*thetaQ2;  % q_I
        XQEVCT2=XQEV2;
        XQEVCT2(E_VContestFUL,:)=[];
        XQEVCTwnxt2=XQEVCT2;
        XQEVCTwnxt2(IND5,:)=[];
    end
    
    
    Cutoff=[-1;3;7;100];
    SRR9=0;
    
     %%%%%%%%   F(Spending|w_I,q_I,unemp,partisan,same,contest=1)  %%%%%%%%%
     
    for i=1:3
%        TAU=[];
        TAU=find((TenureE_VCT<=Cutoff(i+1,1))&(TenureE_VCT>Cutoff(i,1))); % TAU indexes tenure;
%        LOGD_E_VCTi=[];
        LOGD_E_VCTi=LOGD_E_VCT(TAU);
%         LOGLOGD_E_VCTi=LOGLOGD_E_VCT;
%         LOGLOGD_E_VCTi(TAU)=[];
%         LOGLOGTot_E_VCTi=LOGLOGTot_E_VCT;
%         LOGLOGTot_E_VCTi(TAU)=[];
%        XQEVCTi2=[];
        XQEVCTi2=XQEVCT2(TAU);
        
%        LOGW_E_VCTi=[];
        LOGW_E_VCTi=LOGW_E_VCT(TAU);
%        SameE_VCTi=[];
        SameE_VCTi=SameE_VCT(TAU);
%        UnemploymentE_VCTi=[];
        UnemploymentE_VCTi=XS_EVCT_(TAU,1);        
%        PartisanE_VCTi=[];
        PartisanE_VCTi=XS_EVCT_(TAU,2);
        PartyE_VCTi=PartyE_VCT(TAU);        
        
        loglikeE_VDI1iP=[];
        loglikeE_VDI1iM=[];
        loglikeE_VDI2i=[];
        loglikeE_VDI4iP=[];
        loglikeE_VDI4iM=[];
        
        
        E_VCTa01iP=(max(0,1-(E_VCTa(1:6,i).^2)'*gammaCT(1:6,i).^2))^(1/2);
        E_VCTa01iM=(-1)*(max(0,1-(E_VCTa(1:6,i).^2)'*gammaCT(1:6,i).^2))^(1/2);
        loglikeE_VDI1iP(:,1)=2*log(abs(E_VCTa01iP+E_VCTa(1,i)*(LOGD_E_VCTi-E_VCTt(1,i))+E_VCTa(2,i)*(XQEVCTi2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTi-E_VCTt(3,i))+...
            E_VCTa(4,i)*(UnemploymentE_VCTi-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTi.*PartyE_VCTi-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTi-E_VCTt(6,i))));
        loglikeE_VDI1iM(:,1)=2*log(abs(E_VCTa01iM+E_VCTa(1,i)*(LOGD_E_VCTi-E_VCTt(1,i))+E_VCTa(2,i)*(XQEVCTi2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTi-E_VCTt(3,i))+...
            E_VCTa(4,i)*(UnemploymentE_VCTi-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTi.*PartyE_VCTi-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTi-E_VCTt(6,i))));
        
        loglikeE_VDI2i(:,1)=(-1/2)*(((1/gammaCT(1,i))*(LOGD_E_VCTi-E_VCTt(1,i))).^2);
        
        loglikeE_VDI3i=(-1)*length(LOGD_E_VCTi)*log(gammaCT(1,i));
        
        loglikeE_VDI4iP(:,1)=(-1)*log((E_VCTa01iP+E_VCTa(2,i)*(XQEVCTi2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTi-E_VCTt(3,i))+...
            E_VCTa(4,i)*(UnemploymentE_VCTi-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTi.*PartyE_VCTi-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTi-E_VCTt(6,i))).^2+...
        (E_VCTa(1,i)^2)*(gammaCT(1,i)^2));
        
        loglikeE_VDI4iM(:,1)=(-1)*log((E_VCTa01iM+E_VCTa(2,i)*(XQEVCTi2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTi-E_VCTt(3,i))+...
            E_VCTa(4,i)*(UnemploymentE_VCTi-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTi.*PartyE_VCTi-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTi-E_VCTt(6,i))).^2+...
            (E_VCTa(1,i)^2)*(gammaCT(1,i)^2));
        
        loglikeE_VDIi=max(sum(loglikeE_VDI1iP+loglikeE_VDI2i+loglikeE_VDI4iP)+loglikeE_VDI3i,...
            sum(loglikeE_VDI1iM+loglikeE_VDI2i+loglikeE_VDI4iM)+loglikeE_VDI3i);
        SRR9=SRR9+loglikeE_VDIi;
    end
    SRR9=(-1)*SRR9;
    
    
    
    %%%%%%%%   F(total_I|w',q_I,s,Entry,contest=1)  %%%%%%%%%
    %% Here, we want the conditional density of
    % F(total_I|w',q_I,s,tenure,Entry). Use Hermite expansion.
    % denom is the denominator that is common to most of the derivative of the
    % Loglikelihood(of total_I) function.
    
    SRR10=0;
    for i=1:3
%        TAU=[];
        TAU=find((TenureE_VCT<=Cutoff(i+1,1))&(TenureE_VCT>Cutoff(i,1))); % TAU indexes tenure;
%        LOGTotal_E_VCTi=[];
        LOGTotal_E_VCTi=LOGTotal_E_VCT(TAU);
%         LOGLOGD_E_VCTi=LOGLOGD_E_VCT;
%         LOGLOGD_E_VCTi(TAU)=[];
%         LOGLOGTot_E_VCTi=LOGLOGTot_E_VCT;
%         LOGLOGTot_E_VCTi(TAU)=[];
%        XQEVCTi2=[];
        XQEVCTi2=XQEVCT2(TAU);
%        LOGW_E_VCTi=[];
        LOGW_E_VCTi=LOGW_E_VCT(TAU);
%        PartyE_VCTi=[];
        SameE_VCTi=SameE_VCT(TAU);
%        UnemploymentE_VCTi=[];
        UnemploymentE_VCTi=XS_EVCT_(TAU,1);
%        pctWhiteE_VCTi=[];
        PartisanE_VCTi=XS_EVCT_(TAU,2);
        PartyE_VCTi=PartyE_VCT(TAU);        
                
        loglikeE_VTotI1iP=[];
        loglikeE_VTotI1iM=[];
        loglikeE_VTotI2i=[];
        loglikeE_VTotI4iP=[];
        loglikeE_VTotI4iM=[];
        
        E_VCTa02iP=(max(0,1-(E_VCTa(7:12,i).^2)'*gammaCT(7:12,i).^2))^(1/2);%?
        E_VCTa02iM=(-1)*(max(0,1-(E_VCTa(7:12,i).^2)'*gammaCT(7:12,i).^2))^(1/2);
        
        loglikeE_VTotI1iP(:,1)=2*log(abs(E_VCTa02iP+E_VCTa(7,i)*(LOGTotal_E_VCTi-E_VCTt(7,i))+E_VCTa(8,i)*(XQEVCTi2-E_VCTt(8,i))+E_VCTa(9,i)*(LOGW_E_VCTi-E_VCTt(9,i))+...
            E_VCTa(10,i)*(UnemploymentE_VCTi-E_VCTt(10,i))+E_VCTa(11,i)*(PartisanE_VCTi.*PartyE_VCTi-E_VCTt(11,i))+E_VCTa(12,i)*(SameE_VCTi-E_VCTt(12,i))));
        
        loglikeE_VTotI1iM(:,1)=2*log(abs(E_VCTa02iM+E_VCTa(7,i)*(LOGTotal_E_VCTi-E_VCTt(7,i))+E_VCTa(8,i)*(XQEVCTi2-E_VCTt(8,i))+E_VCTa(9,i)*(LOGW_E_VCTi-E_VCTt(9,i))+...
            E_VCTa(10,i)*(UnemploymentE_VCTi-E_VCTt(10,i))+E_VCTa(11,i)*(PartisanE_VCTi.*PartyE_VCTi-E_VCTt(11,i))+E_VCTa(12,i)*(SameE_VCTi-E_VCTt(12,i))));
        
        loglikeE_VTotI2i(:,1)=(-1/2)*(((1/gammaCT(7,i))*(LOGTotal_E_VCTi-E_VCTt(7,i))).^2);
        
        loglikeE_VTotI3i=(-1)*length(LOGTotal_E_VCTi)*log(gammaCT(7,i));
        
        loglikeE_VTotI4iP(:,1)=(-1)*log((E_VCTa02iP+E_VCTa(8,i)*(XQEVCTi2-E_VCTt(8,i))+E_VCTa(9,i)*(LOGW_E_VCTi-E_VCTt(9,i))+...
            E_VCTa(10,i)*(UnemploymentE_VCTi-E_VCTt(10,i))+E_VCTa(11,i)*(PartisanE_VCTi.*PartyE_VCTi-E_VCTt(11,i))+E_VCTa(12,i)*(SameE_VCTi-E_VCTt(12,i))).^2+...
        (E_VCTa(7,i)^2)*(gammaCT(7,i)^2));
        
        loglikeE_VTotI4iM(:,1)=(-1)*log((E_VCTa02iM+E_VCTa(8,i)*(XQEVCTi2-E_VCTt(8,i))+E_VCTa(9,i)*(LOGW_E_VCTi-E_VCTt(9,i))+...
            E_VCTa(10,i)*(UnemploymentE_VCTi-E_VCTt(10,i))+E_VCTa(11,i)*(PartisanE_VCTi.*PartyE_VCTi-E_VCTt(11,i))+E_VCTa(12,i)*(SameE_VCTi-E_VCTt(12,i))).^2+...
        (E_VCTa(7,i)^2)*(gammaCT(7,i)^2));
        
        loglikeE_VTotIi=max(sum(loglikeE_VTotI1iP+loglikeE_VTotI2i+loglikeE_VTotI4iP)+loglikeE_VTotI3i,...
            sum(loglikeE_VTotI1iM+loglikeE_VTotI2i+loglikeE_VTotI4iM)+loglikeE_VTotI3i);
        SRR10=SRR10+loglikeE_VTotIi;
    end
    SRR10=(-1)*SRR10;
    
    
    SRR11=0;
    for i=1:3
 %       TAU=[];
        TAU=find((TenureE_VCTwnxt<=Cutoff(i+1,1))&(TenureE_VCTwnxt>Cutoff(i,1))); % TAU indexes tenure;
 %       LOGW_NXT_E_VCTwnxti=[];
        LOGW_NXT_E_VCTwnxti=LOGW_NXT_E_VCTwnxt(TAU);
%         LOGLOGD_E_VCTwnxti=LOGLOGD_E_VCTwnxt;
%         LOGLOGD_E_VCTwnxti(TAU)=[];
%         LOGLOGTot_E_VCTwnxti=LOGLOGTot_E_VCTwnxt;
%         LOGLOGTot_E_VCTwnxti(TAU)=[];
%        XQEVCTwnxti2=[];
        XQEVCTwnxti2=XQEVCTwnxt2(TAU);
%        LOGW_E_VCTwnxti=[];
        LOGW_E_VCTwnxti=LOGW_E_VCTwnxt(TAU);
%        PartyE_VCTwnxti=[];
        SameE_VCTwnxti=SameE_VCTwnxt(TAU);
%        UnemploymentE_VCTwnxti=[];
        UnemploymentE_VCTwnxti=XS_EVCTwnxt_(TAU,1);
%        pctWhiteE_VCTwnxti=[];
        PartisanE_VCTwnxti=XS_EVCTwnxt_(TAU,2);
        PartyE_VCTwnxti=PartyE_VCTwnxt(TAU);
        loglikeE_VWnxtnxtI1iP=[];
        loglikeE_VWnxtnxtI1iM=[];
        loglikeE_VWnxtnxtI2i=[];
        loglikeE_VWnxtnxtI4iP=[];
        loglikeE_VWnxtnxtI4iM=[];
        
        E_VCTa03iP=(max(0,1-(E_VCTa(13:18,i).^2)'*gammaCT(13:18,i).^2))^(1/2);
        E_VCTa03iM=(-1)*(max(0,1-(E_VCTa(13:18,i).^2)'*gammaCT(13:18,i).^2))^(1/2);

        
        loglikeE_VWnxtnxtI1iP(:,1)=2*log(abs(E_VCTa03iP+E_VCTa(13,i)*(LOGW_NXT_E_VCTwnxti-E_VCTt(13,i))+E_VCTa(14,i)*(XQEVCTwnxti2-E_VCTt(14,i))+E_VCTa(15,i)*(LOGW_E_VCTwnxti-E_VCTt(15,i))+...
            E_VCTa(16,i)*(UnemploymentE_VCTwnxti-E_VCTt(16,i))+E_VCTa(17,i)*(PartisanE_VCTwnxti.*PartyE_VCTwnxti-E_VCTt(17,i))+E_VCTa(18,i)*(SameE_VCTwnxti-E_VCTt(18,i))));
        
        loglikeE_VWnxtnxtI1iM(:,1)=2*log(abs(E_VCTa03iM+E_VCTa(13,i)*(LOGW_NXT_E_VCTwnxti-E_VCTt(13,i))+E_VCTa(14,i)*(XQEVCTwnxti2-E_VCTt(14,i))+E_VCTa(15,i)*(LOGW_E_VCTwnxti-E_VCTt(15,i))+...
            E_VCTa(16,i)*(UnemploymentE_VCTwnxti-E_VCTt(16,i))+E_VCTa(17,i)*(PartisanE_VCTwnxti.*PartyE_VCTwnxti-E_VCTt(17,i))+E_VCTa(18,i)*(SameE_VCTwnxti-E_VCTt(18,i))));
        
        loglikeE_VWnxtnxtI2i(:,1)=(-1/2)*(((1/gammaCT(13,i))*(LOGW_NXT_E_VCTwnxti-E_VCTt(13,i))).^2);
        
        loglikeE_VWnxtnxtI3i=(-1)*length(LOGW_NXT_E_VCTwnxti)*log(gammaCT(13,i));
        
        loglikeE_VWnxtnxtI4iP(:,1)=(-1)*log((E_VCTa03iP+E_VCTa(14,i)*(XQEVCTwnxti2-E_VCTt(14,i))+E_VCTa(15,i)*(LOGW_E_VCTwnxti-E_VCTt(15,i))+...
            E_VCTa(16,i)*(UnemploymentE_VCTwnxti-E_VCTt(16,i))+E_VCTa(17,i)*(PartisanE_VCTwnxti.*PartyE_VCTwnxti-E_VCTt(17,i))+E_VCTa(18,i)*(SameE_VCTwnxti-E_VCTt(18,i))).^2+...
        (E_VCTa(13,i)^2)*(gammaCT(13,i)^2));
        
        loglikeE_VWnxtnxtI4iM(:,1)=(-1)*log((E_VCTa03iM+E_VCTa(14,i)*(XQEVCTwnxti2-E_VCTt(14,i))+E_VCTa(15,i)*(LOGW_E_VCTwnxti-E_VCTt(15,i))+...
            E_VCTa(16,i)*(UnemploymentE_VCTwnxti-E_VCTt(16,i))+E_VCTa(17,i)*(PartisanE_VCTwnxti.*PartyE_VCTwnxti-E_VCTt(17,i))+E_VCTa(18,i)*(SameE_VCTwnxti-E_VCTt(18,i))).^2+...
            (E_VCTa(13,i)^2)*(gammaCT(13,i)^2));
        
        loglikeE_VWnxtnxtIi=max(sum(loglikeE_VWnxtnxtI1iP+loglikeE_VWnxtnxtI2i+loglikeE_VWnxtnxtI4iP)+loglikeE_VWnxtnxtI3i,...
            sum(loglikeE_VWnxtnxtI1iM+loglikeE_VWnxtnxtI2i+loglikeE_VWnxtnxtI4iM)+loglikeE_VWnxtnxtI3i);
        SRR11=SRR11+loglikeE_VWnxtnxtIi;
    end
    SRR11=(-1)*SRR11;
    
end
% if num==5
%     
%     %%%%%%%%   gamma_I=gamma_I(w',q_I,s,Entry,tenure|contest=0)  %%%%%%%%%
%     %% Here, we want the (deterministic) function of spending|contest=0 as a
%     %% function of w',q_I,s,Entry,tenure
%     Xgamma_INCT=[ones(size(LOGD_E_VNCT,1),1),LOGW_E_VNCT,PartyE_VNCT.*XS_EVNCT_(:,1),PartyE_VNCT.*XS_EVNCT_(:,2),TenureE_VNCT,XQEVNCT2,...
%         LOGW_E_VNCT.^2,PartyE_VNCT.*(XS_EVNCT_(:,1)).^2,PartyE_VNCT.*(XS_EVNCT_(:,2)).^2,TenureE_VNCT.^2,XQEVNCT2.^2];
%     SRR12=sum((LOGD_E_VNCT-Xgamma_INCT*gammaNCT).^2);
%     
%     
%     %%%%%%%   total_I=tot_I(w',q_I,s,Entry|contest=0)  %%%%%%%%%
%     % Here, we want the deterministic function of fund-raising conditional on
%     % being uncontested as a function of w',q_I,s,tenure,Entry
%     
%     Xt_INCT=[ones(size(LOGD_E_VNCT,1),1),LOGW_E_VNCT,PartyE_VNCT.*XS_EVNCT_(:,1),PartyE_VNCT.*XS_EVNCT_(:,2),TenureE_VNCT,XQEVNCT2,...
%         LOGW_E_VNCT.^2,PartyE_VNCT.*(XS_EVNCT_(:,1)).^2,PartyE_VNCT.*(XS_EVNCT_(:,2)).^2,TenureE_VNCT.^2,XQEVNCT2.^2];
%     SRR13=sum((LOGTotal_E_VNCT-Xt_INCT*tNCT).^2);
%     
%     
%     %%%%%%%%   w'_(t+1)_I=savings_I(w',q_I,s,Entry|contest=0)  %%%%%%%%%
%     %% Here, we want the deterministic function of savings conditional on being
%     %% uncontested as a function of (w',q_I,s,tenure,Entry.
%     
%     Xw_INCT=[ones(size(LOGD_E_VNCT,1),1),LOGW_E_VNCT,PartyE_VNCT.*XS_EVNCT_(:,1),PartyE_VNCT.*XS_EVNCT_(:,2),TenureE_VNCT,XQEVNCT2,...
%         LOGW_E_VNCT.^2,PartyE_VNCT.*(XS_EVNCT_(:,1)).^2,PartyE_VNCT.*(XS_EVNCT_(:,2)).^2,TenureE_VNCT.^2,XQEVNCT2.^2];
%     SRR14=sum((LOGW_NXT_E_VNCT-Xw_INCT*wNCT).^2);
%     
%     
% end


if num==4
    
    if (((E_VCTa(1:6,1).^2)'*gammaCT(1:6,1).^2)>=1)||(((E_VCTa(1:6,2).^2)'*gammaCT(1:6,2).^2)>=1)||(((E_VCTa(1:6,3).^2)'*gammaCT(1:6,3).^2)>=1)
        for i=1:3
            SRR9=SRR9+10000000*max(0,((E_VCTa(1:6,i).^2)'*gammaCT(1:6,i).^2-1));
        end
    end
    if (((E_VCTa(7:12,1).^2)'*gammaCT(7:12,1).^2)>=1)||(((E_VCTa(7:12,2).^2)'*gammaCT(7:12,2).^2)>=1)||(((E_VCTa(7:12,3).^2)'*gammaCT(7:12,3).^2)>=1)
        for i=1:3
            SRR10=SRR10+10000000*max(0,((E_VCTa(7:12,i).^2)'*gammaCT(7:12,i).^2-1));
        end
    end
    if (((E_VCTa(13:18,1).^2)'*gammaCT(13:18,1).^2)>=1)||(((E_VCTa(13:18,2).^2)'*gammaCT(13:18,2).^2)>=1)||(((E_VCTa(13:18,3).^2)'*gammaCT(13:18,3).^2)>=1)
        for i=1:3
            SRR11=SRR11+10000000*max(0,((E_VCTa(13:18,i).^2)'*gammaCT(13:18,i).^2-1));
        end
    end
end


SRR=SRR9+SRR10+SRR11;


