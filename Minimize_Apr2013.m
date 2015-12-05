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
global PresdumE_VCT
global MidtermE_VCT

global LOGW_NXT_E_VCTwnxt
global LOGW_E_VCTwnxt
global SameE_VCTwnxt
global XS_EVCTwnxt_
global PartyE_VCTwnxt
global PresdumE_VCTwnxt
global MidtermE_VCTwnxt

global RTotDE_VCT
global X_KnotEV1


%Define variables common across all cases
Cutoff=[-1;3;7;100];

thetaQ2=Est2(9:numel(Est2));

if num==2
E_VCTa(:,1)=thetain(1:8,1);
E_VCTa(:,2)=thetain(9:16,1);
E_VCTa(:,3)=thetain(17:24,1);

E_VCTt(:,1)=thetain(25:32,1);
E_VCTt(:,2)=thetain(33:40,1);
E_VCTt(:,3)=thetain(41:48,1);

gammaCT(:,1)=thetain(49:56,1);
gammaCT(:,2)=thetain(57:64,1);
gammaCT(:,3)=thetain(65:72,1);
end

if num>10
E_VCTa(:,1)=thetain(1:8,1);
E_VCTt(:,1)=thetain(9:16,1);
gammaCT(:,1)=thetain(17:24,1);
end

if numel(thetaQ2)==3
    XQEV2=[RTotDE_VCT,RTotDE_VCT.^2,RTotDE_VCT.^3]*thetaQ2; % q_I
    XQEVCT2=XQEV2;
    %XQEVCT2(E_VContestFUL,:)=[]; Already VCT
    XQEVCTwnxt2=XQEVCT2;
    XQEVCTwnxt2(IND5,:)=[];
else
    XQEV2=X_KnotEV1*thetaQ2;  % q_I
    XQEVCT2=XQEV2;
    XQEVCT2(E_VContestFUL,:)=[];
    XQEVCTwnxt2=XQEVCT2;
    XQEVCTwnxt2(IND5,:)=[];
end



if num==11
    SRR9=0;
    
    %%%%%%%%   F(Spending|w_I,q_I,unemp,partisan,same,contest=1)  %%%%%%%%%
    
    i=1;
        
        TAU=find((TenureE_VCT<=Cutoff(2,1))&(TenureE_VCT>Cutoff(1,1))); % TAU indexes tenure;
        
        LOGD_E_VCTi=LOGD_E_VCT(TAU);
        XQEVCTi2=XQEVCT2(TAU);
        LOGW_E_VCTi=LOGW_E_VCT(TAU);
        
        SameE_VCTi=SameE_VCT(TAU);
        
        UnemploymentE_VCTi=XS_EVCT_(TAU,1);
        
        PartisanE_VCTi=XS_EVCT_(TAU,2);
        PartyE_VCTi=PartyE_VCT(TAU);
        
        PresdumE_VCTi=PresdumE_VCT(TAU);
        MidtermE_VCTi=MidtermE_VCT(TAU);
        
        loglikeE_VDI1iP=[];
        loglikeE_VDI1iM=[];
        loglikeE_VDI2i=[];
        loglikeE_VDI4iP=[];
        loglikeE_VDI4iM=[];
        
        
        E_VCTa01iP=(max(0,1-(E_VCTa(1:8,i).^2)'*gammaCT(1:8,i).^2))^(1/2);
        E_VCTa01iM=(-1)*(max(0,1-(E_VCTa(1:8,i).^2)'*gammaCT(1:8,i).^2))^(1/2);
        loglikeE_VDI1iP(:,1)=2*log(abs(E_VCTa01iP+E_VCTa(1,i)*(LOGD_E_VCTi-E_VCTt(1,i))+E_VCTa(2,i)*(XQEVCTi2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTi-E_VCTt(3,i))+...
            E_VCTa(4,i)*(UnemploymentE_VCTi-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTi.*PartyE_VCTi-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTi-E_VCTt(6,i))+...
            E_VCTa(7,i)*(PresdumE_VCTi-E_VCTt(7,i))+E_VCTa(8,i)*(MidtermE_VCTi-E_VCTt(8,i))));
        loglikeE_VDI1iM(:,1)=2*log(abs(E_VCTa01iM+E_VCTa(1,i)*(LOGD_E_VCTi-E_VCTt(1,i))+E_VCTa(2,i)*(XQEVCTi2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTi-E_VCTt(3,i))+...
            E_VCTa(4,i)*(UnemploymentE_VCTi-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTi.*PartyE_VCTi-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTi-E_VCTt(6,i))+...
            E_VCTa(7,i)*(PresdumE_VCTi-E_VCTt(7,i))+E_VCTa(8,i)*(MidtermE_VCTi-E_VCTt(8,i))));
        
        loglikeE_VDI2i(:,1)=(-1/2)*(((1/gammaCT(1,i))*(LOGD_E_VCTi-E_VCTt(1,i))).^2);
        
        loglikeE_VDI3i=(-1)*length(LOGD_E_VCTi)*log(gammaCT(1,i));
        
        loglikeE_VDI4iP(:,1)=(-1)*log((E_VCTa01iP+E_VCTa(2,i)*(XQEVCTi2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTi-E_VCTt(3,i))+...
            E_VCTa(4,i)*(UnemploymentE_VCTi-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTi.*PartyE_VCTi-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTi-E_VCTt(6,i))+...
            E_VCTa(7,i)*(PresdumE_VCTi-E_VCTt(7,i))+E_VCTa(8,i)*(MidtermE_VCTi-E_VCTt(8,i))).^2+...
            (E_VCTa(1,i)^2)*(gammaCT(1,i)^2));
        
        loglikeE_VDI4iM(:,1)=(-1)*log((E_VCTa01iM+E_VCTa(2,i)*(XQEVCTi2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTi-E_VCTt(3,i))+...
            E_VCTa(4,i)*(UnemploymentE_VCTi-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTi.*PartyE_VCTi-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTi-E_VCTt(6,i))+...
            E_VCTa(7,i)*(PresdumE_VCTi-E_VCTt(7,i))+E_VCTa(8,i)*(MidtermE_VCTi-E_VCTt(8,i))).^2+...
            (E_VCTa(1,i)^2)*(gammaCT(1,i)^2));
        
        loglikeE_VDIi=max(sum(loglikeE_VDI1iP+loglikeE_VDI2i+loglikeE_VDI4iP)+loglikeE_VDI3i,...
            sum(loglikeE_VDI1iM+loglikeE_VDI2i+loglikeE_VDI4iM)+loglikeE_VDI3i);
        SRR9=SRR9+loglikeE_VDIi;
    
    SRR9=(-1)*SRR9;
    
    
    if (((E_VCTa(1:8,1).^2)'*gammaCT(1:8,1).^2)>=1)

            SRR9=SRR9+10000000*max(0,((E_VCTa(1:8,i).^2)'*gammaCT(1:8,i).^2-1));

    end
    
    SRR=SRR9;
end

if num==12
    SRR9=0;
    
    %%%%%%%%   F(Spending|w_I,q_I,unemp,partisan,same,contest=1)  %%%%%%%%%
    
    i=1;
        
        TAU=find((TenureE_VCT<=Cutoff(3,1))&(TenureE_VCT>Cutoff(2,1))); % TAU indexes tenure;
        
        LOGD_E_VCTi=LOGD_E_VCT(TAU);
        XQEVCTi2=XQEVCT2(TAU);
        LOGW_E_VCTi=LOGW_E_VCT(TAU);
        
        SameE_VCTi=SameE_VCT(TAU);
        
        UnemploymentE_VCTi=XS_EVCT_(TAU,1);
        
        PartisanE_VCTi=XS_EVCT_(TAU,2);
        PartyE_VCTi=PartyE_VCT(TAU);
        
        PresdumE_VCTi=PresdumE_VCT(TAU);
        MidtermE_VCTi=MidtermE_VCT(TAU);
        
        loglikeE_VDI1iP=[];
        loglikeE_VDI1iM=[];
        loglikeE_VDI2i=[];
        loglikeE_VDI4iP=[];
        loglikeE_VDI4iM=[];
        
        
        E_VCTa01iP=(max(0,1-(E_VCTa(1:8,i).^2)'*gammaCT(1:8,i).^2))^(1/2);
        E_VCTa01iM=(-1)*(max(0,1-(E_VCTa(1:8,i).^2)'*gammaCT(1:8,i).^2))^(1/2);
        loglikeE_VDI1iP(:,1)=2*log(abs(E_VCTa01iP+E_VCTa(1,i)*(LOGD_E_VCTi-E_VCTt(1,i))+E_VCTa(2,i)*(XQEVCTi2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTi-E_VCTt(3,i))+...
            E_VCTa(4,i)*(UnemploymentE_VCTi-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTi.*PartyE_VCTi-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTi-E_VCTt(6,i))+...
            E_VCTa(7,i)*(PresdumE_VCTi-E_VCTt(7,i))+E_VCTa(8,i)*(MidtermE_VCTi-E_VCTt(8,i))));
        loglikeE_VDI1iM(:,1)=2*log(abs(E_VCTa01iM+E_VCTa(1,i)*(LOGD_E_VCTi-E_VCTt(1,i))+E_VCTa(2,i)*(XQEVCTi2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTi-E_VCTt(3,i))+...
            E_VCTa(4,i)*(UnemploymentE_VCTi-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTi.*PartyE_VCTi-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTi-E_VCTt(6,i))+...
            E_VCTa(7,i)*(PresdumE_VCTi-E_VCTt(7,i))+E_VCTa(8,i)*(MidtermE_VCTi-E_VCTt(8,i))));
        
        loglikeE_VDI2i(:,1)=(-1/2)*(((1/gammaCT(1,i))*(LOGD_E_VCTi-E_VCTt(1,i))).^2);
        
        loglikeE_VDI3i=(-1)*length(LOGD_E_VCTi)*log(gammaCT(1,i));
        
        loglikeE_VDI4iP(:,1)=(-1)*log((E_VCTa01iP+E_VCTa(2,i)*(XQEVCTi2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTi-E_VCTt(3,i))+...
            E_VCTa(4,i)*(UnemploymentE_VCTi-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTi.*PartyE_VCTi-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTi-E_VCTt(6,i))+...
            E_VCTa(7,i)*(PresdumE_VCTi-E_VCTt(7,i))+E_VCTa(8,i)*(MidtermE_VCTi-E_VCTt(8,i))).^2+...
            (E_VCTa(1,i)^2)*(gammaCT(1,i)^2));
        
        loglikeE_VDI4iM(:,1)=(-1)*log((E_VCTa01iM+E_VCTa(2,i)*(XQEVCTi2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTi-E_VCTt(3,i))+...
            E_VCTa(4,i)*(UnemploymentE_VCTi-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTi.*PartyE_VCTi-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTi-E_VCTt(6,i))+...
            E_VCTa(7,i)*(PresdumE_VCTi-E_VCTt(7,i))+E_VCTa(8,i)*(MidtermE_VCTi-E_VCTt(8,i))).^2+...
            (E_VCTa(1,i)^2)*(gammaCT(1,i)^2));
        
        loglikeE_VDIi=max(sum(loglikeE_VDI1iP+loglikeE_VDI2i+loglikeE_VDI4iP)+loglikeE_VDI3i,...
            sum(loglikeE_VDI1iM+loglikeE_VDI2i+loglikeE_VDI4iM)+loglikeE_VDI3i);
        SRR9=SRR9+loglikeE_VDIi;

    SRR9=(-1)*SRR9;
    
    
    if (((E_VCTa(1:8,1).^2)'*gammaCT(1:8,1).^2)>=1)

            SRR9=SRR9+10000000*max(0,((E_VCTa(1:8,i).^2)'*gammaCT(1:8,i).^2-1));

    end
    
    SRR=SRR9;
end

if num==13
    SRR9=0;
    
    %%%%%%%%   F(Spending|w_I,q_I,unemp,partisan,same,contest=1)  %%%%%%%%%
    
     i=1;
        
        TAU=find((TenureE_VCT<=Cutoff(4,1))&(TenureE_VCT>Cutoff(3,1))); % TAU indexes tenure;
        
        LOGD_E_VCTi=LOGD_E_VCT(TAU);
        XQEVCTi2=XQEVCT2(TAU);
        LOGW_E_VCTi=LOGW_E_VCT(TAU);
        
        SameE_VCTi=SameE_VCT(TAU);
        
        UnemploymentE_VCTi=XS_EVCT_(TAU,1);
        
        PartisanE_VCTi=XS_EVCT_(TAU,2);
        PartyE_VCTi=PartyE_VCT(TAU);
        
        PresdumE_VCTi=PresdumE_VCT(TAU);
        MidtermE_VCTi=MidtermE_VCT(TAU);
        
        loglikeE_VDI1iP=[];
        loglikeE_VDI1iM=[];
        loglikeE_VDI2i=[];
        loglikeE_VDI4iP=[];
        loglikeE_VDI4iM=[];
        
        
        E_VCTa01iP=(max(0,1-(E_VCTa(1:8,i).^2)'*gammaCT(1:8,i).^2))^(1/2);
        E_VCTa01iM=(-1)*(max(0,1-(E_VCTa(1:8,i).^2)'*gammaCT(1:8,i).^2))^(1/2);
        loglikeE_VDI1iP(:,1)=2*log(abs(E_VCTa01iP+E_VCTa(1,i)*(LOGD_E_VCTi-E_VCTt(1,i))+E_VCTa(2,i)*(XQEVCTi2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTi-E_VCTt(3,i))+...
            E_VCTa(4,i)*(UnemploymentE_VCTi-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTi.*PartyE_VCTi-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTi-E_VCTt(6,i))+...
            E_VCTa(7,i)*(PresdumE_VCTi-E_VCTt(7,i))+E_VCTa(8,i)*(MidtermE_VCTi-E_VCTt(8,i))));
        loglikeE_VDI1iM(:,1)=2*log(abs(E_VCTa01iM+E_VCTa(1,i)*(LOGD_E_VCTi-E_VCTt(1,i))+E_VCTa(2,i)*(XQEVCTi2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTi-E_VCTt(3,i))+...
            E_VCTa(4,i)*(UnemploymentE_VCTi-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTi.*PartyE_VCTi-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTi-E_VCTt(6,i))+...
            E_VCTa(7,i)*(PresdumE_VCTi-E_VCTt(7,i))+E_VCTa(8,i)*(MidtermE_VCTi-E_VCTt(8,i))));
        
        loglikeE_VDI2i(:,1)=(-1/2)*(((1/gammaCT(1,i))*(LOGD_E_VCTi-E_VCTt(1,i))).^2);
        
        loglikeE_VDI3i=(-1)*length(LOGD_E_VCTi)*log(gammaCT(1,i));
        
        loglikeE_VDI4iP(:,1)=(-1)*log((E_VCTa01iP+E_VCTa(2,i)*(XQEVCTi2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTi-E_VCTt(3,i))+...
            E_VCTa(4,i)*(UnemploymentE_VCTi-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTi.*PartyE_VCTi-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTi-E_VCTt(6,i))+...
            E_VCTa(7,i)*(PresdumE_VCTi-E_VCTt(7,i))+E_VCTa(8,i)*(MidtermE_VCTi-E_VCTt(8,i))).^2+...
            (E_VCTa(1,i)^2)*(gammaCT(1,i)^2));
        
        loglikeE_VDI4iM(:,1)=(-1)*log((E_VCTa01iM+E_VCTa(2,i)*(XQEVCTi2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTi-E_VCTt(3,i))+...
            E_VCTa(4,i)*(UnemploymentE_VCTi-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTi.*PartyE_VCTi-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTi-E_VCTt(6,i))+...
            E_VCTa(7,i)*(PresdumE_VCTi-E_VCTt(7,i))+E_VCTa(8,i)*(MidtermE_VCTi-E_VCTt(8,i))).^2+...
            (E_VCTa(1,i)^2)*(gammaCT(1,i)^2));
        
        loglikeE_VDIi=max(sum(loglikeE_VDI1iP+loglikeE_VDI2i+loglikeE_VDI4iP)+loglikeE_VDI3i,...
            sum(loglikeE_VDI1iM+loglikeE_VDI2i+loglikeE_VDI4iM)+loglikeE_VDI3i);
        SRR9=SRR9+loglikeE_VDIi;
   
    SRR9=(-1)*SRR9;
    
    
    if (((E_VCTa(1:8,1).^2)'*gammaCT(1:8,1).^2)>=1)
            SRR9=SRR9+10000000*max(0,((E_VCTa(1:8,i).^2)'*gammaCT(1:8,i).^2-1));
    end
    
    SRR=SRR9;
end


if num==2
    
    
    
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
        
        PresdumE_VCTi=PresdumE_VCT(TAU);
        MidtermE_VCTi=MidtermE_VCT(TAU);
        
        loglikeE_VTotI1iP=[];
        loglikeE_VTotI1iM=[];
        loglikeE_VTotI2i=[];
        loglikeE_VTotI4iP=[];
        loglikeE_VTotI4iM=[];
        
        E_VCTa02iP=(max(0,1-(E_VCTa(1:8,i).^2)'*gammaCT(1:8,i).^2))^(1/2);%?
        E_VCTa02iM=(-1)*(max(0,1-(E_VCTa(1:8,i).^2)'*gammaCT(1:8,i).^2))^(1/2);
        
        loglikeE_VTotI1iP(:,1)=2*log(abs(E_VCTa02iP+E_VCTa(1,i)*(LOGTotal_E_VCTi-E_VCTt(1,i))+E_VCTa(2,i)*(XQEVCTi2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTi-E_VCTt(3,i))+...
            E_VCTa(4,i)*(UnemploymentE_VCTi-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTi.*PartyE_VCTi-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTi-E_VCTt(6,i))+...
            E_VCTa(7,i)*(PresdumE_VCTi-E_VCTt(7,i))+E_VCTa(8,i)*(MidtermE_VCTi-E_VCTt(8,i))));
        
        loglikeE_VTotI1iM(:,1)=2*log(abs(E_VCTa02iM+E_VCTa(1,i)*(LOGTotal_E_VCTi-E_VCTt(1,i))+E_VCTa(2,i)*(XQEVCTi2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTi-E_VCTt(3,i))+...
            E_VCTa(4,i)*(UnemploymentE_VCTi-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTi.*PartyE_VCTi-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTi-E_VCTt(6,i))+...
            E_VCTa(7,i)*(PresdumE_VCTi-E_VCTt(7,i))+E_VCTa(8,i)*(MidtermE_VCTi-E_VCTt(8,i))));
        
        loglikeE_VTotI2i(:,1)=(-1/2)*(((1/gammaCT(1,i))*(LOGTotal_E_VCTi-E_VCTt(1,i))).^2);
        
        loglikeE_VTotI3i=(-1)*length(LOGTotal_E_VCTi)*log(gammaCT(1,i));
        
        loglikeE_VTotI4iP(:,1)=(-1)*log((E_VCTa02iP+E_VCTa(2,i)*(XQEVCTi2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTi-E_VCTt(3,i))+...
            E_VCTa(4,i)*(UnemploymentE_VCTi-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTi.*PartyE_VCTi-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTi-E_VCTt(6,i))+...
            E_VCTa(7,i)*(PresdumE_VCTi-E_VCTt(7,i))+E_VCTa(8,i)*(MidtermE_VCTi-E_VCTt(8,i))).^2+...
            (E_VCTa(1,i)^2)*(gammaCT(1,i)^2));
        
        loglikeE_VTotI4iM(:,1)=(-1)*log((E_VCTa02iM+E_VCTa(2,i)*(XQEVCTi2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTi-E_VCTt(3,i))+...
            E_VCTa(4,i)*(UnemploymentE_VCTi-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTi.*PartyE_VCTi-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTi-E_VCTt(6,i))+...
            E_VCTa(7,i)*(PresdumE_VCTi-E_VCTt(7,i))+E_VCTa(8,i)*(MidtermE_VCTi-E_VCTt(8,i))).^2+...
            (E_VCTa(1,i)^2)*(gammaCT(1,i)^2));
        
        loglikeE_VTotIi=max(sum(loglikeE_VTotI1iP+loglikeE_VTotI2i+loglikeE_VTotI4iP)+loglikeE_VTotI3i,...
            sum(loglikeE_VTotI1iM+loglikeE_VTotI2i+loglikeE_VTotI4iM)+loglikeE_VTotI3i);
        SRR10=SRR10+loglikeE_VTotIi;
    end
    SRR10=(-1)*SRR10;
    
    if (((E_VCTa(1:8,1).^2)'*gammaCT(1:8,1).^2)>=1)||(((E_VCTa(1:8,2).^2)'*gammaCT(1:8,2).^2)>=1)||(((E_VCTa(1:8,3).^2)'*gammaCT(1:8,3).^2)>=1)
        for i=1:3
            SRR10=SRR10+10000000*max(0,((E_VCTa(1:8,i).^2)'*gammaCT(1:8,i).^2-1));
        end
    end
    
    SRR=SRR10;
end


% 
% if num==3
% 
%     
%     SRR11=0;
%     for i=1:3
%         %       TAU=[];
%         TAU=find((TenureE_VCTwnxt<=Cutoff(i+1,1))&(TenureE_VCTwnxt>Cutoff(i,1))); % TAU indexes tenure;
%         %       LOGW_NXT_E_VCTwnxti=[];
%         LOGW_NXT_E_VCTwnxti=LOGW_NXT_E_VCTwnxt(TAU);
%         %         LOGLOGD_E_VCTwnxti=LOGLOGD_E_VCTwnxt;
%         %         LOGLOGD_E_VCTwnxti(TAU)=[];
%         %         LOGLOGTot_E_VCTwnxti=LOGLOGTot_E_VCTwnxt;
%         %         LOGLOGTot_E_VCTwnxti(TAU)=[];
%         %        XQEVCTwnxti2=[];
%         XQEVCTwnxti2=XQEVCTwnxt2(TAU);
%         %        LOGW_E_VCTwnxti=[];
%         LOGW_E_VCTwnxti=LOGW_E_VCTwnxt(TAU);
%         %        PartyE_VCTwnxti=[];
%         SameE_VCTwnxti=SameE_VCTwnxt(TAU);
%         %        UnemploymentE_VCTwnxti=[];
%         UnemploymentE_VCTwnxti=XS_EVCTwnxt_(TAU,1);
%         %        pctWhiteE_VCTwnxti=[];
%         PartisanE_VCTwnxti=XS_EVCTwnxt_(TAU,2);
%         PartyE_VCTwnxti=PartyE_VCTwnxt(TAU);
%         
%         PresdumE_VCTwnxti=PresdumE_VCTwnxt(TAU);
%         MidtermE_VCTwnxti=MidtermE_VCTwnxt(TAU);
%         loglikeE_VWnxtnxtI1iP=[];
%         loglikeE_VWnxtnxtI1iM=[];
%         loglikeE_VWnxtnxtI2i=[];
%         loglikeE_VWnxtnxtI4iP=[];
%         loglikeE_VWnxtnxtI4iM=[];
%         
%         E_VCTa03iP=(max(0,1-(E_VCTa(1:8,i).^2)'*gammaCT(1:8,i).^2))^(1/2);
%         E_VCTa03iM=(-1)*(max(0,1-(E_VCTa(1:8,i).^2)'*gammaCT(1:8,i).^2))^(1/2);
%         
%         
%         loglikeE_VWnxtnxtI1iP(:,1)=2*log(abs(E_VCTa03iP+E_VCTa(1,i)*(LOGW_NXT_E_VCTwnxti-E_VCTt(1,i))+E_VCTa(2,i)*(XQEVCTwnxti2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTwnxti-E_VCTt(3,i))+...
%             E_VCTa(4,i)*(UnemploymentE_VCTwnxti-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTwnxti.*PartyE_VCTwnxti-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTwnxti-E_VCTt(6,i))+...
%             E_VCTa(7,i)*(PresdumE_VCTwnxti-E_VCTt(7,i))+E_VCTa(8,i)*(MidtermE_VCTwnxti-E_VCTt(8,i))));
%         
%         loglikeE_VWnxtnxtI1iM(:,1)=2*log(abs(E_VCTa03iM+E_VCTa(1,i)*(LOGW_NXT_E_VCTwnxti-E_VCTt(1,i))+E_VCTa(2,i)*(XQEVCTwnxti2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTwnxti-E_VCTt(3,i))+...
%             E_VCTa(4,i)*(UnemploymentE_VCTwnxti-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTwnxti.*PartyE_VCTwnxti-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTwnxti-E_VCTt(6,i))+...
%             E_VCTa(7,i)*(PresdumE_VCTwnxti-E_VCTt(7,i))+E_VCTa(8,i)*(MidtermE_VCTwnxti-E_VCTt(8,i))));
%         
%         loglikeE_VWnxtnxtI2i(:,1)=(-1/2)*(((1/gammaCT(1,i))*(LOGW_NXT_E_VCTwnxti-E_VCTt(1,i))).^2);
%         
%         loglikeE_VWnxtnxtI3i=(-1)*length(LOGW_NXT_E_VCTwnxti)*log(gammaCT(1,i));
%         
%         loglikeE_VWnxtnxtI4iP(:,1)=(-1)*log((E_VCTa03iP+E_VCTa(2,i)*(XQEVCTwnxti2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTwnxti-E_VCTt(3,i))+...
%             E_VCTa(4,i)*(UnemploymentE_VCTwnxti-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTwnxti.*PartyE_VCTwnxti-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTwnxti-E_VCTt(6,i))+...
%             E_VCTa(7,i)*(PresdumE_VCTwnxti-E_VCTt(7,i))+E_VCTa(8,i)*(MidtermE_VCTwnxti-E_VCTt(8,i))).^2+...
%             (E_VCTa(1,i)^2)*(gammaCT(1,i)^2));
%         
%         loglikeE_VWnxtnxtI4iM(:,1)=(-1)*log((E_VCTa03iM+E_VCTa(2,i)*(XQEVCTwnxti2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTwnxti-E_VCTt(3,i))+...
%             E_VCTa(4,i)*(UnemploymentE_VCTwnxti-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTwnxti.*PartyE_VCTwnxti-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTwnxti-E_VCTt(6,i))+...
%             E_VCTa(7,i)*(PresdumE_VCTwnxti-E_VCTt(7,i))+E_VCTa(8,i)*(MidtermE_VCTwnxti-E_VCTt(8,i))).^2+...
%             (E_VCTa(1,i)^2)*(gammaCT(1,i)^2));
%         
%         loglikeE_VWnxtnxtIi=max(sum(loglikeE_VWnxtnxtI1iP+loglikeE_VWnxtnxtI2i+loglikeE_VWnxtnxtI4iP)+loglikeE_VWnxtnxtI3i,...
%             sum(loglikeE_VWnxtnxtI1iM+loglikeE_VWnxtnxtI2i+loglikeE_VWnxtnxtI4iM)+loglikeE_VWnxtnxtI3i);
%         SRR11=SRR11+loglikeE_VWnxtnxtIi;
%     end
%     SRR11=(-1)*SRR11;
%     
%     if (((E_VCTa(1:8,1).^2)'*gammaCT(1:8,1).^2)>=1)||(((E_VCTa(1:8,2).^2)'*gammaCT(1:8,2).^2)>=1)||(((E_VCTa(1:8,3).^2)'*gammaCT(1:8,3).^2)>=1)
%         for i=1:3
%             SRR11=SRR11+10000000*max(0,((E_VCTa(1:8,i).^2)'*gammaCT(1:8,i).^2-1));
%         end
%     end
%     SRR=SRR11;
%     
% end


if num==31
    
    SRR11=0;
    i=1;
        %       TAU=[];
        TAU=find((TenureE_VCTwnxt([1:183,187:409,411:521])<=Cutoff(2,1))&(TenureE_VCTwnxt([1:183,187:409,411:521])>Cutoff(1,1))); % TAU indexes tenure;
        %       LOGW_NXT_E_VCTwnxti=[];
        LOGW_NXT_E_VCTwnxti=LOGW_NXT_E_VCTwnxt([1:183,187:409,411:521]);
        LOGW_NXT_E_VCTwnxti=LOGW_NXT_E_VCTwnxti(TAU);
        %         LOGLOGD_E_VCTwnxti=LOGLOGD_E_VCTwnxt;
        %         LOGLOGD_E_VCTwnxti(TAU)=[];
        %         LOGLOGTot_E_VCTwnxti=LOGLOGTot_E_VCTwnxt;
        %         LOGLOGTot_E_VCTwnxti(TAU)=[];
        %        XQEVCTwnxti2=[];
        XQEVCTwnxti2=XQEVCTwnxt2([1:183,187:409,411:521]);
        XQEVCTwnxti2=XQEVCTwnxti2(TAU);
        %        LOGW_E_VCTwnxti=[];
        LOGW_E_VCTwnxti=LOGW_E_VCTwnxt([1:183,187:409,411:521]);
        LOGW_E_VCTwnxti=LOGW_E_VCTwnxti(TAU);
        %        PartyE_VCTwnxti=[];
        SameE_VCTwnxti=SameE_VCTwnxt([1:183,187:409,411:521]);
        SameE_VCTwnxti=SameE_VCTwnxti(TAU);
        %        UnemploymentE_VCTwnxti=[];
        UnemploymentE_VCTwnxti=XS_EVCTwnxt_([1:183,187:409,411:521],1);
        UnemploymentE_VCTwnxti=UnemploymentE_VCTwnxti(TAU);
        %        pctWhiteE_VCTwnxti=[];
        PartisanE_VCTwnxti=XS_EVCTwnxt_([1:183,187:409,411:521],2);
        PartisanE_VCTwnxti=PartisanE_VCTwnxti(TAU);
        PartyE_VCTwnxti=PartyE_VCTwnxt([1:183,187:409,411:521]);
        PartyE_VCTwnxti=PartyE_VCTwnxti(TAU);
        
        PresdumE_VCTwnxti=PresdumE_VCTwnxt([1:183,187:409,411:521]);
        PresdumE_VCTwnxti=PresdumE_VCTwnxti(TAU);
        MidtermE_VCTwnxti=MidtermE_VCTwnxt([1:183,187:409,411:521]);
        MidtermE_VCTwnxti=MidtermE_VCTwnxti(TAU);
        loglikeE_VWnxtnxtI1iP=[];
        loglikeE_VWnxtnxtI1iM=[];
        loglikeE_VWnxtnxtI2i=[];
        loglikeE_VWnxtnxtI4iP=[];
        loglikeE_VWnxtnxtI4iM=[];
        
        E_VCTa03iP=(max(0,1-(E_VCTa(1:8,i).^2)'*gammaCT(1:8,i).^2))^(1/2);
        E_VCTa03iM=(-1)*(max(0,1-(E_VCTa(1:8,i).^2)'*gammaCT(1:8,i).^2))^(1/2);
        
        
        loglikeE_VWnxtnxtI1iP(:,1)=2*log(abs(E_VCTa03iP+E_VCTa(1,i)*(LOGW_NXT_E_VCTwnxti-E_VCTt(1,i))+E_VCTa(2,i)*(XQEVCTwnxti2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTwnxti-E_VCTt(3,i))+...
            E_VCTa(4,i)*(UnemploymentE_VCTwnxti-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTwnxti.*PartyE_VCTwnxti-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTwnxti-E_VCTt(6,i))+...
            E_VCTa(7,i)*(PresdumE_VCTwnxti-E_VCTt(7,i))+E_VCTa(8,i)*(MidtermE_VCTwnxti-E_VCTt(8,i))));
        
        loglikeE_VWnxtnxtI1iM(:,1)=2*log(abs(E_VCTa03iM+E_VCTa(1,i)*(LOGW_NXT_E_VCTwnxti-E_VCTt(1,i))+E_VCTa(2,i)*(XQEVCTwnxti2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTwnxti-E_VCTt(3,i))+...
            E_VCTa(4,i)*(UnemploymentE_VCTwnxti-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTwnxti.*PartyE_VCTwnxti-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTwnxti-E_VCTt(6,i))+...
            E_VCTa(7,i)*(PresdumE_VCTwnxti-E_VCTt(7,i))+E_VCTa(8,i)*(MidtermE_VCTwnxti-E_VCTt(8,i))));
        
        loglikeE_VWnxtnxtI2i(:,1)=(-1/2)*(((1/gammaCT(1,i))*(LOGW_NXT_E_VCTwnxti-E_VCTt(1,i))).^2);
        
        loglikeE_VWnxtnxtI3i=(-1)*length(LOGW_NXT_E_VCTwnxti)*log(gammaCT(1,i));
        
        loglikeE_VWnxtnxtI4iP(:,1)=(-1)*log((E_VCTa03iP+E_VCTa(2,i)*(XQEVCTwnxti2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTwnxti-E_VCTt(3,i))+...
            E_VCTa(4,i)*(UnemploymentE_VCTwnxti-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTwnxti.*PartyE_VCTwnxti-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTwnxti-E_VCTt(6,i))+...
            E_VCTa(7,i)*(PresdumE_VCTwnxti-E_VCTt(7,i))+E_VCTa(8,i)*(MidtermE_VCTwnxti-E_VCTt(8,i))).^2+...
            (E_VCTa(1,i)^2)*(gammaCT(1,i)^2));
        
        loglikeE_VWnxtnxtI4iM(:,1)=(-1)*log((E_VCTa03iM+E_VCTa(2,i)*(XQEVCTwnxti2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTwnxti-E_VCTt(3,i))+...
            E_VCTa(4,i)*(UnemploymentE_VCTwnxti-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTwnxti.*PartyE_VCTwnxti-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTwnxti-E_VCTt(6,i))+...
            E_VCTa(7,i)*(PresdumE_VCTwnxti-E_VCTt(7,i))+E_VCTa(8,i)*(MidtermE_VCTwnxti-E_VCTt(8,i))).^2+...
            (E_VCTa(1,i)^2)*(gammaCT(1,i)^2));
        
        loglikeE_VWnxtnxtIi=max(sum(loglikeE_VWnxtnxtI1iP+loglikeE_VWnxtnxtI2i+loglikeE_VWnxtnxtI4iP)+loglikeE_VWnxtnxtI3i,...
            sum(loglikeE_VWnxtnxtI1iM+loglikeE_VWnxtnxtI2i+loglikeE_VWnxtnxtI4iM)+loglikeE_VWnxtnxtI3i);
        SRR11=SRR11+loglikeE_VWnxtnxtIi;
   
    SRR11=(-1)*SRR11;
    
    if (((E_VCTa(1:8,1).^2)'*gammaCT(1:8,1).^2)>=1) 
            SRR11=SRR11+10000000*max(0,((E_VCTa(1:8,i).^2)'*gammaCT(1:8,i).^2-1));
    end
    SRR=SRR11;
    
end

if num==32
    
    SRR11=0;
    i=1;
        %       TAU=[];
        TAU=find((TenureE_VCTwnxt([1:183,187:409,411:521])<=Cutoff(3,1))&(TenureE_VCTwnxt([1:183,187:409,411:521])>Cutoff(2,1))); % TAU indexes tenure;
        %       LOGW_NXT_E_VCTwnxti=[];
        LOGW_NXT_E_VCTwnxti=LOGW_NXT_E_VCTwnxt([1:183,187:409,411:521]);
        LOGW_NXT_E_VCTwnxti=LOGW_NXT_E_VCTwnxti(TAU);
        %         LOGLOGD_E_VCTwnxti=LOGLOGD_E_VCTwnxt;
        %         LOGLOGD_E_VCTwnxti(TAU)=[];
        %         LOGLOGTot_E_VCTwnxti=LOGLOGTot_E_VCTwnxt;
        %         LOGLOGTot_E_VCTwnxti(TAU)=[];
        %        XQEVCTwnxti2=[];
        XQEVCTwnxti2=XQEVCTwnxt2([1:183,187:409,411:521]);
        XQEVCTwnxti2=XQEVCTwnxti2(TAU);
        %        LOGW_E_VCTwnxti=[];
        LOGW_E_VCTwnxti=LOGW_E_VCTwnxt([1:183,187:409,411:521]);
        LOGW_E_VCTwnxti=LOGW_E_VCTwnxti(TAU);
        %        PartyE_VCTwnxti=[];
        SameE_VCTwnxti=SameE_VCTwnxt([1:183,187:409,411:521]);
        SameE_VCTwnxti=SameE_VCTwnxti(TAU);
        %        UnemploymentE_VCTwnxti=[];
        UnemploymentE_VCTwnxti=XS_EVCTwnxt_([1:183,187:409,411:521],1);
        UnemploymentE_VCTwnxti=UnemploymentE_VCTwnxti(TAU);
        %        pctWhiteE_VCTwnxti=[];
        PartisanE_VCTwnxti=XS_EVCTwnxt_([1:183,187:409,411:521],2);
        PartisanE_VCTwnxti=PartisanE_VCTwnxti(TAU);
        PartyE_VCTwnxti=PartyE_VCTwnxt([1:183,187:409,411:521]);
        PartyE_VCTwnxti=PartyE_VCTwnxti(TAU);
        
        PresdumE_VCTwnxti=PresdumE_VCTwnxt([1:183,187:409,411:521]);
        PresdumE_VCTwnxti=PresdumE_VCTwnxti(TAU);
        MidtermE_VCTwnxti=MidtermE_VCTwnxt([1:183,187:409,411:521]);
        MidtermE_VCTwnxti=MidtermE_VCTwnxti(TAU);
        loglikeE_VWnxtnxtI1iP=[];
        loglikeE_VWnxtnxtI1iM=[];
        loglikeE_VWnxtnxtI2i=[];
        loglikeE_VWnxtnxtI4iP=[];
        loglikeE_VWnxtnxtI4iM=[];
        
        E_VCTa03iP=(max(0,1-(E_VCTa(1:8,i).^2)'*gammaCT(1:8,i).^2))^(1/2);
        E_VCTa03iM=(-1)*(max(0,1-(E_VCTa(1:8,i).^2)'*gammaCT(1:8,i).^2))^(1/2);
        
        
        loglikeE_VWnxtnxtI1iP(:,1)=2*log(abs(E_VCTa03iP+E_VCTa(1,i)*(LOGW_NXT_E_VCTwnxti-E_VCTt(1,i))+E_VCTa(2,i)*(XQEVCTwnxti2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTwnxti-E_VCTt(3,i))+...
            E_VCTa(4,i)*(UnemploymentE_VCTwnxti-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTwnxti.*PartyE_VCTwnxti-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTwnxti-E_VCTt(6,i))+...
            E_VCTa(7,i)*(PresdumE_VCTwnxti-E_VCTt(7,i))+E_VCTa(8,i)*(MidtermE_VCTwnxti-E_VCTt(8,i))));
        
        loglikeE_VWnxtnxtI1iM(:,1)=2*log(abs(E_VCTa03iM+E_VCTa(1,i)*(LOGW_NXT_E_VCTwnxti-E_VCTt(1,i))+E_VCTa(2,i)*(XQEVCTwnxti2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTwnxti-E_VCTt(3,i))+...
            E_VCTa(4,i)*(UnemploymentE_VCTwnxti-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTwnxti.*PartyE_VCTwnxti-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTwnxti-E_VCTt(6,i))+...
            E_VCTa(7,i)*(PresdumE_VCTwnxti-E_VCTt(7,i))+E_VCTa(8,i)*(MidtermE_VCTwnxti-E_VCTt(8,i))));
        
        loglikeE_VWnxtnxtI2i(:,1)=(-1/2)*(((1/gammaCT(1,i))*(LOGW_NXT_E_VCTwnxti-E_VCTt(1,i))).^2);
        
        loglikeE_VWnxtnxtI3i=(-1)*length(LOGW_NXT_E_VCTwnxti)*log(gammaCT(1,i));
        
        loglikeE_VWnxtnxtI4iP(:,1)=(-1)*log((E_VCTa03iP+E_VCTa(2,i)*(XQEVCTwnxti2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTwnxti-E_VCTt(3,i))+...
            E_VCTa(4,i)*(UnemploymentE_VCTwnxti-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTwnxti.*PartyE_VCTwnxti-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTwnxti-E_VCTt(6,i))+...
            E_VCTa(7,i)*(PresdumE_VCTwnxti-E_VCTt(7,i))+E_VCTa(8,i)*(MidtermE_VCTwnxti-E_VCTt(8,i))).^2+...
            (E_VCTa(1,i)^2)*(gammaCT(1,i)^2));
        
        loglikeE_VWnxtnxtI4iM(:,1)=(-1)*log((E_VCTa03iM+E_VCTa(2,i)*(XQEVCTwnxti2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTwnxti-E_VCTt(3,i))+...
            E_VCTa(4,i)*(UnemploymentE_VCTwnxti-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTwnxti.*PartyE_VCTwnxti-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTwnxti-E_VCTt(6,i))+...
            E_VCTa(7,i)*(PresdumE_VCTwnxti-E_VCTt(7,i))+E_VCTa(8,i)*(MidtermE_VCTwnxti-E_VCTt(8,i))).^2+...
            (E_VCTa(1,i)^2)*(gammaCT(1,i)^2));
        
        loglikeE_VWnxtnxtIi=max(sum(loglikeE_VWnxtnxtI1iP+loglikeE_VWnxtnxtI2i+loglikeE_VWnxtnxtI4iP)+loglikeE_VWnxtnxtI3i,...
            sum(loglikeE_VWnxtnxtI1iM+loglikeE_VWnxtnxtI2i+loglikeE_VWnxtnxtI4iM)+loglikeE_VWnxtnxtI3i);
        SRR11=SRR11+loglikeE_VWnxtnxtIi;
   
    SRR11=(-1)*SRR11;
    
    if (((E_VCTa(1:8,1).^2)'*gammaCT(1:8,1).^2)>=1)
            SRR11=SRR11+10000000*max(0,((E_VCTa(1:8,i).^2)'*gammaCT(1:8,i).^2-1));
    end
    SRR=SRR11;
    
end

if num==33
    
    SRR11=0;
    i=1;
        %       TAU=[];
        TAU=find((TenureE_VCTwnxt([1:183,187:409,411:521])<=Cutoff(4,1))&(TenureE_VCTwnxt([1:183,187:409,411:521])>Cutoff(3,1))); % TAU indexes tenure;
        %       LOGW_NXT_E_VCTwnxti=[];
        LOGW_NXT_E_VCTwnxti=LOGW_NXT_E_VCTwnxt([1:183,187:409,411:521]);
        LOGW_NXT_E_VCTwnxti=LOGW_NXT_E_VCTwnxti(TAU);
        %         LOGLOGD_E_VCTwnxti=LOGLOGD_E_VCTwnxt;
        %         LOGLOGD_E_VCTwnxti(TAU)=[];
        %         LOGLOGTot_E_VCTwnxti=LOGLOGTot_E_VCTwnxt;
        %         LOGLOGTot_E_VCTwnxti(TAU)=[];
        %        XQEVCTwnxti2=[];
        XQEVCTwnxti2=XQEVCTwnxt2([1:183,187:409,411:521]);
        XQEVCTwnxti2=XQEVCTwnxti2(TAU);
        %        LOGW_E_VCTwnxti=[];
        LOGW_E_VCTwnxti=LOGW_E_VCTwnxt([1:183,187:409,411:521]);
        LOGW_E_VCTwnxti=LOGW_E_VCTwnxti(TAU);
        %        PartyE_VCTwnxti=[];
        SameE_VCTwnxti=SameE_VCTwnxt([1:183,187:409,411:521]);
        SameE_VCTwnxti=SameE_VCTwnxti(TAU);
        %        UnemploymentE_VCTwnxti=[];
        UnemploymentE_VCTwnxti=XS_EVCTwnxt_([1:183,187:409,411:521],1);
        UnemploymentE_VCTwnxti=UnemploymentE_VCTwnxti(TAU);
        %        pctWhiteE_VCTwnxti=[];
        PartisanE_VCTwnxti=XS_EVCTwnxt_([1:183,187:409,411:521],2);
        PartisanE_VCTwnxti=PartisanE_VCTwnxti(TAU);
        PartyE_VCTwnxti=PartyE_VCTwnxt([1:183,187:409,411:521]);
        PartyE_VCTwnxti=PartyE_VCTwnxti(TAU);
        
        PresdumE_VCTwnxti=PresdumE_VCTwnxt([1:183,187:409,411:521]);
        PresdumE_VCTwnxti=PresdumE_VCTwnxti(TAU);
        MidtermE_VCTwnxti=MidtermE_VCTwnxt([1:183,187:409,411:521]);
        MidtermE_VCTwnxti=MidtermE_VCTwnxti(TAU);
        loglikeE_VWnxtnxtI1iP=[];
        loglikeE_VWnxtnxtI1iM=[];
        loglikeE_VWnxtnxtI2i=[];
        loglikeE_VWnxtnxtI4iP=[];
        loglikeE_VWnxtnxtI4iM=[];
        
        E_VCTa03iP=(max(0,1-(E_VCTa(1:8,i).^2)'*gammaCT(1:8,i).^2))^(1/2);
        E_VCTa03iM=(-1)*(max(0,1-(E_VCTa(1:8,i).^2)'*gammaCT(1:8,i).^2))^(1/2);
        
        
        loglikeE_VWnxtnxtI1iP(:,1)=2*log(abs(E_VCTa03iP+E_VCTa(1,i)*(LOGW_NXT_E_VCTwnxti-E_VCTt(1,i))+E_VCTa(2,i)*(XQEVCTwnxti2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTwnxti-E_VCTt(3,i))+...
            E_VCTa(4,i)*(UnemploymentE_VCTwnxti-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTwnxti.*PartyE_VCTwnxti-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTwnxti-E_VCTt(6,i))+...
            E_VCTa(7,i)*(PresdumE_VCTwnxti-E_VCTt(7,i))+E_VCTa(8,i)*(MidtermE_VCTwnxti-E_VCTt(8,i))));
        
        loglikeE_VWnxtnxtI1iM(:,1)=2*log(abs(E_VCTa03iM+E_VCTa(1,i)*(LOGW_NXT_E_VCTwnxti-E_VCTt(1,i))+E_VCTa(2,i)*(XQEVCTwnxti2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTwnxti-E_VCTt(3,i))+...
            E_VCTa(4,i)*(UnemploymentE_VCTwnxti-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTwnxti.*PartyE_VCTwnxti-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTwnxti-E_VCTt(6,i))+...
            E_VCTa(7,i)*(PresdumE_VCTwnxti-E_VCTt(7,i))+E_VCTa(8,i)*(MidtermE_VCTwnxti-E_VCTt(8,i))));
        
        loglikeE_VWnxtnxtI2i(:,1)=(-1/2)*(((1/gammaCT(1,i))*(LOGW_NXT_E_VCTwnxti-E_VCTt(1,i))).^2);
        
        loglikeE_VWnxtnxtI3i=(-1)*length(LOGW_NXT_E_VCTwnxti)*log(gammaCT(1,i));
        
        loglikeE_VWnxtnxtI4iP(:,1)=(-1)*log((E_VCTa03iP+E_VCTa(2,i)*(XQEVCTwnxti2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTwnxti-E_VCTt(3,i))+...
            E_VCTa(4,i)*(UnemploymentE_VCTwnxti-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTwnxti.*PartyE_VCTwnxti-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTwnxti-E_VCTt(6,i))+...
            E_VCTa(7,i)*(PresdumE_VCTwnxti-E_VCTt(7,i))+E_VCTa(8,i)*(MidtermE_VCTwnxti-E_VCTt(8,i))).^2+...
            (E_VCTa(1,i)^2)*(gammaCT(1,i)^2));
        
        loglikeE_VWnxtnxtI4iM(:,1)=(-1)*log((E_VCTa03iM+E_VCTa(2,i)*(XQEVCTwnxti2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTwnxti-E_VCTt(3,i))+...
            E_VCTa(4,i)*(UnemploymentE_VCTwnxti-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTwnxti.*PartyE_VCTwnxti-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTwnxti-E_VCTt(6,i))+...
            E_VCTa(7,i)*(PresdumE_VCTwnxti-E_VCTt(7,i))+E_VCTa(8,i)*(MidtermE_VCTwnxti-E_VCTt(8,i))).^2+...
            (E_VCTa(1,i)^2)*(gammaCT(1,i)^2));
        
        loglikeE_VWnxtnxtIi=max(sum(loglikeE_VWnxtnxtI1iP+loglikeE_VWnxtnxtI2i+loglikeE_VWnxtnxtI4iP)+loglikeE_VWnxtnxtI3i,...
            sum(loglikeE_VWnxtnxtI1iM+loglikeE_VWnxtnxtI2i+loglikeE_VWnxtnxtI4iM)+loglikeE_VWnxtnxtI3i);
        SRR11=SRR11+loglikeE_VWnxtnxtIi;
   
    SRR11=(-1)*SRR11;
    
    if (((E_VCTa(1:8,1).^2)'*gammaCT(1:8,1).^2)>=1)
            SRR11=SRR11+10000000*max(0,((E_VCTa(1:8,i).^2)'*gammaCT(1:8,i).^2-1));
    end
    SRR=SRR11;
    
end
end



