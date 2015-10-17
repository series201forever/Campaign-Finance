function SRR=Minimize_Apr2013(thetain,num)

% global dF_gamma_ct
% global dF_total_ct
% global dF_nxt_nxt_ct
% global dF_gamma_nct
% global dF_total_nct
% global dF_nxt_nxt_nct
% global Entry
% global Winrnd
% global Ret
% global Betawh
% global Betaump
% global epswh
% global epsump
% global Shockump
% global Shockwh
%global LOGD_NC
%global LOGW_NC
%global LOGWNXT_NC
%global White
%global White_NC
%global Unemployment
%global Unemployment_NC
%global Party
%global X_Knot2

global iterate
global Delt

%Variables used in num=1
global Samplesize
global LOGW_I
global Contest
global RTotD_NC
global Tenure
global Primary_N
global Unempsame
global Partdemo
global X_Knot1



global SamplesizeAC
global LOGD_IAC
global LOGW_IAC
global VSAC
global LOGTot_NCAC
global RTotD_NCAC
global NCAC
global LOGD_CAC
global PartyAC
global TenureAC
global Tenure_NCAC
global WhiteAC
global White_NCAC
global Win_AC
global UnemploymentAC
global Unemployment_NCAC
global LOGLOGTot_NCAC
global LOGLOGD_NCAC
global X_KnotAC1
global X_KnotAC2
global E_VContestFUL
global E_VNContestFUL
global NCE_V
global NCE_VCT
global NCE_VNCT
global NCE_VCTwnxt
global NCE_VTenCT
global NCE_VTenNCT
global NCE_VTenCTwnxt
global PartyE_V
global PartyE_VCT
global PartyE_VNCT
global PartyE_VCTwnxt
global XSNCEVCT_
global XSNCEVNCT_
global TenureE_VCTwnxt
global XS_EVCT_
global XS_EVNCT_
global XS_EVCTwnxt_
global XSNCEVCTwnxt_
global XSEV_
global TenureE_VCT
global TenureE_VNCT
global TenureE_V
global TenureE_VCTwnxt
global LOGTotal_E_VNCT
global LOGD_E_VCT
global LOGD_E_VNCT
global LOGW_NXT_E_VNCT
global LOGW_NXT_E_V
global LOGW_E_VCT
global LOGW_E_VCTwnxt
global LOGW_E_VNCT
global LOGW_NXT_E_VCT
global LOGW_NXT_E_VCTwnxt
global LOGTotal_E_VCT
global LOGLOGD_E_VCT
global LOGLOGD_E_VNCT
global LOGLOGTot_E_VCT
global LOGLOGTot_E_VNCT
global LOGLOGD_E_VCTwnxt
global LOGLOGTot_E_VCTwnxt
% global RTotDE_VCT
% global RTotDE_VNCT
% global RTotDE_V
% global RTotDE_VCTwnxt
global X_KnotEV1
global X_KnotEV2
global Sofarbest
global bestiter
global LOGTot_E_VC
global YearAC
global YEARE_VCT
global YEARE_VNCT
global YEARE_VCTwnxt
global LOGW_CNXTAC
global EPS
global EPS2
global EPS3
global NumSim
global Year
global XXX
global Tenure_NC
global mesh
global IND5
global Est1
global Est2
%thetain=[thetaS;thetaQ;theta;B_I;B_C;Q_C1;Q_C2;B_T;E_VCTa;E_VCTt;gammaCT...
%    E_VNCTa;E_VNCTt;gammaNCT;thetawin;cost1;ben1;sig;cost2];

if num==1
    thetaS=thetain(1:2,1);
    thetaS(1,1)=1;
    thetaQ(1,1)=1;
    % for ii=2:size(X_Knot1,2)
    %     thetaQ(ii,1)=1-sum(abs(thetain(3:(1+ii),1)));
    % end
    thetaQ(2:9,1)=thetain(3:10,1);
    %B_T=thetain(19,1);
    theta=thetain(11:15,1);
    theta2=thetain(16:20,1);
    % B_I=thetain(30,1);
    % B_C=thetain(31,1);
    % Q_C1=thetain(32,1);
    % Q_C2=thetain(33,1);
    % Q_C3=thetain(34,1);
elseif num==2
    thetaS=Est1(1:2,1);
    thetaS(1,1)=1;
    thetaQ(1,1)=1;
    %     for ii=2:size(X_Knot1,2)
    %         thetaQ(ii,1)=1-sum(abs(Est1(3:(ii+1),1)));
    %     end
    thetaQ(2:9,1)=Est1(3:10,1);
    thetaS2=thetain(1:2,1); % do not normalize.
    thetaQ2(1,1)=thetain(3,1);
    %     for ii=2:size(X_Knot1,2)
    %         thetaQ2(ii,1)=thetaQ2(1,1)-sum(abs(thetain(4:(ii+2),1)));
    %     end
    thetaQ2(2:9,1)=thetain(4:11,1);
    B_T=thetain(12,1);
    theta=Est1(11:15,1);
    theta2=Est1(16:20,1);
    B_I=thetain(13,1);
    B_C=thetain(14,1);
    Q_C1=thetain(15,1);
    Q_C2=thetain(16,1);
    Q_C3=thetain(17,1);
elseif num==3
    thetawin=thetain(1:7,1);
    thetaS2=Est2(1:2,1); % do not normalize.
    thetaQ2(1,1)=Est2(3,1);
    %     for ii=2:size(X_Knot1,2)
    %         thetaQ2(ii,1)=thetaQ2(1,1)-sum(abs(thetain(4:(ii+2),1)));
    %     end
    thetaQ2(2:9,1)=Est2(4:11,1);
    
elseif num==4
    thetaS2=Est2(1:2,1); % do not normalize.
    thetaQ2(1,1)=Est2(3,1);
    %     for ii=2:size(X_Knot1,2)
    %         thetaQ2(ii,1)=thetaQ2(1,1)-sum(abs(thetain(4:(ii+2),1)));
    %     end
    thetaQ2(2:9,1)=Est2(4:11,1);
    menc=-5;
    E_VCTa(:,1)=thetain(menc+6:menc+20,1);
    E_VCTa(:,2)=thetain(menc+21:menc+35,1);
    E_VCTa(:,3)=thetain(menc+36:menc+50,1);
    E_VCTt(:,1)=thetain(menc+51:menc+65,1);
    E_VCTt(:,2)=thetain(menc+66:menc+80,1);
    E_VCTt(:,3)=thetain(menc+81:menc+95,1);
    gammaCT(:,1)=thetain(menc+96:menc+110,1);
    gammaCT(:,2)=thetain(menc+111:menc+125,1);
    gammaCT(:,3)=thetain(menc+126:menc+140,1);
elseif num==5
    thetaS2=Est2(1:2,1); % do not normalize.
    thetaQ2(1,1)=Est2(3,1);
    %     for ii=2:size(X_Knot1,2)
    %         thetaQ2(ii,1)=thetaQ2(1,1)-sum(abs(thetain(4:(ii+2),1)));
    %     end
    thetaQ2(2:9,1)=Est2(4:11,1);
    gammaNCT(:,1)=thetain(1:11,1);
    tNCT(:,1)=thetain(12:22,1);
    wNCT(:,1)=thetain(23:33,1);
end



% beta_1=thetain(333,1);
% beta_2=thetain(334,1);
% sig=thetain(335,1);






if num<=2
    XS=[Unempsame,Partdemo]*thetaS;
    % XSNC=[Unemployment_NC,White_NC]*thetaS;
    % XSNC2=[Unemployment_NC,White_NC]*thetaS2;
    XQ=X_Knot1*thetaQ;  % q_I
    
    XSAC=[UnemploymentAC,WhiteAC]*thetaS;   %State in this period.
    % XSNCAC=[Unemployment_NCAC,White_NCAC]*thetaS;   %State in uncontested periods.
    % XSNCAC2=[Unemployment_NCAC,White_NCAC]*thetaS2;   %State in uncontested periods.
    XQAC=X_KnotAC1*thetaQ;  % q_I
    
    PRENTAC=normcdf([ones(SamplesizeAC,1),LOGW_IAC,XQAC,XSAC.*PartyAC,log(TenureAC+1)]*theta,0,1); % Entry prob.
    E_Primary_NAC=[ones(SamplesizeAC,1),LOGW_IAC,XQAC,XSAC.*PartyAC,log(TenureAC+1)]*theta2; %Number of expected entrants.
    % XSNCEVCT=XSNCEVCT_*thetaS;
    % XSNCEVNCT=XSNCEVNCT_*thetaS;
    %
    % XSNCEVCTwnxt=XSNCEVCTwnxt_*thetaS;
    
    XS_EVCT=XS_EVCT_*thetaS;
    
    XS_EVNCT=XS_EVNCT_*thetaS;
    
    % XSEV=XSEV_*thetaS;
    % XS_EVCTwnxt=XS_EVCTwnxt_*thetaS;
    % lenE_VCT=length(XS_EVCT);
    % lenE_VNCT=length(XS_EVNCT);
    
    XQEV=X_KnotEV1*thetaQ;  % q_I
    XQEVCT=XQEV;
    XQEVCT(E_VContestFUL,:)=[];
end
if num>=2
    XS2=[Unemployment,White]*thetaS2;
    XQ2=X_Knot1*thetaQ2;  % q_I
    XSAC2=[UnemploymentAC,WhiteAC]*thetaS2;   %State in this period.
    XQAC2=X_KnotAC1*thetaQ2;  % q_I
    XS_EVCT2=XS_EVCT_*thetaS2;
    XS_EVNCT2=XS_EVNCT_*thetaS2;
    XQEV2=X_KnotEV1*thetaQ2;  % q_I
    XQEVCT2=XQEV2;
    XQEVCT2(E_VContestFUL,:)=[];
    XQEVCTwnxt2=XQEVCT2;
    XQEVCTwnxt2(IND5,:)=[];
    XQEVNCT2=XQEV2;
    XQEVNCT2(E_VNContestFUL,:)=[];
end

% XQEVNCT=XQEV;
% XQEVNCT(E_VNContestFUL,:)=[];


% XXX=XQEVCT;
%
%
% XXX=XQEVCT;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Estimation of Probability     %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SOLVE
%% min(theta) sum({contest-PHI(X.*theta)}.^2)
if num==1
    dL=(Contest-[ones(Samplesize,1),LOGW_I,XQ,XS,log(Tenure+1)]*theta).^2;
    SRR1=sum(dL,1);
    SRR1=SRR1/var(Contest);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Estimation of E[primary_N]     %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SOLVE
%% min(theta) sum({primary_N-X.*theta)}.^2)
if num==1
    dL=(Primary_N-[ones(Samplesize,1),LOGW_I,XQ,XS.*Party,log(Tenure+1)]*theta2).^2;
    SRR3=sum(dL,1);
    SRR3=SRR3/var(Primary_N);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%      Ai & Chen     %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% [rho(z_1,alpha),rho(z_2,alpha),...,rho(z_n,alpha)]'
if num==2
    RHO=VSAC-0.5-B_I*LOGD_IAC-B_C*LOGD_CAC-(XSAC2-mean(XSAC2)).*PartyAC-XQAC2-PRENTAC*Q_C1-E_Primary_NAC*Q_C2-PRENTAC.*E_Primary_NAC*Q_C3-B_T*log(TenureAC+1);   %% sign of -[PRENTAC,PRENTAC.^2]*[Q_C1,Q_C2]' unlike in paper.
    P=[ones(length(VSAC),1),log(TenureAC+1),UnemploymentAC.*PartyAC,WhiteAC.*PartyAC,X_KnotAC1(:,1:8),LOGW_IAC,(log(TenureAC+1)).^2,X_KnotAC1(:,1:8).^2,LOGW_IAC.^2,...
        log(TenureAC+1).*PartyAC, log(TenureAC+1).*LOGW_IAC, PartyAC.*LOGW_IAC]; %% Basis function
    
    SRR2=RHO'*P*inv(P'*P)*P'*RHO;  %% EQ (6) in Ai and Chen p1799
    SRR2=SRR2/var(VSAC);
end

if num==4
    Cutoff=[-1;3;7;100];
    SRR9=0;
    for i=1:3
        TAU=[];
        TAU=find((TenureE_VCT<=Cutoff(i+1,1))&(TenureE_VCT>Cutoff(i,1))); % TAU indexes tenure;
        LOGD_E_VCTi=[];
        LOGD_E_VCTi=LOGD_E_VCT(TAU);
%         LOGLOGD_E_VCTi=LOGLOGD_E_VCT;
%         LOGLOGD_E_VCTi(TAU)=[];
%         LOGLOGTot_E_VCTi=LOGLOGTot_E_VCT;
%         LOGLOGTot_E_VCTi(TAU)=[];
        XQEVCTi2=[];
        XQEVCTi2=XQEVCT2(TAU);
        
        LOGW_E_VCTi=[];
        LOGW_E_VCTi=LOGW_E_VCT(TAU);
        PartyE_VCTi=[];
        PartyE_VCTi=PartyE_VCT(TAU);
        UnemploymentE_VCTi=[];
        UnemploymentE_VCTi=XS_EVCT_(TAU,1);
        pctWhiteE_VCTi=[];
        pctWhiteE_VCTi=XS_EVCT_(TAU,2);
        loglikeE_VDI1iP=[];
        loglikeE_VDI1iM=[];
        loglikeE_VDI2i=[];
        loglikeE_VDI4iP=[];
        loglikeE_VDI4iM=[];
        
        
        E_VCTa01iP=(1-(E_VCTa(1:5,i).^2)'*gammaCT(1:5,i).^2)^(1/2);
        E_VCTa01iM=(-1)*(1-(E_VCTa(1:5,i).^2)'*gammaCT(1:5,i).^2)^(1/2);
        loglikeE_VDI1iP(:,1)=2*log(abs(E_VCTa01iP+E_VCTa(1,i)*(LOGD_E_VCTi-E_VCTt(1,i))+E_VCTa(2,i)*(XQEVCTi2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTi-E_VCTt(3,i))+...
            E_VCTa(4,i)*(PartyE_VCTi.*UnemploymentE_VCTi-E_VCTt(4,i))+E_VCTa(5,i)*(PartyE_VCTi.*pctWhiteE_VCTi-E_VCTt(5,i))));
        loglikeE_VDI1iM(:,1)=2*log(abs(E_VCTa01iM+E_VCTa(1,i)*(LOGD_E_VCTi-E_VCTt(1,i))+E_VCTa(2,i)*(XQEVCTi2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTi-E_VCTt(3,i))+...
            E_VCTa(4,i)*(PartyE_VCTi.*UnemploymentE_VCTi-E_VCTt(4,i))+E_VCTa(5,i)*(PartyE_VCTi.*pctWhiteE_VCTi-E_VCTt(5,i))));
        
        loglikeE_VDI2i(:,1)=(-1/2)*(((1/gammaCT(1,i))*(LOGD_E_VCTi-E_VCTt(1,i))).^2);
        
        loglikeE_VDI3i=(-1)*length(LOGD_E_VCTi)*log(gammaCT(1,i));
        
        loglikeE_VDI4iP(:,1)=(-1)*log((E_VCTa01iP+E_VCTa(2,i)*(XQEVCTi2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTi-E_VCTt(3,i))+...
            E_VCTa(4,i)*(PartyE_VCTi.*UnemploymentE_VCTi-E_VCTt(4,i))+E_VCTa(5,i)*(PartyE_VCTi.*pctWhiteE_VCTi-E_VCTt(5,i))).^2+...
        (E_VCTa(1,i)^2)*(gammaCT(1,i)^2));
        
        loglikeE_VDI4iM(:,1)=(-1)*log((E_VCTa01iM+E_VCTa(2,i)*(XQEVCTi2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTi-E_VCTt(3,i))+...
            E_VCTa(4,i)*(PartyE_VCTi.*UnemploymentE_VCTi-E_VCTt(4,i))+E_VCTa(5,i)*(PartyE_VCTi.*pctWhiteE_VCTi-E_VCTt(5,i))).^2+...
            (E_VCTa(1,i)^2)*(gammaCT(1,i)^2));
        
        loglikeE_VDIi=max(sum(loglikeE_VDI1iP+loglikeE_VDI2i+loglikeE_VDI4iP)+loglikeE_VDI3i,...
            sum(loglikeE_VDI1iM+loglikeE_VDI2i+loglikeE_VDI4iM)+loglikeE_VDI3i);
        SRR9=SRR9+loglikeE_VDIi;
    end
    SRR9=(-1)*SRR9;
    
    
    
    %%%%%%%%   F(total_I|w',q_I,s,Entry,contest=1)  %%%%%%%%%
    %% Here, we want the conditional density of
    %% F(total_I|w',q_I,s,tenure,Entry). Use Hermite expansion.
    %% denom is the denominator that is common to most of the derivative of the
    %% Loglikelihood(of total_I) function.
    
    SRR10=0;
    for i=1:3
        TAU=[];
        TAU=find((TenureE_VCT<=Cutoff(i+1,1))&(TenureE_VCT>Cutoff(i,1))); % TAU indexes tenure;
        LOGTotal_E_VCTi=[];
        LOGTotal_E_VCTi=LOGTotal_E_VCT(TAU);
%         LOGLOGD_E_VCTi=LOGLOGD_E_VCT;
%         LOGLOGD_E_VCTi(TAU)=[];
%         LOGLOGTot_E_VCTi=LOGLOGTot_E_VCT;
%         LOGLOGTot_E_VCTi(TAU)=[];
        XQEVCTi2=[];
        XQEVCTi2=XQEVCT2(TAU);
        LOGW_E_VCTi=[];
        LOGW_E_VCTi=LOGW_E_VCT(TAU);
        PartyE_VCTi=[];
        PartyE_VCTi=PartyE_VCT(TAU);
        UnemploymentE_VCTi=[];
        UnemploymentE_VCTi=XS_EVCT_(TAU,1);
        pctWhiteE_VCTi=[];
        pctWhiteE_VCTi=XS_EVCT_(TAU,2);
        loglikeE_VTotI1iP=[];
        loglikeE_VTotI1iM=[];
        loglikeE_VTotI2i=[];
        loglikeE_VTotI4iP=[];
        loglikeE_VTotI4iM=[];
        
        E_VCTa02iP=(1-(E_VCTa(6:10,i).^2)'*gammaCT(6:10,i).^2)^(1/2);
        E_VCTa02iM=(-1)*(1-(E_VCTa(6:10,i).^2)'*gammaCT(6:10,i).^2)^(1/2);
        
        loglikeE_VTotI1iP(:,1)=2*log(abs(E_VCTa02iP+E_VCTa(6,i)*(LOGTotal_E_VCTi-E_VCTt(6,i))+E_VCTa(7,i)*(XQEVCTi2-E_VCTt(7,i))+E_VCTa(8,i)*(LOGW_E_VCTi-E_VCTt(8,i))+...
            E_VCTa(9,i)*(PartyE_VCTi.*UnemploymentE_VCTi-E_VCTt(9,i))+E_VCTa(10,i)*(PartyE_VCTi.*pctWhiteE_VCTi-E_VCTt(10,i))));
        
        loglikeE_VTotI1iM(:,1)=2*log(abs(E_VCTa02iM+E_VCTa(6,i)*(LOGTotal_E_VCTi-E_VCTt(6,i))+E_VCTa(7,i)*(XQEVCTi2-E_VCTt(7,i))+E_VCTa(8,i)*(LOGW_E_VCTi-E_VCTt(8,i))+...
            E_VCTa(9,i)*(PartyE_VCTi.*UnemploymentE_VCTi-E_VCTt(9,i))+E_VCTa(10,i)*(PartyE_VCTi.*pctWhiteE_VCTi-E_VCTt(10,i))));
        
        loglikeE_VTotI2i(:,1)=(-1/2)*(((1/gammaCT(6,i))*(LOGTotal_E_VCTi-E_VCTt(6,i))).^2);
        
        loglikeE_VTotI3i=(-1)*length(LOGTotal_E_VCTi)*log(gammaCT(6,i));
        
        loglikeE_VTotI4iP(:,1)=(-1)*log((E_VCTa02iP+E_VCTa(7,i)*(XQEVCTi2-E_VCTt(7,i))+E_VCTa(8,i)*(LOGW_E_VCTi-E_VCTt(8,i))+...
            E_VCTa(9,i)*(PartyE_VCTi.*UnemploymentE_VCTi-E_VCTt(9,i))+E_VCTa(10,i)*(PartyE_VCTi.*pctWhiteE_VCTi-E_VCTt(10,i))).^2+...
        (E_VCTa(6,i)^2)*(gammaCT(6,i)^2));
        
        loglikeE_VTotI4iM(:,1)=(-1)*log((E_VCTa02iM+E_VCTa(7,i)*(XQEVCTi2-E_VCTt(7,i))+E_VCTa(8,i)*(LOGW_E_VCTi-E_VCTt(8,i))+...
            E_VCTa(9,i)*(PartyE_VCTi.*UnemploymentE_VCTi-E_VCTt(9,i))+E_VCTa(10,i)*(PartyE_VCTi.*pctWhiteE_VCTi-E_VCTt(10,i))).^2+...
        (E_VCTa(6,i)^2)*(gammaCT(6,i)^2));
        
        loglikeE_VTotIi=max(sum(loglikeE_VTotI1iP+loglikeE_VTotI2i+loglikeE_VTotI4iP)+loglikeE_VTotI3i,...
            sum(loglikeE_VTotI1iM+loglikeE_VTotI2i+loglikeE_VTotI4iM)+loglikeE_VTotI3i);
        SRR10=SRR10+loglikeE_VTotIi;
    end
    SRR10=(-1)*SRR10;
    
    
    SRR11=0;
    for i=1:3
        TAU=[];
        TAU=find((TenureE_VCTwnxt<=Cutoff(i+1,1))&(TenureE_VCTwnxt>Cutoff(i,1))); % TAU indexes tenure;
        LOGW_NXT_E_VCTwnxti=[];
        LOGW_NXT_E_VCTwnxti=LOGW_NXT_E_VCTwnxt(TAU);
%         LOGLOGD_E_VCTwnxti=LOGLOGD_E_VCTwnxt;
%         LOGLOGD_E_VCTwnxti(TAU)=[];
%         LOGLOGTot_E_VCTwnxti=LOGLOGTot_E_VCTwnxt;
%         LOGLOGTot_E_VCTwnxti(TAU)=[];
        XQEVCTwnxti2=[];
        XQEVCTwnxti2=XQEVCTwnxt2(TAU);
        LOGW_E_VCTwnxti=[];
        LOGW_E_VCTwnxti=LOGW_E_VCTwnxt(TAU);
        PartyE_VCTwnxti=[];
        PartyE_VCTwnxti=PartyE_VCTwnxt(TAU);
        UnemploymentE_VCTwnxti=[];
        UnemploymentE_VCTwnxti=XS_EVCTwnxt_(TAU,1);
        pctWhiteE_VCTwnxti=[];
        pctWhiteE_VCTwnxti=XS_EVCTwnxt_(TAU,2);
        E_VCTa03iP=(1-(E_VCTa(11:15,i).^2)'*gammaCT(11:15,i).^2)^(1/2);
        E_VCTa03iM=(-1)*(1-(E_VCTa(11:15,i).^2)'*gammaCT(11:15,i).^2)^(1/2);
        loglikeE_VWnxtnxtI1iP=[];
        loglikeE_VWnxtnxtI1iM=[];
        loglikeE_VWnxtnxtI2i=[];
        loglikeE_VWnxtnxtI4iP=[];
        loglikeE_VWnxtnxtI4iM=[];
        
        loglikeE_VWnxtnxtI1iP(:,1)=2*log(abs(E_VCTa03iP+E_VCTa(11,i)*(LOGW_NXT_E_VCTwnxti-E_VCTt(11,i))+E_VCTa(12,i)*(XQEVCTwnxti2-E_VCTt(12,i))+E_VCTa(13,i)*(LOGW_E_VCTwnxti-E_VCTt(13,i))+...
            E_VCTa(14,i)*(PartyE_VCTwnxti.*UnemploymentE_VCTwnxti-E_VCTt(14,i))+E_VCTa(15,i)*(PartyE_VCTwnxti.*pctWhiteE_VCTwnxti-E_VCTt(15,i))));
        
        loglikeE_VWnxtnxtI1iM(:,1)=2*log(abs(E_VCTa03iM+E_VCTa(11,i)*(LOGW_NXT_E_VCTwnxti-E_VCTt(11,i))+E_VCTa(12,i)*(XQEVCTwnxti2-E_VCTt(12,i))+E_VCTa(13,i)*(LOGW_E_VCTwnxti-E_VCTt(13,i))+...
            E_VCTa(14,i)*(PartyE_VCTwnxti.*UnemploymentE_VCTwnxti-E_VCTt(14,i))+E_VCTa(15,i)*(PartyE_VCTwnxti.*pctWhiteE_VCTwnxti-E_VCTt(15,i))));
        
        loglikeE_VWnxtnxtI2i(:,1)=(-1/2)*(((1/gammaCT(11,i))*(LOGW_NXT_E_VCTwnxti-E_VCTt(11,i))).^2);
        
        loglikeE_VWnxtnxtI3i=(-1)*length(LOGW_NXT_E_VCTwnxti)*log(gammaCT(11,i));
        
        loglikeE_VWnxtnxtI4iP(:,1)=(-1)*log((E_VCTa03iP+E_VCTa(12,i)*(XQEVCTwnxti2-E_VCTt(12,i))+E_VCTa(13,i)*(LOGW_E_VCTwnxti-E_VCTt(13,i))+...
            E_VCTa(14,i)*(PartyE_VCTwnxti.*UnemploymentE_VCTwnxti-E_VCTt(14,i))+E_VCTa(15,i)*(PartyE_VCTwnxti.*pctWhiteE_VCTwnxti-E_VCTt(15,i))).^2+...
        (E_VCTa(11,i)^2)*(gammaCT(11,i)^2));
        
        loglikeE_VWnxtnxtI4iM(:,1)=(-1)*log((E_VCTa03iP+E_VCTa(12,i)*(XQEVCTwnxti2-E_VCTt(12,i))+E_VCTa(13,i)*(LOGW_E_VCTwnxti-E_VCTt(13,i))+...
            E_VCTa(14,i)*(PartyE_VCTwnxti.*UnemploymentE_VCTwnxti-E_VCTt(14,i))+E_VCTa(15,i)*(PartyE_VCTwnxti.*pctWhiteE_VCTwnxti-E_VCTt(15,i))).^2+...
            (E_VCTa(11,i)^2)*(gammaCT(11,i)^2));
        
        loglikeE_VWnxtnxtIi=max(sum(loglikeE_VWnxtnxtI1iP+loglikeE_VWnxtnxtI2i+loglikeE_VWnxtnxtI4iP)+loglikeE_VWnxtnxtI3i,...
            sum(loglikeE_VWnxtnxtI1iM+loglikeE_VWnxtnxtI2i+loglikeE_VWnxtnxtI4iM)+loglikeE_VWnxtnxtI3i);
        SRR11=SRR11+loglikeE_VWnxtnxtIi;
    end
    SRR11=(-1)*SRR11;
    
end
if num==5
    
    %%%%%%%%   gamma_I=gamma_I(w',q_I,s,Entry,tenure|contest=0)  %%%%%%%%%
    %% Here, we want the (deterministic) function of spending|contest=0 as a
    %% function of w',q_I,s,Entry,tenure
    Xgamma_INCT=[ones(size(LOGD_E_VNCT,1),1),LOGW_E_VNCT,PartyE_VNCT.*XS_EVNCT_(:,1),PartyE_VNCT.*XS_EVNCT_(:,2),TenureE_VNCT,XQEVNCT2,...
        LOGW_E_VNCT.^2,PartyE_VNCT.*(XS_EVNCT_(:,1)).^2,PartyE_VNCT.*(XS_EVNCT_(:,2)).^2,TenureE_VNCT.^2,XQEVNCT2.^2];
    SRR12=sum((LOGD_E_VNCT-Xgamma_INCT*gammaNCT).^2);
    
    
    %%%%%%%   total_I=tot_I(w',q_I,s,Entry|contest=0)  %%%%%%%%%
    % Here, we want the deterministic function of fund-raising conditional on
    % being uncontested as a function of w',q_I,s,tenure,Entry
    
    Xt_INCT=[ones(size(LOGD_E_VNCT,1),1),LOGW_E_VNCT,PartyE_VNCT.*XS_EVNCT_(:,1),PartyE_VNCT.*XS_EVNCT_(:,2),TenureE_VNCT,XQEVNCT2,...
        LOGW_E_VNCT.^2,PartyE_VNCT.*(XS_EVNCT_(:,1)).^2,PartyE_VNCT.*(XS_EVNCT_(:,2)).^2,TenureE_VNCT.^2,XQEVNCT2.^2];
    SRR13=sum((LOGTotal_E_VNCT-Xt_INCT*tNCT).^2);
    
    
    %%%%%%%%   w'_(t+1)_I=savings_I(w',q_I,s,Entry|contest=0)  %%%%%%%%%
    %% Here, we want the deterministic function of savings conditional on being
    %% uncontested as a function of (w',q_I,s,tenure,Entry.
    
    Xw_INCT=[ones(size(LOGD_E_VNCT,1),1),LOGW_E_VNCT,PartyE_VNCT.*XS_EVNCT_(:,1),PartyE_VNCT.*XS_EVNCT_(:,2),TenureE_VNCT,XQEVNCT2,...
        LOGW_E_VNCT.^2,PartyE_VNCT.*(XS_EVNCT_(:,1)).^2,PartyE_VNCT.*(XS_EVNCT_(:,2)).^2,TenureE_VNCT.^2,XQEVNCT2.^2];
    SRR14=sum((LOGW_NXT_E_VNCT-Xw_INCT*wNCT).^2);
    
    
end

%%%%%%%%   Pr(win|w'_(t+1),w'_t,q_I,s,Entry)  %%%%%%%%%
if num==3
    dL2=[ones(length(Win_AC),1),LOGW_IAC,XQAC2,XSAC2.*PartyAC,LOGW_IAC.^2,XQAC2.^2,log(TenureAC+1)]*thetawin;
    SRR15=sum((Win_AC-dL2).^2)/var(Win_AC);
    
end
if num==4
    
    if (((E_VCTa(1:5,1).^2)'*gammaCT(1:5,1).^2)>=1)|(((E_VCTa(1:5,2).^2)'*gammaCT(1:5,2).^2)>=1)|(((E_VCTa(1:5,3).^2)'*gammaCT(1:5,3).^2)>=1)
        SRR9=0;
        for i=1:3
            SRR9=SRR9+10000000*max(0,((E_VCTa(1:5,i).^2)'*gammaCT(1:5,i).^2-1));
        end
    end
    if (((E_VCTa(6:10,1).^2)'*gammaCT(6:10,1).^2)>=1)|(((E_VCTa(6:10,2).^2)'*gammaCT(6:10,2).^2)>=1)|(((E_VCTa(6:10,3).^2)'*gammaCT(6:10,3).^2)>=1)
        SRR10=0;
        for i=1:3
            SRR10=SRR10+10000000*max(0,((E_VCTa(6:10,i).^2)'*gammaCT(6:10,i).^2-1));
        end
    end
    if (((E_VCTa(11:15,1).^2)'*gammaCT(11:15,1).^2)>=1)|(((E_VCTa(11:15,2).^2)'*gammaCT(11:15,2).^2)>=1)|(((E_VCTa(11:15,3).^2)'*gammaCT(11:15,3).^2)>=1)
        SRR11=0;
        for i=1:3
            SRR11=SRR11+10000000*max(0,((E_VCTa(11:15,i).^2)'*gammaCT(11:15,i).^2-1));
        end
    end
end
if num==1
    SRR=SRR1+SRR3;
elseif num==2
    SRR=SRR2;
elseif num==3
    SRR=SRR15;
elseif num==4
    SRR=SRR9+SRR10+SRR11;
elseif num==5
    SRR=SRR12+SRR13+SRR14;
end
if SRR<Sofarbest
    Sofarbest=SRR;
    bestiter=iterate;
end

iterate=iterate+1;
thetain_=thetain;

if bestiter==iterate-1
    if num==1
        save Best1.txt thetain_ SRR bestiter -ASCII
    elseif num==2
        save Best2.txt thetain_ SRR bestiter -ASCII
    elseif num==3
        save Best3.txt thetain_ SRR bestiter -ASCII
    elseif num==4
        save Best4.txt thetain_ SRR bestiter -ASCII
    elseif num==5
        save Best5.txt thetain_ SRR bestiter -ASCII
    end
end
