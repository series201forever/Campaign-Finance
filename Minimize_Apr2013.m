function [SRR]=Minimize_Apr2013(thetain,num)

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
global LOGD_E_VCTwnxt
global LOGTotal_E_VCTwnxt
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

E_VCTa(:,1)=thetain(1:8,1);
E_VCTt(:,1)=thetain(9:16,1);
gammaCT(:,1)=thetain(17:24,1);

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



 i=1;
 
 
 if num<30
     if num==11||num==21
         TAU=find((TenureE_VCT<=Cutoff(2,1))&(TenureE_VCT>Cutoff(1,1))); % TAU indexes tenure;
     end
     if num==12||num==22
         TAU=find((TenureE_VCT<=Cutoff(3,1))&(TenureE_VCT>Cutoff(2,1))); % TAU indexes tenure;
     end
     if num==13||num==23
         TAU=find((TenureE_VCT<=Cutoff(4,1))&(TenureE_VCT>Cutoff(3,1))); % TAU indexes tenure;
     end
     
     
     
     
     if num==11||num==12||num==13
         %%%%%%%%   F(Spending|w_I,q_I,unemp,partisan,same,contest=1)  %%%%%%%%%
         X=LOGD_E_VCT(TAU);
     end
     
     if num==21||num==22||num==23
         %%%%%%%%   F(total_I|w',q_I,s,Entry,contest=1)  %%%%%%%%%
         X=LOGTotal_E_VCT(TAU);
     end
     SRR9=0;
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
     
     
     
     E_VCTa01iP=(max(0,1/(((2*pi)^4)*prod(gammaCT(1:8,i)))-(E_VCTa(1:8,i).^2)'*gammaCT(1:8,i).^2))^(1/2);%/(((2*pi)^4)*prod(gammaCT(1:8,i)))
     E_VCTa01iM=(-1)*(max(0,1/(((2*pi)^4)*prod(gammaCT(1:8,i)))-(E_VCTa(1:8,i).^2)'*gammaCT(1:8,i).^2))^(1/2);
     loglikeE_VDI1iP(:,1)=log((E_VCTa01iP+E_VCTa(1,i)*(X-E_VCTt(1,i))+E_VCTa(2,i)*(XQEVCTi2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTi-E_VCTt(3,i))+...
         E_VCTa(4,i)*(UnemploymentE_VCTi-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTi.*PartyE_VCTi-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTi-E_VCTt(6,i))+...
         E_VCTa(7,i)*(PresdumE_VCTi-E_VCTt(7,i))+E_VCTa(8,i)*(MidtermE_VCTi-E_VCTt(8,i))).^2);
     loglikeE_VDI1iM(:,1)=log((E_VCTa01iM+E_VCTa(1,i)*(X-E_VCTt(1,i))+E_VCTa(2,i)*(XQEVCTi2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTi-E_VCTt(3,i))+...
         E_VCTa(4,i)*(UnemploymentE_VCTi-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTi.*PartyE_VCTi-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTi-E_VCTt(6,i))+...
         E_VCTa(7,i)*(PresdumE_VCTi-E_VCTt(7,i))+E_VCTa(8,i)*(MidtermE_VCTi-E_VCTt(8,i))).^2);
     
     loglikeE_VDI2i(:,1)=(-1/2)*(((1/gammaCT(1,i))*(X-E_VCTt(1,i))).^2);
     
%     loglikeE_VDI3i=(-1)*length(X)*log(gammaCT(1,i));
     
     loglikeE_VDI4iP(:,1)=(-1)*log((E_VCTa01iP+E_VCTa(2,i)*(XQEVCTi2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTi-E_VCTt(3,i))+...
         E_VCTa(4,i)*(UnemploymentE_VCTi-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTi.*PartyE_VCTi-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTi-E_VCTt(6,i))+...
         E_VCTa(7,i)*(PresdumE_VCTi-E_VCTt(7,i))+E_VCTa(8,i)*(MidtermE_VCTi-E_VCTt(8,i))).^2+...
         (E_VCTa(1,i)^2)*(gammaCT(1,i)^2))-log(sqrt(2*pi)*gammaCT(1,i));%sqrt(2*pi)*
     
     loglikeE_VDI4iM(:,1)=(-1)*log((E_VCTa01iM+E_VCTa(2,i)*(XQEVCTi2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTi-E_VCTt(3,i))+...
         E_VCTa(4,i)*(UnemploymentE_VCTi-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTi.*PartyE_VCTi-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTi-E_VCTt(6,i))+...
         E_VCTa(7,i)*(PresdumE_VCTi-E_VCTt(7,i))+E_VCTa(8,i)*(MidtermE_VCTi-E_VCTt(8,i))).^2+...
         (E_VCTa(1,i)^2)*(gammaCT(1,i)^2))-log(sqrt(2*pi)*gammaCT(1,i));%sqrt(2*pi)*
     
     loglikeE_VDIi=max(sum(loglikeE_VDI1iP+loglikeE_VDI2i+loglikeE_VDI4iP),...
         sum(loglikeE_VDI1iM+loglikeE_VDI2i+loglikeE_VDI4iM));
     SRR9=SRR9+loglikeE_VDIi;
     
     SRR9=(-1)*SRR9;
     
     
     if (((E_VCTa(1:8,1).^2)'*gammaCT(1:8,1).^2)>=1/(((2*pi)^4)*prod(gammaCT(1:8,i))))
         
         SRR9=SRR9+10000000*max(0,((E_VCTa(1:8,i).^2)'*gammaCT(1:8,i).^2-1/(((2*pi)^4)*prod(gammaCT(1:8,i)))));
         
     end
%      if min(E_VCTa01iP+E_VCTa(1,i)*(X-E_VCTt(1,i))+E_VCTa(2,i)*(XQEVCTi2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTi-E_VCTt(3,i))+...
%              E_VCTa(4,i)*(UnemploymentE_VCTi-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTi.*PartyE_VCTi-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTi-E_VCTt(6,i))+...
%              E_VCTa(7,i)*(PresdumE_VCTi-E_VCTt(7,i))+E_VCTa(8,i)*(MidtermE_VCTi-E_VCTt(8,i)))<10^(-5)
%          SRR9=SRR9+10^10*max(0,min(E_VCTa01iP+E_VCTa(1,i)*(X-E_VCTt(1,i))+E_VCTa(2,i)*(XQEVCTi2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTi-E_VCTt(3,i))+...
%              E_VCTa(4,i)*(UnemploymentE_VCTi-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTi.*PartyE_VCTi-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTi-E_VCTt(6,i))+...
%              E_VCTa(7,i)*(PresdumE_VCTi-E_VCTt(7,i))+E_VCTa(8,i)*(MidtermE_VCTi-E_VCTt(8,i))));
%      end
     
     SRR=SRR9;
 end

    


if num>30
    
    SRR11=0;
    
    if num==31||num==311||num==321
        TAU=find((TenureE_VCTwnxt([1:183,187:409,411:521])<=Cutoff(2,1))&(TenureE_VCTwnxt([1:183,187:409,411:521])>Cutoff(1,1))); % TAU indexes tenure;
    end
    if num==32||num==312||num==322
        TAU=find((TenureE_VCTwnxt([1:183,187:409,411:521])<=Cutoff(3,1))&(TenureE_VCTwnxt([1:183,187:409,411:521])>Cutoff(2,1))); % TAU indexes tenure;
    end
    if num==33||num==313||num==323
        TAU=find((TenureE_VCTwnxt([1:183,187:409,411:521])<=Cutoff(4,1))&(TenureE_VCTwnxt([1:183,187:409,411:521])>Cutoff(3,1))); % TAU indexes tenure;
    end
    LOGW_NXT_E_VCTwnxti=LOGW_NXT_E_VCTwnxt([1:183,187:409,411:521]);
    LOGD_E_VCTwnxti=LOGD_E_VCTwnxt([1:183,187:409,411:521]);   
    LOGTotal_E_VCTwnxti=LOGTotal_E_VCTwnxt([1:183,187:409,411:521]);       
    if num<40
    X=LOGW_NXT_E_VCTwnxti(TAU);
    end
    if num==311||num==312||num==313
     X=LOGD_E_VCTwnxti(TAU);   
    end
    if num==321||num==322||num==323
     X=LOGTotal_E_VCTwnxti(TAU);
    end


        XQEVCTwnxti2=XQEVCTwnxt2([1:183,187:409,411:521]);
        XQEVCTwnxti2=XQEVCTwnxti2(TAU);
        LOGW_E_VCTwnxti=LOGW_E_VCTwnxt([1:183,187:409,411:521]);
        LOGW_E_VCTwnxti=LOGW_E_VCTwnxti(TAU);
        SameE_VCTwnxti=SameE_VCTwnxt([1:183,187:409,411:521]);
        SameE_VCTwnxti=SameE_VCTwnxti(TAU);
        UnemploymentE_VCTwnxti=XS_EVCTwnxt_([1:183,187:409,411:521],1);
        UnemploymentE_VCTwnxti=UnemploymentE_VCTwnxti(TAU);
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
        
        E_VCTa03iP=(max(0,1/(((2*pi)^4)*prod(gammaCT(1:8,i)))-(E_VCTa(1:8,i).^2)'*gammaCT(1:8,i).^2))^(1/2);
        E_VCTa03iM=(-1)*(max(0,1/(((2*pi)^4)*prod(gammaCT(1:8,i)))-(E_VCTa(1:8,i).^2)'*gammaCT(1:8,i).^2))^(1/2);
        
        loglikeE_VWnxtnxtI1iP(:,1)=log((E_VCTa03iP+E_VCTa(1,i)*(X-E_VCTt(1,i))+E_VCTa(2,i)*(XQEVCTwnxti2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTwnxti-E_VCTt(3,i))+...
            E_VCTa(4,i)*(UnemploymentE_VCTwnxti-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTwnxti.*PartyE_VCTwnxti-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTwnxti-E_VCTt(6,i))+...
            E_VCTa(7,i)*(PresdumE_VCTwnxti-E_VCTt(7,i))+E_VCTa(8,i)*(MidtermE_VCTwnxti-E_VCTt(8,i))).^2);
        
         loglikeE_VWnxtnxtI1iM(:,1)=log((E_VCTa03iM+E_VCTa(1,i)*(X-E_VCTt(1,i))+E_VCTa(2,i)*(XQEVCTwnxti2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTwnxti-E_VCTt(3,i))+...
             E_VCTa(4,i)*(UnemploymentE_VCTwnxti-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTwnxti.*PartyE_VCTwnxti-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTwnxti-E_VCTt(6,i))+...
             E_VCTa(7,i)*(PresdumE_VCTwnxti-E_VCTt(7,i))+E_VCTa(8,i)*(MidtermE_VCTwnxti-E_VCTt(8,i))).^2);
         
        loglikeE_VWnxtnxtI2i(:,1)=(-1/2)*(((1/gammaCT(1,i))*(X-E_VCTt(1,i))).^2);
        
%        loglikeE_VWnxtnxtI3i=(-1)*length(X)*log(gammaCT(1,i));
        
        loglikeE_VWnxtnxtI4iP(:,1)=(-1)*log((E_VCTa03iP+E_VCTa(2,i)*(XQEVCTwnxti2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTwnxti-E_VCTt(3,i))+...
            E_VCTa(4,i)*(UnemploymentE_VCTwnxti-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTwnxti.*PartyE_VCTwnxti-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTwnxti-E_VCTt(6,i))+...
            E_VCTa(7,i)*(PresdumE_VCTwnxti-E_VCTt(7,i))+E_VCTa(8,i)*(MidtermE_VCTwnxti-E_VCTt(8,i))).^2+...
            (E_VCTa(1,i)^2)*(gammaCT(1,i)^2))-log(sqrt(2*pi)*gammaCT(1,i));
        
         loglikeE_VWnxtnxtI4iM(:,1)=(-1)*log((E_VCTa03iM+E_VCTa(2,i)*(XQEVCTwnxti2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTwnxti-E_VCTt(3,i))+...
             E_VCTa(4,i)*(UnemploymentE_VCTwnxti-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTwnxti.*PartyE_VCTwnxti-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTwnxti-E_VCTt(6,i))+...
             E_VCTa(7,i)*(PresdumE_VCTwnxti-E_VCTt(7,i))+E_VCTa(8,i)*(MidtermE_VCTwnxti-E_VCTt(8,i))).^2+...
             (E_VCTa(1,i)^2)*(gammaCT(1,i)^2))-log(sqrt(2*pi)*gammaCT(1,i));
         
        loglikeE_VWnxtnxtIi=max(sum(loglikeE_VWnxtnxtI1iP+loglikeE_VWnxtnxtI2i+loglikeE_VWnxtnxtI4iP),...
            sum(loglikeE_VWnxtnxtI1iM+loglikeE_VWnxtnxtI2i+loglikeE_VWnxtnxtI4iM));
        SRR11=SRR11+loglikeE_VWnxtnxtIi;
   
    SRR11=(-1)*SRR11;
    
    if (((E_VCTa(1:8,1).^2)'*gammaCT(1:8,1).^2)>=1/(((2*pi)^4)*prod(gammaCT(1:8,i)))) 
            SRR11=SRR11+10000000*max(0,((E_VCTa(1:8,i).^2)'*gammaCT(1:8,i).^2-1/(((2*pi)^4)*prod(gammaCT(1:8,i)))));
    end
%     if min(E_VCTa03iP+E_VCTa(1,i)*(X-E_VCTt(1,i))+E_VCTa(2,i)*(XQEVCTwnxti2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTwnxti-E_VCTt(3,i))+...
%             E_VCTa(4,i)*(UnemploymentE_VCTwnxti-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTwnxti.*PartyE_VCTwnxti-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTwnxti-E_VCTt(6,i))+...
%             E_VCTa(7,i)*(PresdumE_VCTwnxti-E_VCTt(7,i))+E_VCTa(8,i)*(MidtermE_VCTwnxti-E_VCTt(8,i)))<10^(-5)
%          SRR11=SRR11+10^10*max(0,min(E_VCTa03iP+E_VCTa(1,i)*(X-E_VCTt(1,i))+E_VCTa(2,i)*(XQEVCTwnxti2-E_VCTt(2,i))+E_VCTa(3,i)*(LOGW_E_VCTwnxti-E_VCTt(3,i))+...
%             E_VCTa(4,i)*(UnemploymentE_VCTwnxti-E_VCTt(4,i))+E_VCTa(5,i)*(PartisanE_VCTwnxti.*PartyE_VCTwnxti-E_VCTt(5,i))+E_VCTa(6,i)*(SameE_VCTwnxti-E_VCTt(6,i))+...
%             E_VCTa(7,i)*(PresdumE_VCTwnxti-E_VCTt(7,i))+E_VCTa(8,i)*(MidtermE_VCTwnxti-E_VCTt(8,i))));
%     end
    SRR=SRR11;

            
        
end

end

%%

%test1=[XQEVCTi2,LOGW_E_VCTi,SameE_VCTi,UnemploymentE_VCTi,PartisanE_VCTi.*PartyE_VCTi,PresdumE_VCTi,MidtermE_VCTi];
%testreg1=[test\LOGD_E_VCT(TAU),test\LOGTotal_E_VCT(TAU)]
%  test2=[XQEVCTwnxti2,LOGW_E_VCTwnxti,SameE_VCTwnxti,UnemploymentE_VCTwnxti,PartisanE_VCTwnxti.*PartyE_VCTwnxti,PresdumE_VCTwnxti, MidtermE_VCTwnxti];
%  testreg2=[test2\LOGD_E_VCTwnxti(TAU),test2\LOGTotal_E_VCTwnxti(TAU),test2\LOGW_NXT_E_VCTwnxti(TAU)]
% 
% [Est411,Est4311,Est421,Est4321]
% [Est412,Est4312,Est422,Est4322]
% [Est413,Est4313,Est423,Est4323]
% [Est431,Est432,Est433]
% 
% 
%  Est411(1)*Est411(3)-Est421(1)*Est421(3)
%  Est4311(1)*Est4311(3)-Est4321(1)*Est4321(3)
%  Est412(1)*Est412(3)-Est422(1)*Est422(3)
%  Est4312(1)*Est4312(3)-Est4322(1)*Est4322(3)
%  Est413(1)*Est413(3)-Est423(1)*Est423(3)
%  Est4313(1)*Est4313(3)-Est4323(1)*Est4323(3)