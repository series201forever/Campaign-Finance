function C=Actions(num,ST,Same, q_I, X_Knot1,rtotd, Ten,w_I,demo,Presdum,Midterm,shockpartisan, shockump, presseq, Entry, coefentry, E_VCTa, E_VCTt, gammact, coefspend, coeffund,coefsave, dF_gamma_ct, dF_total_ct, dF_nxt_nxt,Winrnd,coefprobwin,Ret,Betapartisan,Betaump)
%options = optimset();
%%%%%%%%%%%%%%%%%%%%%%%%%
%Make sure that specification chosen at the first stage and "entryprob",
%"X_(contested)", "X_(uncontested)", "win" matches.

global N;  % the number of simulations
global interest

options=optimoptions('fsolve','display','none');


C=zeros(5,10,N);
for i=1:N
    j=1;
    ST_=ST;
    Ten_=Ten;
    Same_=Same;
    w_I_=log(1+interest)+w_I;
    while j<Ret(i)
%       Update state
        ST_=[Betaump'*[1;ST_(1,1)]+shockump(j,i);Betapartisan'*[1;ST_(2,1)]+shockpartisan(j,i)];
        %Update tenure
        Ten_=Ten_+1;
        %Update same
        if Same_==1
            Same_=Same_-2*presseq(j,i);
        else
            Same_=Same_+2*presseq(j,i);
        end
        
        %Classify tenure
        if (Ten_<=3)
            iK=1;
        elseif (Ten_>3)&&(Ten_<=7)
            iK=2;
        else
            iK=3;
        end
        
        %Calculate probability of entry given updated state
        %%%%%%%%%%%%%%%Need to match with first stage specification
        if length(coefentry)==13
            entryprob=[1,w_I_,ST_(1,1)*Same_,ST_(2,1)*demo,log(Ten_+1),X_Knot1]*coefentry;
        else
            entryprob=[1,w_I_,ST_(1,1)*Same_,ST_(2,1)*demo,log(Ten_+1),X_Knot1,X_Knot1.*(w_I_*ones(1,8))]*coefentry;
        end
            
        if (Entry(j,i)<entryprob)  %% CONTEST.
            C(1,j,i)=1;  % Contest==1
            qtile1=dF_gamma_ct(j,i);
            %%%%%%%%%%%%%%%Need to match with first stage specification
            X_=[q_I ,w_I_,ST_(1,1),ST_(2,1).*demo,Same_,Presdum,Midterm]';
            %%%%%%%%%%%%%%%
            A1sp=(sqrt(1-(E_VCTa(1:8,iK).^2)'*(gammact(1:8,iK).^2))+E_VCTa(2:8,iK)'*(X_-E_VCTt(2:8,iK)))^2;
            A2sp=2*E_VCTa(1,iK)*(sqrt(1-(E_VCTa(1:8,iK).^2)'*(gammact(1:8,iK).^2))+E_VCTa(2:8,iK)'*(X_-E_VCTt(2:8,iK)));

            Jsp=(sqrt(1-(E_VCTa(1:8,iK).^2)'*(gammact(1:8,iK).^2))+E_VCTa(2:8,iK)'*(X_-E_VCTt(2:8,iK)))^2+E_VCTa(1,iK)^2*gammact(1,iK)^2;
%             spending=fminsearch(@(X) PHI([E_VCTa(1:5,iK);sqrt(1-(E_VCTa(1:5,iK).^2)'*(gammact(1:5,iK).^2))],E_VCTt(1:5,iK)...
%                    ,gammact(1:5,iK),[q_I,w_I_,ST_'*thetaS*Party,log(1+Ten_)]',X,qtile1),E_VCTt(1,iK));
            spending=fsolve(@(X) PHI_improve(A1sp,A2sp,Jsp,(1/sqrt(2*pi)),E_VCTa(1,iK),E_VCTt(1,iK)...
                   ,gammact(1,iK),X,qtile1),11.5,options); 

               
              %% Computes F^(-1)(dF_gamma_ct) where
                                                             %% F^(-1) is the inverse cdf of
                                                             %% LOGD_NXT_E_VCT
                                                                                                        
    %         spending2=fminunc(@(X) PHI([E_VCTa(1:5,iK);sqrt(1-(E_VCTa(1:5,iK).^2)'*(gamma(1:5).^2))],E_VCTt(1:5,iK)...
    %             ,gamma(1:5,iK),[q_I,w_I_,ST_'*thetaS*Party,log(1+Ten_)]',X,dF_gamma_ct(j,i)),6,options);    %% Computes F^(-1)(dF_gamma_ct) where
    %                                                                                                %% F^(-1) is the inverse cdf of
    %                                                                                                %% LOGD_NXT_E_VCT
    %                                                                                                
    %          Check1=PHI([E_VCTa(1:5,iK);sqrt(1-(E_VCTa(1:5,iK).^2)'*(gammact(1:5).^2))],E_VCTt(1:5,iK)...
    %              ,gammact(1:5,iK),[q_I,w_I_,ST_'*thetaS*Party,log(1+Ten_)]',spending,dF_gamma_ct(j,i))   %% Check
    %         
    %        
            qtile2=dF_total_ct(j,i);
            A1tot=(sqrt(1-(E_VCTa(9:16,iK).^2)'*(gammact(9:16,iK).^2))+E_VCTa(10:16,iK)'*(X_-E_VCTt(10:16,iK)))^2;
            A2tot=2*E_VCTa(9,iK)*(sqrt(1-(E_VCTa(9:16,iK).^2)'*(gammact(9:16,iK).^2))+E_VCTa(10:16,iK)'*(X_-E_VCTt(10:16,iK)));

            Jtot=(sqrt(1-(E_VCTa(9:16,iK).^2)'*(gammact(9:16,iK).^2))+E_VCTa(10:16,iK)'*(X_-E_VCTt(10:16,iK)))^2+E_VCTa(9,iK)^2*gammact(9,iK)^2;
%             total=fminsearch(@(X) PHI([E_VCTa(6:10,iK);sqrt(1-(E_VCTa(6:10,iK).^2)'*(gammact(6:10,iK).^2))],E_VCTt(6:10,iK)...
%                     ,gammact(6:10,iK),[q_I,w_I_,ST_'*thetaS*Party,log(1+Ten_)]',X,qtile2),E_VCTt(6,iK));   %% Computes F^(-1)(dF_total_ct) where
            total=fsolve(@(X) PHI_improve(A1tot,A2tot,Jtot,(1/sqrt(2*pi)),E_VCTa(9,iK),E_VCTt(9,iK)...
                   ,gammact(9,iK),X,qtile2),11.5,options);     %% Computes F^(-1)(dF_gamma_ct) where
                                                             %% F^(-1) is the inverse cdf of
                                                             %% LOGTotal_NXT_E_VCT
    %         total2=fminunc(@(X) PHI([E_VCTa(6:10,iK);sqrt(1-(E_VCTa(6:10,iK).^2)'*(gamma(6:10).^2))],E_VCTt(6:10,iK)...
    %             ,gamma(6:10,iK),[q_I,w_I_,ST_'*thetaS*Party,log(1+Ten_)]',X,dF_total_ct(j,i)),7,options);                                                                                               
    %          CHeck2=PHI([E_VCTa(6:10,iK);sqrt(1-(E_VCTa(6:10,iK).^2)'*(gamma(6:10).^2))],E_VCTt(6:10,iK)...
    %              ,gamma(6:10,iK),[q_I,w_I_,ST_'*thetaS*Party,log(1+Ten_)]',total,dF_total_ct(j,i))

            qtile3=dF_nxt_nxt(j,i);
            A1nxt=(sqrt(1-(E_VCTa(17:24,iK).^2)'*(gammact(17:24,iK).^2))+E_VCTa(18:24,iK)'*(X_-E_VCTt(18:24,iK)))^2;
            A2nxt=2*E_VCTa(17,iK)*(sqrt(1-(E_VCTa(17:24,iK).^2)'*(gammact(17:24,iK).^2))+E_VCTa(18:24,iK)'*(X_-E_VCTt(18:24,iK)));

            Jnxt=(sqrt(1-(E_VCTa(17:24,iK).^2)'*(gammact(17:24,iK).^2))+E_VCTa(18:24,iK)'*(X_-E_VCTt(18:24,iK)))^2+E_VCTa(17,iK)^2*gammact(17,iK)^2;
            
%             savings=fminsearch(@(X) PHI([E_VCTa(11:15,iK);sqrt(1-(E_VCTa(11:15,iK).^2)'*(gammact(11:15,iK).^2))],E_VCTt(11:15,iK),...
%                  gammact(11:15,iK),[q_I,w_I_,ST_'*thetaS*Party,log(1+Ten_)]',X,qtile3),E_VCTt(11,iK));     %% Computes F^(-1)(dF_nxt_nxt) where
                                                                                                         %% F^(-1) is the inverse cdf of
                                                                                                         %% LOGWNXT_NXT_E_VCT
            savings=fsolve(@(X) PHI_improve(A1nxt,A2nxt,Jnxt,(1/sqrt(2*pi)),E_VCTa(17,iK),E_VCTt(17,iK)...
                   ,gammact(17,iK),X,qtile3),E_VCTt(17,iK),options);
    %         savings2=fminunc(@(X) PHI([E_VCTa(11:15,iK);sqrt(1-(E_VCTa(11:15,iK).^2)'*(gamma(11:15).^2))],E_VCTt(11:15,iK),...
    %             gamma(11:15,iK),[q_I,w_I_,ST_'*thetaS*Party,log(1+Ten_)]',X,dF_nxt_nxt(j,i)),6,options);
    %            Check3=PHI([E_VCTa(11:15,iK);sqrt(1-(E_VCTa(11:15,iK).^2)'*(gammact(11:15).^2))],E_VCTt(11:15,iK),...
    %               gammact(11:15,iK),[q_I,w_I_,ST_'*thetaS*Party,log(1+Ten_)]',savings,dF_nxt_nxt(j,i))    
    
%     thetawin(1,1)
%     thetawin(2,1)*w_I_
%     thetawin(3,1)*q_I
%     thetawin(4,1)*ST_'*thetaS*Party
%     thetawin(5,1)*w_I_^2
%     thetawin(6,1)*q_I^2
%     thetawin(7,1)*log(1+Ten_)
%             [1,w_I_,q_I,ST_'*thetaS*Party,w_I_^2,q_I^2,log(1+Ten_)]*thetawin
%             if num==34
%                 show1=[1,w_I_,q_I,ST_'*thetaS*Party,w_I_^2,q_I^2,log(1+Ten_)]
%                 show2=Winrnd(j,i)
%                 show3=thetawin
%                 ai=i
%                 jay=j
%             end



%%%%%%%%%%%%%%%Need to match with first stage specification 
             win=(Winrnd(j,i)<=[1,w_I_,ST_(1,1)*Same_,ST_(1,1).^2.*Same_,ST_(2,1)*demo,w_I_.^2,X_Knot1,log(Ten_+1),(log(Ten_+1)).^2,...
                 X_Knot1.*((ST_(1,1)*Same_)*ones(1,8)),X_Knot1.*((ST_(2,1)*demo)*ones(1,8)),X_Knot1.*(log(Ten_+1)*ones(1,8)),...
                 w_I_.^3]*coefprobwin); %% win if Winrnd<probability of winning                 
             C(2,j,i)=spending;
             C(3,j,i)=total;
             C(4,j,i)=savings;
             C(5,j,i)=win;
             w_I_=log(1+interest)+savings;
             j=j+1;
            % j=(1-win)*Ret(1,i)+win*j;
        else
            C(1,j,i)=0;
            %%%%%%%%%%%%%%%Need to match with first stage specification
            X_1=[1,w_I_,ST_(1,1),ST_(2,1).*demo,Same_,Same_*ST_(1,1),Ten_,X_Knot1,w_I_^2,(ST_(1,1))^2,(ST_(1,1))^2.*Same_, (ST_(2,1)^2)*demo,Ten_^2,...
                 rtotd.*w_I_,X_Knot1.*((ST_(1,1)*Same_)*ones(1,8)),X_Knot1.*((ST_(2,1)*demo)*ones(1,8)),X_Knot1.*(Ten_*ones(1,8)),...
                 Presdum,Presdum.*Same_,ST_(1,1).*Presdum.*Same_,Midterm, w_I_.*Ten_,rtotd.^2.*w_I_];
   
            X_2=[1,w_I_,ST_(1,1),ST_(2,1).*demo,Same_,Same_*ST_(1,1),Ten_,X_Knot1,w_I_^2,(ST_(1,1))^2,(ST_(1,1))^2.*Same_, (ST_(2,1)^2)*demo,Ten_^2,...
                 rtotd.*w_I_,X_Knot1.*((ST_(1,1)*Same_)*ones(1,8)),X_Knot1.*((ST_(2,1)*demo)*ones(1,8)),X_Knot1.*(Ten_*ones(1,8)),...
                 Presdum,Presdum.*Same_,ST_(1,1).*Presdum.*Same_,Midterm,rtotd.^2.*w_I_];
              

            spending=X_1*coefspend;
            total=X_1*coeffund;
            savings=X_2*coefsave;
             C(2,j,i)=spending;
             C(3,j,i)=total;
             C(4,j,i)=savings;
             C(5,j,i)=1;
            w_I_=log(1+interest)+savings;
            j=j+1;
        end
  
    end
end
