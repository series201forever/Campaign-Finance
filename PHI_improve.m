function [Phi]=PHI_improve(A1,A2,J,C,alpha, tau, gamma, X,qtile)

%% PHI first computes the conditional cdf Pr(x<X)=F(X | X_) given the Hermite pdf approximation
%% where alphas, taus and gammas are the parameters of the approximation.
%% X_ is the conditioning variable, and X is the upper end of the integral.
%% Then we compute (qtile-F(X|X_))^2.
%% If we solve min_(@)X PHI(.) we get the value of X s.t. F(X|X_)=qtile.

%% X_ is 6 X 1 vector of conditioning variables
%% X is the scalar, the point to which integrate over the pdf.

%% taus and gammas are 6 X 1 vector with the first element being the
%% variable that we take conditional expectation over.

%% alphas is 6 X 1 vec. with the last element being alpha0

%% get the value of alpha0=[1-(alphas(1:5,1).^2)'*gammas(1:5,1).^2]^(1/2) into
%% alphas(6,1).

A3=exp(-(X-tau)^2/(2*gamma^2));
A4=cdf('norm',(X-tau)/gamma,0,1);

Phi=A1*A4+A2*...
    (-gamma)*A3*C+alpha^2*((-1)*C*(X-tau)*gamma^2*...
    A3+gamma^3*A4);

    Phi=qtile-Phi/J;

end

%dPhi=(-2)*sqrt(Phi)*(1/(sqrt(2*pi)*gammas(1,1)))*((alphas(6,1)+alphas(1:5,
%1)'*([X;X_]-taus))^2)*exp(((X-taus(1,1))^2)*(-1/(2*gammas(1,1)^2)))/J;