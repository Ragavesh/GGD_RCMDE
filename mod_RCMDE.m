function Out_RCMDE=mod_RCMDE(x,m,c,tau,Scale,beta)
%
% This function calculates the refined composite multiscale dispersion entropy (RCMDE) of a univariate signal x
%
% Inputs:
%
% x: univariate signal - a vector of size 1 x N (the number of sample points)
% m: embedding dimension
% c: number of classes (it is usually equal to a number between 3 and 9 - we used c=6 in our studies)
% tau: time lag (it is usually equal to 1)
% Scale: number of scale factors
%
%Outputs:
%
% Out_RCMDE: a vector of size 1 * Scale - the RCMDE of x
%
%
% Ref:
% [1] H. Azami, M. Rostaghi, D. Abasolo, and J. Escudero, "Refined Composite Multiscale Dispersion Entropy and its Application to Biomedical
% Signals", IEEE Transactions on Biomedical Engineering, 2017.
% [2] M. Rostaghi and H. Azami, "Dispersion Entropy: A Measure for Time-Series Analysis", IEEE Signal Processing Letters. vol. 23, n. 5, pp. 610-614, 2016.
%
% If you use the code, please make sure that you cite references [1] and [2].
%
% Hamed Azami and Javier Escudero Rodriguez
% Emails: hamed.azami@ed.ac.uk and javier.escudero@ed.ac.uk
%
%  20-January-2017
%%

Out_RCMDE=NaN*ones(1,Scale);

Out_RCMDE(1)=DisEn_GCDF(x,m,c,tau,beta);

sigma=std(x);
mu=mean(x);

for j=2:Scale
    pdf=[];
    for jj=1:j
        xs = Multi(x(jj:end),j);
        [DE, T_pdf]=mod_DisEn_NCDF_ms(xs,m,c,mu,sigma,tau,beta);
        pdf=[pdf ; T_pdf];
    end
    pdf=mean(pdf,1);
    pdf=pdf(pdf~=0);
    Out_RCMDE(j)=-sum(pdf .* log(pdf));
end


function M_Data = Multi(Data,S)

%  generate the consecutive coarse-grained time series
%  Input:   Data: time series;
%           S: the scale factor
% Output:
%           M_Data: the coarse-grained time series at the scale factor S

L = length(Data);
J = fix(L/S);

for i=1:J
    M_Data(i) = mean(Data((i-1)*S+1:i*S));
end