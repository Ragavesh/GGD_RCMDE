function [Out_DisEn, npdf]=DisEn_GCDF(x,m,nc,tau,beta)
%
% This function calculates dispersion entropy (DisEn) using normal cumulative distribution function (NCDF) with defined mean (mu) and standard deviation (sigma) values.
%
% Inputs:
%
% x: univariate signal - a vector of size 1 x N (the number of sample points)
% m: embedding dimension
% nc: number of classes (it is usually equal to a number between 3 and 9 - we used c=6 in our studies)
% tau: time lag (it is usually equal to 1)
%
% Outputs:
%
% Out_DisEn: scalar quantity - the DisEn of x
% npdf: a vector of length nc^m, showing the normalized number of disersion patterns of x
%
% Ref:
%
% [1] H. Azami, M. Rostaghi, D. Abasolo, and J. Escudero, "Refined Composite Multiscale Dispersion Entropy and its Application to Biomedical
% Signals", IEEE Transactions on Biomedical Engineering, 2017.
% [2] M. Rostaghi and H. Azami, "Dispersion Entropy: A Measure for Time-Series Analysis", IEEE Signal Processing Letters. vol. 23, n. 5, pp. 610-614, 2016.
%
% If you use the code, please make sure that you cite references [1] and [2].
%
% Hamed Azami, Mostafa Rostaghi, and Javier Escudero Rodriguez
% hamed.azami@ed.ac.uk, rostaghi@yahoo.com, and javier.escudero@ed.ac.uk
%
%  20-January-2017
%%

N=length(x);

mu_x=mean(x);
sigma_x=std(x);

rho=sqrt(gamma(1/beta)/gamma(3/beta))*sigma_x;
%pdf   
fy=(beta/(2*rho*gamma(1/beta)))*exp(-(abs(x - mu_x)/rho).^beta);
%cdf
c = 0.5*sign(x-mu_x);

Fy=(c+0.5) - c.*(gammainc((1.0/beta),(abs(x-mu_x)/rho)).^beta); 
%/(2*gamma(1/beta)); % not included 

y=Fy;
%y=gengauss_cdf(x,mu_x,sigma_x,beta); % verified with the current cdf,
% it matches the output

for i_N=1:N
    if y(i_N)==1
        y(i_N)=1-(1e-10);
    end
    
     if y(i_N)==0
        y(i_N)=(1e-10);
     end
end

z=round(y*nc+0.5);

all_patterns=[1:nc]';

for f=2:m
    temp=all_patterns;
    all_patterns=[];
    j=1;
    for w=1:nc
        [a,b]=size(temp);
        all_patterns(j:j+a-1,:)=[temp,w*ones(a,1)];
        j=j+a;
    end
end

for i=1:nc^m
    key(i)=0;
    for ii=1:m
        key(i)=key(i)*10+all_patterns(i,ii);
    end
end
embd2=zeros(N-(m-1)*tau,1);
for i = 1:m, 
   embd2=[z(1+(i-1)*tau:N-(m-i)*tau)]'*10^(m-i)+embd2;
   
end

pdf=zeros(1,nc^m);

for id=1:nc^m
    [R,C]=find(embd2==key(id));
    pdf(id)=length(R);
end

npdf=pdf/(N-(m-1)*tau);
p=npdf(npdf~=0);
Out_DisEn = -sum(p .* log(p));
