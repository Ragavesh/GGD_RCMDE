% gengauss_cdf.m - compute Generalized Gaussian Cumulative Distribution Function.
%   See "A Practical Procedure to Estimate the Shape Parameter in the 
%   Generalized Gaussian Distribution", J. Armando Dominguez-Molina and 
%   R. M. Rodriguez-Dagnino, U. de Guanajuato.               
%
%               Vector Form of CDF !!!
%
%  Created by:  Jim Huntley,  08/28/06.
%
function [cdf] = gengauss_cdf(x,mu,sigma,beta)
tol = 1e-8;
trace = [];
warning off MATLAB:quad:MinStepSize;
minx = min(x);

A = sqrt(exp(2*log(sigma) + gammaln(1/beta) - gammaln(3/beta)));
coef = exp(-(log(2) + gammaln(1+1/beta) + log(A)));
%end
pdf = @(x)(coef .* exp(-abs((x-mu)./A).^beta)); 

% Integrate PDF to get CDF.
warning off MATLAB:quad:MinStepSize;
sz = size(x,2);
for jz = 1:sz
    %cdf(jz) = quad(@gengauss_pdf,minx,x(jz),tol,trace,mu,sigma,beta);
    cdf(jz) = integral(pdf,minx,x(jz),'AbsTol',1e-6,'RelTol',1e-3,'ArrayValued',true);
end
return
