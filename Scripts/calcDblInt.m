function sol = calcDblInt(h, k, nu)

% Truncated form of infinite interval
innerintegrand = @(x,y) normcdf(h./sqrt(nu*(1./x + 1./y)),0,1).*chi2pdf(x,nu);
outerintegrand = @(y) integral(@(x)innerintegrand(x,y), 0, Inf)^(k-1).*chi2pdf(y,nu);
sol = integral(outerintegrand, 0, Inf);

% Change of variables: x = s/(1-s), y = t/(1-t)
%innerintegrand = @(s,t) normcdf(h./sqrt(nu.*((1-s)./s + (1-t)./t)),0,1).*chi2pdf(s./(1-s),nu).*(1./((1-s).^2));
%outerintegrand = @(t) integral(@(s)innerintegrand(s,t),0,1)^(k-1).*chi2pdf(t./(1-t),nu).*(1./((1-t).^2));
%sol = integral(outerintegrand, 0, 1);