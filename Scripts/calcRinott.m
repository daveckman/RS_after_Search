function h = calcRinott(k, alpha, n0)

% Use simulation to estimate Rinott h

% df for chi-squared
nu = n0 - 1;

% Initialize for bisection search
sol = 0;
hlower = 0;
hupper = 30;
hmid = (hlower + hupper)/2;
prec = 0.0001;

% Do bisection until double integral is within prec of 1-alpha
while abs(sol - (1-alpha)) > prec  
    h = hmid;
    sol = calcDblInt(hmid, k, nu);
    if sol > 1-alpha
        hupper = hmid;
        hmid = (hupper + hlower)/2;
    elseif sol <= 1-alpha
        hlower = hmid;
        hmid = (hupper + hlower)/2;
    end   
end