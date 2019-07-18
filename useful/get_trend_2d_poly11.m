function [ fit_data ] = get_trend_2d_poly11( X1, X2, Y, robust_flag, Y_weights )
%A funcion to compare a number of data sets with one other.

% *INPUT*
%           X1: 2D ARRAY - rows contain the independant data vectors
%
%           X2: 2D ARRAY - rows contain the independant data vectors
%
%           X1: 2D ARRAY - rows contain the dependant data vectors
%
% *OUTPUT*
%           fit_data: STRUCTURE -  has the fitting parameters and the
%           goodnes of fit statistics, as output from fit.m.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set up the output structure
ndata = length(Y(:, 1)); % the number of data vectors
ldata = length(Y(1, :)); % the length of the data vectors

fit_data.type = 'poly11';
if nargin < 4 || robust_flag == 0
    fit_data.robust = 0;
else
    fit_data.robust = 'bisquare';
end
if nargin < 5 || isempty(Y_weights)
    fit_data.weights_flag = 0;
else
    fit_data.weights_flag = 1;
end
fit_data.p00 = nan(ndata, 3); % for the zeroth-order coefficient and the confidence intervals
fit_data.p10 = nan(ndata, 3); % for the X1 coefficiant and confidence intervals
fit_data.p01 = nan(ndata, 3); % for the X2 coefficient and confidence intervals
fit_data.sse = nan(ndata, 1);
fit_data.rsquare = nan(ndata, 1);
fit_data.dfe = nan(ndata, 1);
fit_data.adjrsquare = nan(ndata, 1);
fit_data.rmse = nan(ndata, 1);

%% loop over the data and get the fits

for i = 1: ndata
    if isvector(X1)
        X1i = X1';
    else
        X1i = X1(i, :)';
    end
    if isvector(X2)
        X2i = X2';
    else
        X2i = X2(i, :)';
    end
    jbad = find(isnan(X1i) | isnan(X2i) | isnan(Y(i, :)') | isnan(Y_weights(i, :)') ); % find any nan data. want to exclude it
    if length(jbad) < ldata - 3 % need at least 4 points for a fit
        
        if fit_data.robust ~= 0
            if fit_data.weights_flag == 0
                [cf, gof] = fit([X1i, X2i], Y(i, :)', 'poly11', 'Exclude', jbad, 'Robust','Bisquare');
            else
%                 disp('using weights')
                [cf, gof] = fit([X1i, X2i], Y(i, :)', 'poly11', 'Exclude', jbad, 'Robust','Bisquare', 'Weights', Y_weights(i, :)');
            end
        else
            if fit_data.weights_flag == 0
                [cf, gof] = fit([X1i, X2i], Y(i, :)', 'poly11', 'Exclude', jbad);
            else
                [cf, gof] = fit([X1i, X2i], Y(i, :)', 'poly11', 'Exclude', jbad, 'Weights', Y_weights(i, :)');
            end
        end
        coeffi = coeffvalues(cf); % the fit coefficients
        confi = confint(cf); % the confidence intervals of the fit coefficients
        fit_data.p00(i, :) = [coeffi(1), confi(:, 1)'];
        fit_data.p10(i, :) = [coeffi(2), confi(:, 2)'];
        fit_data.p01(i, :) = [coeffi(3), confi(:, 3)'];
        fit_data.sse(i) = gof.sse;
        fit_data.rsquare(i) = gof.rsquare;
        fit_data.dfe(i) = gof.dfe;
        fit_data.adjrsquare(i) = gof.adjrsquare;
        fit_data.rmse(i) = gof.rmse;
    else
        warning('not enough data for a fit')
    end
end

%
end

