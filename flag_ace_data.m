function [tanstruct_out, TRND] = flag_ace_data(tanstruct_in)
% flags the data and saves individual quality flag fields in the matlab
% structs in DIR.
%% this code has been modified from Patrick's orignial flagging code

%%
DIR = 'C:\Users\ryann\ACE\matdata\flag_testing\';
tanstruct = tanstruct_in;
% tanstruct_out = tanstruct_in;
% % % delete(gcp('nocreate')); parpool;
ALT = 0.5:149.5; Na = length(ALT);

DNU = [14390:14541 22060:25491 29680]; % these are the v2.2/v3.0 'DO NOT USE' occultations on the validation website

%-----------------------------problems found in v3.5-----------------------
DNU = [DNU 436321:446471]; disp('Hey Judith! NOV 2011 have been excludth. For some reason.')
BAR = three_bars_filter(DIR);
baddies = [DNU BAR];
baddies = [];
%--------------------------------------------------------------------------


% % %     if ISO == 0
% % %         species = {'H2O' 'O3' 'N2O' 'CO' 'CH4' 'NO' 'NO2' 'HNO3' 'HF' 'HCl' 'OCS' 'N2O5' 'ClONO2' 'HCN' 'CH3Cl' 'CF4' 'CCl2F2' 'CCl3F' 'COF2' 'COCl2' 'COClF' 'C2H6' 'C2H2' 'CHF2Cl' 'HCOOH' 'SF6' 'HNO4' 'H2O2' 'H2CO' 'CH3OH' 'CCl4' 'N2' 'O2' 'CFC113' 'HCFC141b' 'HCFC142b' 'CO2' 'T' 'P'};
% % %         atta = load([DIR 'attitude/attitude.mat'], 'data'); atta = atta.data;
% % %     else
% % %         species = {'T_iso' 'P_iso' 'n_iso' 'H2O_181' 'H2O_171' 'H2O_162' 'H2O_182' 'CO2_636' 'CO2_628' 'CO2_627' 'CO2_638' 'CO2_637' 'O3_668' 'O3_686' 'O3_667' 'O3_676' 'N2O_456' 'N2O_546' 'N2O_448' 'N2O_447' 'CO_36' 'CO_28' 'CO_27' 'CO_38' 'CH4_311' 'CH4_212' 'OCS_624' 'OCS_632' 'OCS_623'};
% % %         atta = load([DIR 'attitude_iso/attitude_iso.mat'], 'data'); atta = atta.data;
% % %     end
% % %     Ns = length(species);

% % %
% % %     atta = load([DIR 'ACE_v3p6_CO2.mat'],           'tanstruct'); atta = atta.tanstruct;
% % %
% % %     occ = atta.occultation; [~, I_bad] = intersect(occ, baddies);
% % %     srs = rem(occ,10);
% % %
% % %     mjd = atta.date_mjd; N = length(mjd); lat = atta.lat_tangent;
% % %     mjds = 0:ceil(max(mjd)-min(mjd));

occ = tanstruct.occultation; [~, I_bad] = intersect(occ, baddies);
% % % srs = rem(occ,10);
mjd = tanstruct.date_mjd; N = length(mjd); lat = tanstruct.lat_tangent;
mjds = 0:ceil(max(mjd)-min(mjd));


LATS = {[-90 -60]; [-60 0]; [0 60]; [60 90]}; Nl = length(LATS);

flag_baddies = zeros(Na, N);
flag_fill_v9 = flag_baddies;
flag_fill_v8 = flag_baddies;
flag_out_ext = flag_baddies;
flag_out_run = flag_baddies;
flag_no_stat = flag_baddies; flag_baddies(:, I_bad) = 1;

% % %         data = load([DIR X '/' X '.mat'], 'data'); data = data.data; if ~isequal(data.occ, atta.occ), error('Trouble on the occs.'), end
X = tanstruct.gas;
data = tanstruct; % changed to use input instead of loading files
if strcmp(X, 'T') || strcmp(X, 'P') || strcmp(X, 'n') || strcmp(X, 'T_iso')
    data.vmr_error( data.vmr_error == 0 ) = -888;
end
x = data.vmr; E = data.vmr_error;
t_F = x == -999;
x( t_F ) = nan; E( t_F ) = nan; flag_baddies( t_F ) = 0;
flag_fill_v9(t_F) = 1;

t_F = E == -888;
x( t_F ) = nan; flag_baddies( t_F ) = 0;
flag_fill_v8(t_F) = 1;

if strcmp(X, 'C2H6')
    t = any(data.vmr(21:150,:) > -999) & any(data.vmr_error(21:150,:) > -888);
    flag_baddies(:, t) = 1;
end

x(:, I_bad) = nan;
if strcmp(X, 'P') || strcmp(X, 'n')
    x = nan(size(x));
end

for h = 1:Na
    
    FOE = zeros(1, N); FOR = FOE; FNS = FOE; FYS = FOE;
    
    
    t = ~isnan(x(h,:)); Nt = sum(t);
    
    if Nt > 10
        %--------------------------------------------------------------
        % EXTREME
        o = occ(t); xt = x(h,t); MJD = mjd(t)-min(mjd); LAT = lat(t); %#ok<*PFBNS>
        
        t = abs(xt) < 1e4*median(abs(xt)); % med(abs(xt)) for data that has strong negative bias, or just noise (med can be very close to zero)
        [~, I] = intersect(occ, o(~t)); FOE(I) = 1;
        o = o(t); xt = xt(t); MJD = MJD(t); LAT = LAT(t);
        
        XXi = min(xt); XXf = max(xt);
        
        [FLAG, I_st] = extreme_outliers(xt, Nl, LAT, LATS, 0.025, XXi, XXf);
        [~, I] = intersect(occ, o(I_st)); FNS(I) = 1;
        
        t = FLAG == 1;
        [~, I] = intersect(occ, o(t)); FOE(I) = 1;
        flag_out_ext(h,:) = FOE;
        t = ~t;
        %--------------------------------------------------------------
        % REGULAR
        there_is_rejection = true; rej = true(1,Nl); first = true;
        RRR = 0;
        while there_is_rejection
            o = o(t); xt = xt(t); MJD = MJD(t); LAT = LAT(t);
            RRR = RRR+ 1;
            if first
                [REJ, FLAG, I_st] = regular_outliers(xt, rej, mjds, MJD, Nl, LAT, LATS);
                [~, I] = intersect(occ, o(I_st)); FYS(I) = 1;
                first = false;
            else
                [REJ, FLAG] = regular_outliers(xt, rej, mjds, MJD, Nl, LAT, LATS);
            end
            
            rej( REJ == 0 ) = false;
            t = FLAG == 1;
            
            if any(t)
                [~, I] = intersect(occ, o(t)); FOR(I) = 1;
                flag_out_run(h,:) = FOR;
                t = ~t;
            else
                %                             there_is_rejection = false;
            end
            there_is_rejection = false;
        end
        FNS( FYS == 1 ) = 0;
        %--------------------------------------------------------------
    else
        FNS(t) = 1;
    end
    flag_no_stat(h,:) = FNS;
end

F = flag_out_ext + flag_out_run; F(F>1) = 1;
F = flag_autocheck_cluster(x, F, occ, mjd, lat);

t = F == 0;
flag_out_ext(t) = 0;
flag_out_run(t) = 0;

data.flag_baddies = flag_baddies;
data.flag_fill_v9 = flag_fill_v9;
data.flag_fill_v8 = flag_fill_v8;
data.flag_out_ext = flag_out_ext;
data.flag_out_run = flag_out_run;
data.flag_no_stat = flag_no_stat;

%------------------------------------------------------
if strcmp(X, 'P') || strcmp(X, 'GLC')
    f_all = zeros(size(data.X));
else
    fbd = data.flag_baddies;
    
    if strcmp(X, 'N2O')
        ncid = read_netcdf([DIR_net 'ACEFTS_L2_v3p6_CH4.nc']);
        x_N = ncid.CH4; f_N = ncid.quality_flag;
        x_N(f_N>7) = nan; x_N(:, any(f_N >= 4 & f_N <= 6)) = nan;
        c_N = x_N(20,:);
        T_N = ncid.temperature(20,:);
        n_N = x(20,:);
        
        t_N = (c_N<10.6e-7 & n_N>1.6e-7) | (c_N<9.5e-7 & n_N>1.4e-7) & abs(lat) > 55 & T_N<196;
        fbd(:,t_N) = 1;
    end
    
    ff9 = data.flag_fill_v9;
    ff8 = data.flag_fill_v8;
    fox = data.flag_out_ext;
    foR = data.flag_out_run;
    fns = data.flag_no_stat;
    
    f_all = zeros(size(fbd));
    
    f_all(fns == 1) = 2;
    
    f_all(foR == 1) = 4;
    f_all(fox == 1) = 5;
    
    f_all(fbd == 1) = 7;
    f_all(ff8 == 1) = 8;
    f_all(ff9 == 1) = 9;
end
data.flag = f_all;
tanstruct_out = tanstruct_in;
tanstruct_out.quality_flags = data.flag;

% % %     x = data.X(:,Id); E = data.X_err(:,Id); f_all = f_all(:,Id);


% % %         save([DIR X '/' X '.mat'], 'data');

TRND = flag_autocheck(X, x, F, mjd, lat, DIR);
% % % delete(gcp);

function [REJ, FLAG, I_st] = regular_outliers(xt, rej, mjds, MJD, Nl, LAT, LATS)

pm = 7.5;
LIM = 10;

FLAG = zeros(size(xt)); ind = 1:length(FLAG);
REJ = zeros(size(rej));

I_st = []; % where data is being analyzed
for lt = 1:Nl
    if rej(lt)
        tlat = LAT > LATS{lt}(1) & LAT <= LATS{lt}(2);
        
        mn = nan(size(mjds)); sd = mn;
        mi = 0;
        for m = mjds, mi = mi+1;
            
            t = MJD > m-pm & MJD < m+pm & tlat;
            if sum(t) > 8
                xtt = xt(t);
                mn(mi) = median(xtt); sd(mi) = mean  ( abs( xtt - median(xtt) ) );
            end
        end
        
        xtt = xt(tlat); it = ind(tlat);
        
        M = interp1(mjds, mn, MJD(tlat));
        S = interp1(mjds, sd, MJD(tlat));
        
        f_ind = it(xtt > M + LIM*S | xtt < M - LIM*S);
        FLAG( f_ind ) = FLAG( f_ind ) + 1;
        
        I_st = unique([I_st it(~isnan(M))]);
        REJ(lt) = sum(f_ind);
    end
end

function [FLAG, I_st] = extreme_outliers(xt, Nl, LAT, LATS, LIM, XXi, XXf)

MNT = 0;
if XXi < 0
    MNT = abs(XXi);
elseif XXi == 0
    MNT = min(xt(xt~=0));
end
xt = log(xt + MNT*2); XXi = log(XXi+MNT*2); XXf = log(XXf+MNT*2);

stp = 1e-5;
xs = XXi-stp : stp : XXf+stp;

N = length(xt);

FLAG = zeros(size(xt)); ind = 1:N;
I_st = []; % where data is not being analyzed

for lt = 1:Nl
    tlat = LAT > LATS{lt}(1) & LAT <= LATS{lt}(2);
    if sum(tlat) > 150
        XX = xt(tlat); IT = ind(tlat);
        
        [XX, I] = sort(XX); Nt = length(XX);  Nx = length(XX); IT = IT(I);
        Iext = 6:Nt-5; % do the fit assuming that the highest/lowest 5 are extreme values (doesn't matter if they are or not)
        
        xx = XX(Iext);
        
        failure = true; deg = 4;
        while failure && deg > 0, warning off all
            try
                obj = gmdistribution.fit(xx', deg);
                failure = false;
            catch
                deg = deg-1; % reduce the number of gaussians you're fitting to
            end
        end, warning on all
        m = obj.mu;
        S = sqrt(squeeze(obj.Sigma(1,1,:)));
        p = obj.PComponents';
        
        y = gaussian(xs, m, S, p);
        z = sum(y,1) / (sum(sum(y,1)*stp));
        
        Z = cumsum(z) * stp * Nx;
        
        mmx = max( xs(Nx-Z >= LIM)); if isempty(mmx), mmx = max(xs); end
        mmn = min( xs(Z >= LIM)); if isempty(mmn), mmn = min(xs); end
        
        f_ind = IT(XX < mmn | XX > mmx);
        FLAG( f_ind ) = 1;
    else
        I_st = unique([I_st ind(tlat)]);
    end
end

function baddies = three_bars_filter(DIR)

d = load([DIR 'ACE_v3p6_CO2.mat'],           'tanstruct'); d = d.tanstruct;
% % % a = load([DIR 'attitude/attitude.mat'], 'data'); a = a.data;

% 2009 bar, mjds 55104-55148

m = d.date_mjd; o = d.occultation; x = d.vmr; x(x == -999) = nan; % changed to use CO2 instead of attitude
t = m > 55104 & m < 55148;

xt = x(42,t); ot = o(t); tt = xt < 3.55e-4 | xt > 3.95e-4;
ott = ot(tt);

% 2008 bar, mjds 54732-54798

m = d.date_mjd; o = d.occultation; x = d.vmr; x(x == -999) = nan;
t = m > 54732 & m < 54798;

xt = x(42,t); ot = o(t); tt = xt < 3.5e-4 | xt > 4.1e-4;
ott = [ot(tt) ott];

% 2006 bar, mjds 53730-53757

m = d.date_mjd; o = d.occultation; x = d.vmr; x(x == -999) = nan;
t = m > 53730 & m < 53757;

xt = x(42,t); ot = o(t); tt = xt < 3.45e-4 | xt > 3.9e-4;

ott = [ot(tt) ott];

knownfine = [127920,128481,130190,276740,277851,278351,279961,280691,280840,282151,282600,282670,334040,334221,335100,335301,336340];
knowndoubfed = [25990,87851,94121,96891,102880,111340,111911,114160,121640,124510,129761,129901,145581,157410,162471,169550,172351,176381,178351,178451,181850,181960,181980,181990,182060,223090,240540,255780,265451,276820,276860,276900,277140,277190,277250,277280,281381,281561,281761,286241,286251,286271,286291,286311,286331,286341,286351,286391,286401,286411,289770,289800,289820,289880,330830,330940,331070,331380,331390,334851,334911,335231,335300,335721,405780];

baddies = setdiff([ott knowndoubfed], knownfine);
%--------------------------------------------------------------------------

function y = gaussian(x, m, s, p)

if length(m) > 1 && all(size(m) == size(s))
    x = repmat(x, length(m), 1);
    m = repmat(m, 1, length(x));
    s = repmat(s, 1, length(x));
    if nargin == 4
        p = repmat(p, 1, length(x));
    else
        p = ones(size(x));
    end
end
y = p./s/sqrt(2*pi) .* exp(-0.5 * ( (x-m)./s ).^2 );
%--------------------------------------------------------------------------

function F = flag_autocheck_cluster(x, Fo, occ, mjd, lat)

ALT = 0.5:149.5;

mjd = mjd'; lat = lat'; occ = occ'; F = Fo;

ta = any(~isnan(x), 2);
x = x(ta,:); F = F(ta,:); alt = ALT(ta);
f = any(F);

mmax = floor(max(mjd)); mmin = ceil(min(mjd));

mpm = 7; lpm = 10;
mjds = mmin+mpm : mpm : mmax-mpm; Nm = length(mjds);
lats = -90-lpm : lpm : 90+lpm;    Nl = length(lats);

clust = nan(Nl, Nm);
for i = 1:Nm
    t = mjd >= mjds(i)-mpm & mjd < mjds(i)+mpm;
    for j = 1:Nl
        tt = t & lat >= lats(j)-lpm & lat < lats(j) + lpm;
        
        clust(j,i) = sum(f(tt));
    end
end
clust(clust == 0) = nan;

CLUST = clust(~isnan(clust));
m = nanmean(CLUST); s = nanstd(CLUST); lim = m+3*s;

mc = find(any(clust > lim)); Nmc = length(mc);
if Nmc == 1
    mc = [mc mc]; Nmc1 = true;
else
    Nmc1 = false;
end

if Nmc > 0
    mcd = [diff(mc) 0];
    I = find(mcd ~= 1); NI = length(I); ii = 1;
    if Nmc1, NI = 1; end
    
    for i = 1:NI
        mrange = mc(ii:I(i));
        ii = I(i)+1;
        
        t = mjd >= mjds(min(mrange))-mpm & mjd <= mjds(max(mrange))+mpm;
        
        lc = find( any(clust(:,mrange) > lim, 2) );
        if length(lc) == 1
            lc = repmat(lc,1,2); lc1 = true;
        else
            lc1 = false;
        end
        lcd = [diff(lc); 0];
        J = find(lcd ~= 1); NJ = length(J); jj = 1;
        if lc1, NJ = 1; end
        
        for j = 1:NJ
            lrange = lc(jj:J(j));
            jj = J(j)+1;
            
            tt = t & lat >= lats(min(lrange))-lpm & lat <= lats(max(lrange))+lpm;
            
            fr = F(:,tt); H = find(any(fr,2));
            xr = x(:,tt);
            for h = 1:2:length(H)
                Hh = alt >= alt(H(h))-2 & alt <= alt(H(h))+2;
                xh = xr(Hh,:);
                fh = fr(Hh,:);
                
                oc = repmat(occ(tt)',sum(Hh),1);
                at = repmat(alt(Hh)',1,length(xh));
                
                ttt = ~isnan(xh);
                xh = xh(ttt);
                fh = fh(ttt);
                oc = oc(ttt);
                at = at(ttt);
                
                m = median(xh); mad = median(abs(xh - m));
                
                th = fh > 0 & abs(xh-m) < 10*mad;
                oc = oc(th); at = at(th); N = length(oc);
                for n = 1:N
                    F( alt == at(n), occ == oc(n) ) = 0;
                end
            end
        end
    end
end
Fo(ta,:) = F;
F = Fo;
%--------------------------------------------------------------------------

function TRND = flag_autocheck(X, x, f, mjd, lat, DIR)

% check monthly mean, std, outlier density, and data count (after seasonal
% trend removed)

alt = 0.5:149.5; Na = length(alt);

LIM = 2.5;
fLIM = 2;

Ts.n2 = []; Ts.T = []; iT = 0;

year1 = 2004; lpm = 15;

warning('off','stats:statrobustfit:IterationLimit');

[y, m] = mjd2utc(mjd);

tf = any(f > 0);
x(f>0) = nan; x(:,tf) = nan;

ppv = 10^floor( nanmin(nanmin( log10(abs(x)) )) );
x = x/ppv;

[ymax, month] = mjd2utc(max(mjd));
doy = mjd2doy(mjd);
ystr = num2str(ymax);
if month < 10
    mstr = ['0' num2str(month)];
else
    mstr = num2str(month);
end

fid = fopen([DIR 'L2_check_log-' ystr '-' mstr '.txt'], 'a');
if ftell(fid) == 0
    fprintf(fid,'species\tissue\tlat\talt\r\n');
end

y = y - year1;
yp = y + doy'/365.24;
years = min(y):max(y); Ny = length(years);

lats = -90+lpm : 2*lpm : 90-lpm; Nl = length(lats);

fden = nan(1, Ny); dcnt = fden;
for j = 1:Ny
    t = m == month & y == years(j);
    if sum(t) > 0
        fden(j) = nanmean(tf(t)) * 100;
        dcnt(j) = sum(~tf(t));
    end
end
mfd = nanmean(fden); sfd = nanstd(fden); f_new = fden(Ny);
mdc = nanmean(dcnt); sdc = nanstd(dcnt); d_new = dcnt(Ny);

if f_new > mfd+fLIM*sfd
    fprintf(fid, [X '\tflag den\t-\t-\r\n']);
end
if d_new < mdc-fLIM*sdc
    fprintf(fid, [X '\tdata cnt\t-\t-\r\n']);
end

slp = nan(Na, Nl); slpe = slp;
for i = 1:Nl
    tl = lat > lats(i)-lpm & lat <= lats(i)+lpm & m == month;
    
    xslp = nan(Na, Ny); yps = nan(1,Ny);
    for j = 1:Ny
        ty = tl & y == years(j);
        xslp(:,j) = nanmean(x(:,ty), 2); yps(:,j) = nanmean(yp(ty));
    end
    
    for h = 1:Na
        t = ~isnan(xslp(h,:)); n = sum(t);
        
        if n > 10
            n2 = n-2;
            tT = Ts.n2 == n2;
            if any(tT)
                T = Ts.T(tT);
            else
                iT = iT+1;
                T = t_test_calc(n2, 0.99);
                Ts.n2(iT) = n2; Ts.T(iT) = T;
            end
            
            [p,stat] = robustfit(yps(t), xslp(h,t)); slp(h,i) = p(2)*ppv; slpe(h,i) = T*stat.se(2)*ppv;
            p = flip( polyfit(yps(t), xslp(h,t), 1) );
            
            x(h,tl) = x(h,tl) - polyval(flip(p), yp(tl));
        elseif n > 0
            x(h,tl) = x(h,tl) - nanmean(x(h,tl));
        end
    end
    
    xslp = nan(Na, Ny); sslp = xslp;
    for j = 1:Ny
        ty = tl & y == years(j);
        
        xslp(:,j) = nanmean(x(:,ty), 2);
        sslp(:,j) = nanstd(x(:,ty), 0, 2);
    end
    
    xt = x(:,tl);
    m_all = nanmean(xt,2); s_all = nanstd(xt,0,2); m_new = xslp(:,Ny);
    
    smean = nanmean(sslp,2); sstd = nanstd(sslp,0,2); s_new = sslp(:,Ny);
    
    if any( abs(m_new - m_all) > LIM*s_all )
        I = find( abs(m_new - m_all) > LIM*s_all )';
        for hi = I
            fprintf(fid, [X '\tmean lim\t' num2str(lats(i)) '\t' num2str(alt(hi)) '\r\n']);
        end
    end
    
    if any( abs(s_new - smean) > LIM*sstd )
        I = find( abs(s_new - smean) > LIM*sstd )';
        for hi = I
            fprintf(fid, [X '\tstdv lim\t' num2str(lats(i)) '\t' num2str(alt(hi)) '\r\n']);
        end
    end
end

slp(abs(slp) < slpe) = nan;
if any(any(~isnan(slp)))
    figure, contourf(lats, alt, slp), colorbar, title([X ', trend above error'])
    saveas(gcf,[DIR X '-trend-' ystr '-' mstr '.fig']), close all
end
TRND.linear_trend = slp;
TRND.linear_trend_err = slpe;
TRND.lat = lats;
TRND.alt = alt;
% % save([DIR X '-trend-' ystr '-' mstr '.mat'], 'TRND');

fclose(fid);
%--------------------------------------------------------------------------

function T = t_test_calc(df, conf, t_step)

if nargin < 3
    t_step = 0.001;
end
t = -10: t_step : 10;
t2 = t.^2;

if df <= 300
    if ~rem(df,2)
        Gn = prod( df-1 : -2 : 3 ); Gd = prod( [2 sqrt(df) df-2 : -2 : 2] );
    else
        Gn = prod( df-1 : -2 : 2 ); Gd = prod( [pi sqrt(df) df-2 : -2 : 3] );
    end
    G = Gn / Gd;
    f = G * (1+t2/df).^-( (df+1) / 2 );
else
    f = 1/sqrt(2*pi)*exp(-t2/2);
end

F = nan(size(t));
for i = 1:length(t)
    F(i) = sum(f(1:i)) * t_step;
end
[~,I] = min(abs(F-conf));

T = t(I);
%--------------------------------------------------------------------------