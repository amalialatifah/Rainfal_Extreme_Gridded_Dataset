clear; clc;
C.stationDir = 'D:\Paper_revision\CH_Gabungan_Harian_Update';
C.paths.ERA5  = 'D:\Paper_revision\ERA5_daily_INDONESIA_ALL.nc';
C.paths.GPCP  = 'D:\Paper_revision\gabungan_GPCP_1998_2024.nc';
C.paths.SAOBS = 'D:\Paper_revision\rr_0.25deg_reg_v2.0_saobs.nc';
C.paths.CHIRPS = 'D:\2026_analysis_ENSO\CHIRPS';
C.paths.MSWEP = 'D:\Reanalysis_data\MSWEP';
C.outDir = 'D:\Paper_revision\Reviewer_3Tests_FINAL_Output';
C.startDate = datetime(2010,1,1);
C.endDate = datetime(2017,12,31);
C.minCompleteness = 0.70;  % STRICT > 70%, not >=
C.wetThreshold = 1;        % wet day if rain >= 1 mm/day
C.nBootstrap = 5000;       % block bootstrap; preliminary with 8 years
C.seed = 42;
C.tryMSWEPzip = true;
C.cleanUnzippedTemp = true;
C.preferNetcdfTime = true;
C.maxReasonableRain = 1500; % discard suspicious values above this; inspect first
C.preflightOnly = false; % already preflighted; true for metadata-only inspection   % set true to inspect files without expensive processing
% IMPORTANT: all supplied NetCDF products must already represent DAILY totals
% in mm/day (ERA5 in m/day is converted automatically).
% Check station timing in BMKG metadata; +/-1 day is a sensitivity test,
% not proof of the original observation-hour accumulation window.
C.products = {'ERA5','GPCP','SAOBS','CHIRPS','MSWEP'};
C.ERA5mmAlready = true; % data sample 7.93 suggests mm/day; VERIFY preprocessing
C.minAnnualCompleteness = 0.70;
C.minPairedStations = 3;
C.failOnIncompleteProduct = true;
C.useCache = true;  % reuse station, product, and per-ZIP MSWEP caches
C.forceRebuildStations = false;  % true once if station Excel values changed unexpectedly
C.forceRebuildProducts = {};     % e.g. {'MSWEP'} or {'GPCP'}
C.cacheVersion = 'reviewer_cache_v3_20261105';
C.cacheDir = fullfile(C.outDir,'CACHE');
if ~exist(C.outDir,'dir'), mkdir(C.outDir); end
if ~exist(C.cacheDir,'dir'), mkdir(C.cacheDir); end
rng(C.seed);
logfile = fullfile(C.outDir,'run_log.txt'); if exist(logfile,'file'),delete(logfile);end; diary(logfile); diary on;
fprintf('Started %s\n',datestr(now));
fprintf('CACHE folder: %s\n',C.cacheDir);
fprintf('Evaluation period: %s to %s\n',datestr(C.startDate),datestr(C.endDate));
% Preflight metadata. NetCDF variable names can vary across distributors.
files = struct();
for j=1:numel(C.products)
    p=C.products{j}; files.(p) = findInputFiles(C.paths.(p),p);
    fprintf('%s: %d data files identified\n',p,numel(files.(p)));
    if isempty(files.(p)), error('No source files found for %s',p); end
    if endsWith(lower(files.(p){1}),'.zip')
        fprintf('  ZIP contents: sample archive %s\n', files.(p){1});
        try, showZipContents(files.(p){1}); catch ME, warning('ZIP listing failed: %s',ME.message); end
    else
        try, inspectNetcdf(files.(p){1}); catch ME, warning('%s preflight: %s',p,ME.message); end
    end
end
if C.preflightOnly
    fprintf('PREFLIGHT COMPLETE (station files were not imported in this mode).\n');
    diary off; return;
end


% Station import cache: checkpoint immediately after the slow Excel parsing.
[stations,stationLog] = readStationsCached(C);
writetable(stationLog,fullfile(C.outDir,'station_import_log.csv'));
if isempty(stations), error('No stations with usable daily RR and coordinates. The provided example has headers but no daily observations.'); end
allDates = (C.startDate:caldays(1):C.endDate)';
N = numel(allDates); S = numel(stations);
obs = NaN(N,S);
for k=1:S
    [tf,loc] = ismember(dateshift(stations(k).dates,'start','day'),allDates);
    ok = tf & isfinite(stations(k).rain);
    % Duplicates are averaged at import; the unique dates below are stable.
    obs(loc(ok),k) = stations(k).rain(ok);
end
lat = [stations.lat]'; lon = [stations.lon]';
stid = string({stations.id})'; stname = string({stations.name})';
locInfo = table(stid,stname,lat,lon,'VariableNames',{'StationID','StationName','Latitude','Longitude'});
fullCoverage = sum(isfinite(obs),1)'/N;
locInfo.CompletenessFullPeriod = fullCoverage;
writetable(locInfo,fullfile(C.outDir,'station_coverage_full.csv'));
fprintf('Imported %d stations; %d with >70%% of %s..%s days.\n',S,sum(fullCoverage>C.minCompleteness),datestr(C.startDate),datestr(C.endDate));
% Exclude ineligible gauge stations BEFORE gridded extraction.
% Full 168-station diagnostics remain in station_coverage_full.csv.
keep = fullCoverage>C.minCompleteness;
if ~any(keep),error('No gauges satisfy strict >70%% completeness.');end
obs=obs(:,keep);
stations=stations(keep);
locInfo=locInfo(keep,:);
fullCoverage=fullCoverage(keep);
lat=lat(keep);lon=lon(keep);stid=stid(keep);stname=stname(keep);
S=sum(keep);
fprintf('Only %d qualifying gauges will be sampled from gridded datasets.\n',S);


% File-aware caches: include source size/timestamp and station coordinates.
% A new MSWEP ZIP invalidates ONLY the final MSWEP product cache;
% completed per-ZIP cache files are reused.
D = struct();
stationSignature = [string(stid) string(stname) string(lat) string(lon)];
for j=1:numel(C.products)
    p=C.products{j};
    cache=fullfile(C.cacheDir,['cache_' p '.mat']);
    signature = inputSignature(files.(p));
    useThis = C.useCache && ~any(strcmp(p,C.forceRebuildProducts));
    gotCache = false;
    if useThis && isfile(cache)
        try
            tmp=load(cache,'X','meta');
            gotCache = isfield(tmp,'meta') && ...
                isequal(tmp.meta.version,C.cacheVersion) && ...
                isequal(tmp.meta.files,signature) && ...
                isequal(tmp.meta.stations,stationSignature) && ...
                isequal(tmp.meta.dates,[C.startDate C.endDate]) && ...
                isequal(tmp.meta.size,[N S]) && ...
                isequal(tmp.meta.ERA5mmAlready,C.ERA5mmAlready) && ...
                isequal(size(tmp.X),[N S]);
            if gotCache
                X=tmp.X;
                fprintf('CACHE HIT: %s (%d valid station-days)\n',p,nnz(isfinite(X)));
            else
                fprintf('CACHE INVALIDATED: %s source, dates, or stations changed\n',p);
            end
        catch ME
            fprintf('CACHE unreadable for %s: %s\n',p,ME.message);
        end
    end
    if ~gotCache
        X=NaN(N,S);
        if strcmp(p,'MSWEP')
            X=loadMSWEPCached(C,files.(p),allDates,lat,lon,X,stationSignature);
        else
            for f=1:numel(files.(p))
                if mod(f,10)==1 || f==numel(files.(p))
                    fprintf('%s file %d/%d\n',p,f,numel(files.(p)));
                end
                X=insertNetcdf(files.(p){f},allDates,lat,lon,X,p,C);
            end
        end
        if C.useCache
            meta=struct('version',C.cacheVersion,'files',signature, ...
                'stations',stationSignature,'dates',[C.startDate C.endDate], ...
                'size',[N S],'ERA5mmAlready',C.ERA5mmAlready);
            saveCacheAtomically(cache,X,meta);
            fprintf('CACHE SAVED: %s (%d valid station-days)\n',p,nnz(isfinite(X)));
        end
    end
    D.(p)=X;
end
save(fullfile(C.outDir,'evaluation_config.mat'),'C','allDates','locInfo','-v7.3');

% #1 COMMON PERIOD: calendar 2010-2017, not "days where at least one grid is valid".
% First filter gauge independently of products, then require same station/day
% for all five products. This prevents different samples per product.
commonDay=true(N,1); % entire prescribed 2010-01-01..2017-12-31
stationGaugeAllowed=(sum(isfinite(obs),1)/N)>C.minCompleteness;
if ~any(stationGaugeAllowed),error('No gauges pass strict >70%% over 2010-2017');end
for j=1:numel(C.products)
    p=C.products{j}; z=D.(p);
    fprintf('%s: finite daily station-cells = %d, stations with ANY data = %d\n', ...
         p,nnz(isfinite(z)),nnz(any(isfinite(z),1)));
    if C.failOnIncompleteProduct && ~any(isfinite(z(:)))
        error('%s has ZERO valid extracted daily data. Check time/units/ZIP structure; NO comparisons produced.',p);
    end
end
commonMask=isfinite(obs);
for j=1:numel(C.products)
    commonMask=commonMask & isfinite(D.(C.products{j}));
end
pairedCoverage=sum(commonMask,1)/N;
stationAllowed=stationGaugeAllowed & pairedCoverage>C.minCompleteness;
if ~any(stationAllowed)
    z=locInfo;z.GaugeCoverage=fullCoverage;z.AllProductsPairedCoverage=pairedCoverage(:);
    writetable(z,fullfile(C.outDir,'diagnostic_station_coverage.csv'));
    error('No stations have >70%% identical six-source daily coverage. See diagnostic_station_coverage.csv. Do not change threshold silently.');
end
if sum(stationAllowed)<C.minPairedStations
    warning('Only %d stations pass; national summary will be fragile',sum(stationAllowed));
end
commonMask(:,~stationAllowed)=false;
fprintf('COMMON PERIOD: 2010-01-01 to 2017-12-31; %d days, %d qualifying stations.\n',N,sum(stationAllowed));
stationSummary=locInfo;
stationSummary.GaugeCoverage=fullCoverage;
stationSummary.AllProductsPairedCoverage=pairedCoverage(:);
stationSummary.Included=stationAllowed(:);
writetable(stationSummary,fullfile(C.outDir,'common_station_selection.csv'));
writetable(table(allDates,'VariableNames',{'CommonCalendarDates'}), ...
    fullfile(C.outDir,'common_calendar_dates.csv'));

% A single fixed mask for all products and all shifts, so shift comparisons
% cannot be explained by differences in included station-days.
shiftMask=commonMask;
for j=1:numel(C.products)
    xx=D.(C.products{j});
    neighbours=false(N,S);
    neighbours(2:end-1,:)=isfinite(xx(1:end-2,:)) & isfinite(xx(3:end,:));
    shiftMask=shiftMask & neighbours;
end
fprintf('Strict sensitivity days (all products +/-1): %d station-days\n',nnz(shiftMask));
if nnz(shiftMask)<100
    warning('Strict +/-1 common sample very small; interpret timing test with caution.');
end
scoreRows=cell(0,1); shiftRows=cell(0,1); annualRows=cell(0,1);
for j=1:numel(C.products)
    p=C.products{j}; X=D.(p);
    for k=find(stationAllowed)
        m=commonMask(:,k);
        v=skill(obs(m,k),X(m,k),C.wetThreshold);
        v.RX1dayBias=annualRX1Bias(obs(:,k),X(:,k),allDates,m);
        scoreRows{end+1,1}={p,char(stid(k)),sum(m),v.R,v.KGE,v.RMSE,v.Bias,v.POD,v.FAR,v.CSI,v.RX1dayBias}; %#ok<AGROW>
        % #5: shift -1, 0, +1 of product DATE LABEL, not station date.
        % Shift is assessed only over dates valid for all 3 shifts, for fairness.
        sh=NaN(N,3);
        sh(2:end,1)=X(1:end-1,k);  % shift -1 label: product from previous day
        sh(:,2)=X(:,k);
        sh(1:end-1,3)=X(2:end,k);  % shift +1 label: product from following day
        fixed=shiftMask(:,k); % globally matched station-days for all products/offsets
        for t=1:3
            vv=skill(obs(fixed,k),sh(fixed,t),C.wetThreshold);
            vv.RX1dayBias=annualRX1Bias(obs(:,k),sh(:,t),allDates,fixed);
            shiftRows{end+1,1}={p,char(stid(k)),t-2,sum(fixed),vv.R,vv.KGE,vv.RMSE,vv.POD,vv.FAR,vv.CSI,vv.RX1dayBias}; %#ok<AGROW>
        end
        % Annual mean rainfall on identical six-product paired days.
        yy=year(allDates);
        for y=unique(yy(commonDay))'
            a=m & yy==y;
            if sum(a)<ceil(C.minAnnualCompleteness*sum(commonDay & yy==y)), continue; end
            annualRows{end+1,1}={p,char(stid(k)),y,sum(a),mean(obs(a,k)),mean(X(a,k))}; %#ok<AGROW>
        end
    end
end
if isempty(scoreRows),error('No valid station scores calculated');end
scores=cell2table(vertcat(scoreRows{:}),'VariableNames', ...
    {'Product','StationID','Ndays','R','KGE','RMSE','Bias','POD','FAR','CSI','RX1dayBias'});
scores=normalizeTableCells(scores);
if isempty(shiftRows),error('No common valid days for +/-1 test');end
shifts=cell2table(vertcat(shiftRows{:}),'VariableNames', ...
    {'Product','StationID','ShiftDays','Ndays','R','KGE','RMSE','POD','FAR','CSI','RX1dayBias'});
shifts=normalizeTableCells(shifts);
if isempty(annualRows),error('No annual observations satisfy annual coverage');end
annual=cell2table(vertcat(annualRows{:}),'VariableNames', ...
    {'Product','StationID','Year','Ndays','GaugeAnnualDailyMean','ProductAnnualDailyMean'});
annual=normalizeTableCells(annual);
writetable(scores,fullfile(C.outDir,'common_station_scores.csv'));
writetable(shifts,fullfile(C.outDir,'day_shift_station_scores.csv'));
writetable(annual,fullfile(C.outDir,'annual_station_values.csv'));

% Product summary: median across stations; NOT thousands of days as independent replicates.
summaryRows=cell(numel(C.products),1);
for j=1:numel(C.products)
    p=C.products{j}; z=scores(strcmp(scores.Product,p),:);
    summaryRows{j}={p,height(z),median(z.R,'omitnan'),median(z.KGE,'omitnan'), ...
        median(z.RMSE,'omitnan'),median(z.Bias,'omitnan'),median(z.POD,'omitnan'), ...
        median(z.FAR,'omitnan'),median(z.CSI,'omitnan')};
end
summaryT=cell2table(vertcat(summaryRows{:}),'VariableNames', ...
    {'Product','Nstations','MedianR','MedianKGE','MedianRMSE','MedianBias','MedianPOD','MedianFAR','MedianCSI'});
summaryT=normalizeTableCells(summaryT);
writetable(summaryT,fullfile(C.outDir,'common_product_summary.csv'));

% #4: paired year-block bootstrap comparison of YEARLY ANOMALY correlations.
% Each annual pair uses the same station/year intersection as all products.
years=unique(year(allDates(commonDay))); Y=numel(years); P=numel(C.products);
annualGauge=NaN(Y,1); annualProduct=NaN(Y,P);
for q=1:Y
    yy=(year(allDates)==years(q));
    % Equal station coverage: commonMask guarantees paired day values.
    m=commonMask & repmat(yy,N,1);
    numValid=sum(m,1);
    good=stationAllowed & numValid>=ceil(0.7*sum(commonDay & yy));
    if ~any(good), continue; end
    O=obs(:,good); mm=m(:,good);
    annualGauge(q)=mean(O(mm),'omitnan');
    for j=1:P
        Z=D.(C.products{j}); Z=Z(:,good);
        annualProduct(q,j)=mean(Z(mm),'omitnan');
    end
end
validYears=isfinite(annualGauge) & all(isfinite(annualProduct),2);
years=years(validYears); annualGauge=annualGauge(validYears); annualProduct=annualProduct(validYears,:);
writetable(array2table([years annualGauge annualProduct], ...
    'VariableNames',[{'Year','Gauge'},C.products]),fullfile(C.outDir,'annual_national_paired.csv'));
if numel(years)<5, error('Only %d complete annual observations: cannot estimate reliable annual correlations.',numel(years));end
if numel(years)<8,warning('Only %d annual observations, correlation differences have wide uncertainty.',numel(years));end
corrOut=NaN(P,1); for j=1:P, corrOut(j)=pearsonLocal(annualGauge,annualProduct(:,j)); end
writetable(table(string(C.products)',corrOut,'VariableNames',{'Product','AnnualCorrelation'}), ...
    fullfile(C.outDir,'annual_correlations.csv'));

% Bootstrap all products using the SAME sampled years in each replicate.
B=C.nBootstrap; bootR=NaN(B,P); yCount=numel(years);
for b=1:B
    idx=randi(yCount,[yCount,1]);
    for j=1:P
        bootR(b,j)=pearsonLocal(annualGauge(idx),annualProduct(idx,j));
    end
end
cmpRows=cell(0,1);
for j=1:P-1
    for k=j+1:P
        delta=corrOut(j)-corrOut(k);
        draws=bootR(:,j)-bootR(:,k); draws=draws(isfinite(draws));
        if numel(draws)<100, ci=[NaN NaN]; pv=NaN;
        else
            ci=[quantileLocal(draws,0.025),quantileLocal(draws,0.975)];
            pv=NaN; % percentile bootstrap is CI estimation, not a valid null p test
        end
        cmpRows{end+1,1}={C.products{j},C.products{k},delta,ci(1),ci(2),pv, ...
            ci(1)>0 || ci(2)<0,yCount}; %#ok<AGROW>
    end
end
cmp=cell2table(vertcat(cmpRows{:}),'VariableNames', ...
    {'ProductA','ProductB','DifferenceAnnualR','CI95_Low','CI95_High','BootstrapP_NotEstimated','CIExcludesZero','Nyears'});
cmp=normalizeTableCells(cmp);
writetable(cmp,fullfile(C.outDir,'paired_annual_correlation_bootstrap.csv'));

% Supplement #4: paired comparisons of daily skill (KGE and RMSE),
% first calculated independently within each station-year, then aggregated
% by YEAR (median station skill). The bootstrap resamples whole years.
% This does not treat individual daily values as independent replicates.
Yall=unique(year(allDates(commonDay))); K=NaN(numel(Yall),P); E=K;
for iy=1:numel(Yall)
    ym=year(allDates)==Yall(iy);
    for j=1:P
        tmpK=NaN(S,1);tmpE=NaN(S,1);
        X=D.(C.products{j});
        for s=find(stationAllowed)
            mask=commonMask(:,s) & ym;
            if sum(mask)<max(15,ceil(C.minAnnualCompleteness*sum(commonDay & ym))),continue;end
            v=skill(obs(mask,s),X(mask,s),C.wetThreshold);
            tmpK(s)=v.KGE;tmpE(s)=v.RMSE;
        end
        K(iy,j)=median(tmpK,'omitnan');
        E(iy,j)=median(tmpE,'omitnan');
    end
end
writetable(array2table([Yall K E], ...
    'VariableNames',[{'Year'},strcat(C.products,'_KGE'),strcat(C.products,'_RMSE')]), ...
    fullfile(C.outDir,'annual_median_station_skills.csv'));
metricNames={'KGE','RMSE'}; matrices={K,E}; cmpSkillRows=cell(0,1);
for mtype=1:numel(matrices)
    Z=matrices{mtype};
    for j=1:P-1
        for k=j+1:P
            paired=isfinite(Z(:,j)) & isfinite(Z(:,k));
            zj=Z(paired,j);zk=Z(paired,k);n=numel(zj);
            if n<5,continue;end
            observed=mean(zj-zk);
            draws=NaN(B,1);
            for b=1:B
                ii=randi(n,n,1);draws(b)=mean(zj(ii)-zk(ii));
            end
            lower=quantileLocal(draws,0.025);upper=quantileLocal(draws,0.975);
            pApprox=NaN; % CI inference only; not a calibrated hypothesis-test p-value
            cmpSkillRows{end+1,1}={metricNames{mtype},C.products{j},C.products{k}, ...
                observed,lower,upper,pApprox,lower>0||upper<0,n}; %#ok<AGROW>
        end
    end
end
if ~isempty(cmpSkillRows)
    compSkill=cell2table(vertcat(cmpSkillRows{:}),'VariableNames', ...
        {'Metric','ProductA','ProductB','DifferenceAminusB','CI95_Low','CI95_High', ...
        'BootstrapP_NotEstimated','CIExcludesZero','Nyears'});
compSkill=normalizeTableCells(compSkill);
    writetable(compSkill,fullfile(C.outDir,'paired_skill_year_bootstrap.csv'));
end

% Supplementary plot of shift-dependent median KGE and annual correlations.
fig=figure('Color','w','Position',[100 100 1000 570]); hold on;
for j=1:P
    p=C.products{j}; mm=NaN(3,1);
    for ii=1:3
        z=shifts(strcmp(shifts.Product,p) & shifts.ShiftDays==ii-2,:);
        mm(ii)=median(z.KGE,'omitnan');
    end
    plot(-1:1,mm,'-o','LineWidth',1.8,'MarkerSize',7);
end
set(gca,'FontSize',15,'XTick',-1:1); xlabel('Product-day offset (days)'); ylabel('Median station KGE');
legend(C.products,'Location','best'); grid on; title('Daily timestamp sensitivity (common period)');
exportgraphics(fig,fullfile(C.outDir,'supp_daily_shift_KGE.png'),'Resolution',300); close(fig);
fig=figure('Color','w','Position',[100 100 1000 570]);
bar(corrOut); set(gca,'XTick',1:P,'XTickLabel',C.products,'FontSize',14);
ylabel('Correlation of paired annual national mean rainfall'); ylim([-1 1]); grid on;
exportgraphics(fig,fullfile(C.outDir,'supp_annual_correlation.png'),'Resolution',300); close(fig);

fprintf('\nCOMPLETED: %d stations, %d common dates, %d years for annual bootstrap.\n', ...
    sum(stationAllowed),sum(commonDay),numel(years));
fprintf('Results at %s\n',C.outDir); diary off;

%% LOCAL FUNCTIONS
function [stations,logs]=readStationsCached(C)
    files=stationFileList(C.stationDir);
    sig=inputSignature(files);
    dest=fullfile(C.cacheDir,'cache_station_import.mat');
    if C.useCache && ~C.forceRebuildStations && isfile(dest)
        try
            Z=load(dest,'stations','logs','signature','version');
            if isequal(Z.signature,sig) && strcmp(Z.version,C.cacheVersion)
                stations=Z.stations;logs=Z.logs;
                fprintf('CACHE HIT: station Excel/CSV import (%d stations)\n',numel(stations));
                return;
            end
        catch ME
            fprintf('Station cache unreadable: %s\n',ME.message);
        end
    end
    fprintf('Reading station spreadsheets (first run only)...\n');
    [stations,logs]=readStations(C);
    if C.useCache
        signature=sig;version=C.cacheVersion;
        temp=[dest '.writing.mat'];
        save(temp,'stations','logs','signature','version','-v7.3');
        movefile(temp,dest,'f');
        fprintf('CACHE SAVED: station import (%d stations)\n',numel(stations));
    end
end

function F=stationFileList(folder)
    f=[dir(fullfile(folder,'*.xlsx'));dir(fullfile(folder,'*.xls'));dir(fullfile(folder,'*.csv'))];
    f=f(~startsWith(string({f.name}),'~$'));
    F=sort(fullfile({f.folder},{f.name})');
end

function sig=inputSignature(F)
    % Sorted list of (full path, bytes, modification time) for invalidation.
    F=sort(F(:));
    sig=strings(numel(F),3);
    for i=1:numel(F)
        st=dir(F{i});
        if isempty(st)
            sig(i,:)=[string(F{i}) "MISSING" "MISSING"];
        else
            sig(i,:)=[string(F{i}),string(st(1).bytes),string(st(1).datenum)];
        end
    end
end

function saveCacheAtomically(destination,X,meta)
    % Do not replace a good cache with a partial one after an interruption.
    temp=[destination '.writing.mat'];
    save(temp,'X','meta','-v7.3');
    movefile(temp,destination,'f');
end

function X=loadMSWEPCached(C,F,allDates,lat,lon,X,stationSignature)
    N=numel(allDates);S=numel(lat);
    zipDir=fullfile(C.cacheDir,'MSWEP_BY_ARCHIVE');
    if ~isfolder(zipDir),mkdir(zipDir);end
    for j=1:numel(F)
        f=F{j};
        if endsWith(lower(f),'.zip')
            if ~zipContainsEvaluationYear(f,C.startDate,C.endDate)
                fprintf('MSWEP %d/%d: skip outside period\n',j,numel(F));
                continue;
            end
        end
        [~,basename,ext]=fileparts(f);
        % Unique stable key based on full filename, without third-party hash.
        idx=sprintf('%03d',j);
        chunkCache=fullfile(zipDir,['mswep_' idx '_' basename '.mat']);
        sig=inputSignature({f});
        hit=false;
        if C.useCache && ~any(strcmp('MSWEP',C.forceRebuildProducts)) && isfile(chunkCache)
            try
                Z=load(chunkCache,'X','meta');
                hit=isfield(Z,'meta') && isequal(Z.meta.version,C.cacheVersion) && ...
                    isequal(Z.meta.files,sig) && ...
                    isequal(Z.meta.stations,stationSignature) && ...
                    isequal(Z.meta.dates,[C.startDate C.endDate]) && ...
                    isequal(size(Z.X),[N S]);
                if hit
                    if any(isfinite(Z.X(:)))
                        v=isfinite(Z.X); X(v)=Z.X(v);
                    end
                    fprintf('MSWEP %d/%d: CACHE HIT (%d station-days)\n',j,numel(F),nnz(isfinite(Z.X)));
                end
            catch ME
                fprintf('MSWEP per-archive cache read error: %s\n',ME.message);
            end
        end
        if hit,continue;end
        fprintf('MSWEP %d/%d: processing %s\n',j,numel(F),f);
        Xchunk=NaN(N,S);
        if endsWith(lower(f),'.zip')
            temp=fullfile(C.outDir,'TEMP_MSWEP_UNZIP');
            if ~isfolder(temp),mkdir(temp);end
            % Extracting an archive can be expensive, but only the FIRST time.
            % Clear abandoned files in temporary folder from earlier crashes.
            leftover=dir(fullfile(temp,'**','*.nc'));
            if ~isempty(leftover)
                fprintf('  Found old temporary NetCDF files; clearing temporary directory...\n');
                rmdir(temp,'s');mkdir(temp);
            end
            try
                unpacked=unzip(f,temp);
                nNC=0;
                for i=1:numel(unpacked)
                    if ~(endsWith(lower(string(unpacked{i})),'.nc') || ...
                         endsWith(lower(string(unpacked{i})),'.nc4')),continue;end
                    d=filenameDates(unpacked{i},'MSWEP',0);
                    if isempty(d) || d<C.startDate || d>C.endDate,continue;end
                    nNC=nNC+1;
                    Xchunk=insertNetcdf(unpacked{i},allDates,lat,lon,Xchunk,'MSWEP',C);
                end
                if nNC==0,error('ZIP has no date-matching .nc files: %s',f);end
                fprintf('  MSWEP %d daily grids imported from archive\n',nNC);
            catch ME
                error('MSWEP ZIP %s failed: %s',f,ME.message);
            end
            if C.cleanUnzippedTemp && isfolder(temp),rmdir(temp,'s');end
        else
            Xchunk=insertNetcdf(f,allDates,lat,lon,Xchunk,'MSWEP',C);
        end
        if C.useCache
            meta=struct('version',C.cacheVersion,'files',sig, ...
                'stations',stationSignature,'dates',[C.startDate C.endDate]);
            saveCacheAtomically(chunkCache,Xchunk,meta);
            fprintf('  MSWEP checkpoint SAVED: %d station-days\n',nnz(isfinite(Xchunk)));
        end
        v=isfinite(Xchunk); X(v)=Xchunk(v);
    end
end

function [stations,logs]=readStations(C)
    ff=[dir(fullfile(C.stationDir,'*.xlsx'));dir(fullfile(C.stationDir,'*.xls'));dir(fullfile(C.stationDir,'*.csv'))];
    ff=ff(~startsWith(string({ff.name}),'~$'));
    stations=struct('id',{},'name',{},'lat',{},'lon',{},'dates',{},'rain',{});
    R=cell(0,1);
    for i=1:numel(ff)
        f=fullfile(ff(i).folder,ff(i).name);
        try
            if endsWith(lower(f),'.csv'), raw=readcell(f); else, raw=readcell(f,'Sheet',1); end
            [id,nm,la,lo]=stationMetadata(raw,ff(i).name);
            h=0;
            for r=1:min(35,size(raw,1))
                a=upper(strtrim(string(raw{r,1})));
                if a=="TANGGAL" || a=="DATE",h=r;break;end
            end
            if h==0, error('Header TANGGAL not found'); end
            nRows=size(raw,1)-h;
            dates=NaT(nRows,1); rain=NaN(nRows,1);
            nValid=0;
            for r=h+1:size(raw,1)
                dt=parseGaugeDateFast(raw{r,1});
                if isnat(dt),continue;end
                rr=tonum(raw{r,2});
                if rr<0 || rr>C.maxReasonableRain, rr=NaN; end
                nValid=nValid+1;
                dates(nValid)=dateshift(dt,'start','day');
                rain(nValid)=rr;
            end
            dates=dates(1:nValid);rain=rain(1:nValid);
            if isempty(dates), error('No daily observations after header'); end
            [ud,~,ic]=unique(dates); vr=NaN(numel(ud),1);
            for q=1:numel(ud), z=rain(ic==q);vr(q)=mean(z,'omitnan');end
            if ~isfinite(la)||~isfinite(lo), error('Latitude/longitude not read'); end
            s.id=id;s.name=nm;s.lat=la;s.lon=lo;s.dates=ud;s.rain=vr;
            stations(end+1)=s; %#ok<AGROW>
            status='OK';
        catch ME
            status=['SKIP: ' ME.message];
        end
        R{end+1,1}={ff(i).name,status}; %#ok<AGROW>
    end
    if isempty(R), logs=table(strings(0,1),strings(0,1),'VariableNames',{'File','Status'});
    else, logs=cell2table(vertcat(R{:}),'VariableNames',{'File','Status'});end
end
function [id,nm,lat,lon]=stationMetadata(raw,fname)
    [~,id]=fileparts(fname);nm=id;lat=NaN;lon=NaN;
    for k=1:min(12,size(raw,1))
        label=upper(replace(strtrim(string(raw{k,1}))," ",""));
        vv=raw{k,min(3,size(raw,2))};
        if contains(label,'WMO'),id=char(string(vv));end
        if contains(label,'NAMA'),nm=char(string(vv));end
        if contains(label,'LINTANG')||contains(label,'LATITUDE'),lat=tonum(vv);end
        if contains(label,'BUJUR')||contains(label,'LONGITUDE'),lon=tonum(vv);end
    end
end
function n=tonum(v)
    if isnumeric(v)&&isscalar(v),n=double(v);return;end
    s=strtrim(string(v));
    if ismissing(s)||s==""||strcmpi(s,'-')||strcmpi(s,'8888')||strcmpi(s,'9999')||strcmpi(s,'TT')
        n=NaN;return;
    end
    n=str2double(strrep(s,',','.'));
end
function d=parseGaugeDateFast(x)
    % Fast deterministic date parsing. Avoid trying 10 datetime formats for
    % every gauge row (potentially millions of rows across Excel files).
    d=NaT;
    if isdatetime(x), d=x; return; end
    if isnumeric(x) && isscalar(x) && isfinite(x)
        if x>20000 && x<100000
            d=datetime(x,'ConvertFrom','excel');
        elseif x>=19000101 && x<=21001231 && floor(x)==x
            y=floor(x/10000); m=floor(mod(x,10000)/100); dd=mod(x,100);
            if m>=1 && m<=12 && dd>=1 && dd<=31
                try,d=datetime(y,m,dd);catch,end
            end
        end
        return;
    end
    s=strtrim(char(string(x)));
    if isempty(s) || strcmpi(s,'NaN'),return;end
    % Common BMKG dates: DD/MM/YYYY, DD-MM-YYYY, YYYY-MM-DD.
    tok=regexp(s,'^(\d{1,4})[\/\-.](\d{1,2})[\/\-.](\d{2,4})','tokens','once');
    if ~isempty(tok)
        a=str2double(tok{1});b=str2double(tok{2});c=str2double(tok{3});
        if numel(tok{1})==4
            y=a;m=b;dd=c;
        else
            dd=a;m=b;y=c;
            if y<100,if y<50,y=y+2000;else,y=y+1900;end;end
        end
        if y>=1900 && y<=2100 && m>=1 && m<=12 && dd>=1 && dd<=31
            try
                d=datetime(y,m,dd);
                if year(d)~=y || month(d)~=m || day(d)~=dd,d=NaT;end
            catch,d=NaT;end
        end
        return;
    end
    if ~isempty(regexp(s,'^\d{8}$','once'))
        xnum=str2double(s);y=floor(xnum/10000);
        m=floor(mod(xnum,10000)/100);dd=mod(xnum,100);
        if y>=1900 && y<=2100 && m>=1 && m<=12 && dd>=1 && dd<=31
            try,d=datetime(y,m,dd);catch,end
        end
        return;
    end
    % Rare textual month names, e.g. 01-Jan-2010.
    if ~isempty(regexp(s,'[A-Za-z]','once'))
        try,d=datetime(s,'InputFormat','dd-MMM-yyyy','Locale','en_US');catch,end
    end
end
function F=findInputFiles(path,p)
    F={};
    if isfile(path),F={path};return;end
    if ~isfolder(path),warning('%s path missing: %s',p,path);return;end
    a=[dir(fullfile(path,'**','*.nc'));dir(fullfile(path,'**','*.nc4'));dir(fullfile(path,'**','*.zip'))];
    if ~strcmp(p,'MSWEP'),a=a(~endsWith(lower(string({a.name})),'.zip'));end
    a=a(~[a.isdir]);
    F=fullfile({a.folder},{a.name})';
    F=sort(F);
end
function inspectNetcdf(f)
    I=ncinfo(f);fprintf('\nPREFLIGHT %s\n',f);
    for k=1:numel(I.Variables)
        v=I.Variables(k);
        if isempty(v.Dimensions), dims=string.empty(1,0); else, dims=string({v.Dimensions.Name}); end
        fprintf('  %s [%s]\n',v.Name,strjoin(dims,', '));
    end
end
function X=loadMSWEP(C,F,allDates,lat,lon,X)
    if isempty(F),warning('No MSWEP files found');return;end
    for j=1:numel(F)
        f=F{j};
        fprintf('MSWEP archive/data %d/%d: %s\n',j,numel(F),f);
        if endsWith(lower(f),'.zip')
            if ~C.tryMSWEPzip,continue;end
            if ~zipContainsEvaluationYear(f,C.startDate,C.endDate)
                fprintf('  Skipping ZIP (no files from evaluation years): %s\n',f);
                continue;
            end
            temp=fullfile(C.outDir,'TEMP_MSWEP_UNZIP');
            if ~exist(temp,'dir'),mkdir(temp);end
            % Only one zip is unpacked at a time; requires free space.
            try
                unpacked=unzip(f,temp);
                nNC=0;
                for i=1:numel(unpacked)
                    if endsWith(lower(string(unpacked{i})),{'.nc','.nc4'})
                        d=filenameDates(unpacked{i},'MSWEP',0);
                        if isempty(d) || d<C.startDate || d>C.endDate,continue;end
                        nNC=nNC+1;
                        X=insertNetcdf(unpacked{i},allDates,lat,lon,X,'MSWEP',C);
                    end
                end
                if nNC==0, error('ZIP contains no .nc/.nc4 files. Need a dedicated GeoTIFF reader or original MSWEP NetCDF files: %s',f);end
            catch ME
                error('Cannot process MSWEP archive %s: %s',f,ME.message);
            end
            if C.cleanUnzippedTemp && isfolder(temp),rmdir(temp,'s');end
        else
            X=insertNetcdf(f,allDates,lat,lon,X,'MSWEP',C);
        end
    end
end
function X=insertNetcdf(f,allDates,stationLat,stationLon,X,p,C)
    try
        info=ncinfo(f);
        [v,latV,lonV,timeV]=detectVars(info,p);
        if isempty(v)||isempty(latV)||isempty(lonV),error('No precipitation/lat/lon variables');end
        gv=ncinfo(f,v);
        dn=string({gv.Dimensions.Name}); sz=double(gv.Size); nd=numel(sz);
        L=double(ncread(f,latV));G=double(ncread(f,lonV));
        if ~isvector(L)||~isvector(G),error('Curvilinear 2D lon/lat unsupported; prepare regular grid');end
        L=L(:);G=G(:);
        latDim=find(contains(lower(dn),"lat")|strcmpi(dn,latV),1);
        lonDim=find(contains(lower(dn),"lon")|strcmpi(dn,lonV),1);
        if isempty(latDim),latDim=find(sz==numel(L),1);end
        if isempty(lonDim),lonDim=find(sz==numel(G) & (1:nd)~=latDim,1);end
        if isempty(latDim)||isempty(lonDim),error('Cannot identify lat/lon dimensions');end
        tDim=find(contains(lower(dn),'time'),1);
        [timeDates,hasTime]=getTime(f,timeV,C);
        if isempty(tDim)
            other=setdiff(1:nd,[latDim lonDim]);
            if ~isempty(other),tDim=other(1);end
        end
        if ~isempty(tDim) && hasTime && numel(timeDates)~=sz(tDim)
            warning('Time length mismatch: %s',f);hasTime=false;
        end
        if ~hasTime
            timeDates=filenameDates(f,p,~isempty(tDim)*sz(max(1,nonemptyDim(tDim))));
        end
        if isempty(timeDates),error('Date cannot be inferred for %s',f);end
        if isempty(tDim),timeDates=timeDates(1);end
        [ok,ix]=ismember(dateshift(timeDates(:),'start','day'),allDates);
        goodT=find(ok);
        if isempty(goodT),return;end
        if any(G<0) && any(stationLon>180),stationLon=mod(stationLon+180,360)-180;end
        if all(G>=0),qLon=mod(stationLon,360);else,qLon=stationLon;end
        latI=zeros(numel(stationLat),1);lonI=latI;
        for s=1:numel(stationLat)
            [~,latI(s)]=min(abs(L-stationLat(s)));
            [~,lonI(s)]=min(abs(G-qLon(s)));
        end
        % NetCDF point/time hyperslabs prevent reading enormous grids.
        fill=extractAttrs(gv);
        for s=1:numel(stationLat)
            start=ones(1,nd);count=ones(1,nd);
            start(latDim)=latI(s);start(lonDim)=lonI(s);
            if ~isempty(tDim)
                % Contiguous date range read, indexing selected target dates.
                start(tDim)=min(goodT);count(tDim)=max(goodT)-min(goodT)+1;
            end
            a=double(ncread(f,v,start,count));a=a(:);
            if isempty(tDim),a=a(1); else,a=a(goodT-min(goodT)+1);end
            a=cleanGridValues(a,fill,p,gv);
            % Product daily accumulation must be in mm/day here.
            % ERA5 total precipitation is conventionally in metres; convert.
            if strcmp(p,'ERA5')
                units=lower(string(attrValue(gv,'units')));
                if units=="m" || units=="metres" || units=="meter"
                    a=a*1000;
                elseif units=="" && ~C.ERA5mmAlready
                    error('ERA5 missing units and ERA5mmAlready=false; inspect preprocessing');
                end
            end
            % Never silently convert rates to amounts without temporal context.
            X(ix(goodT),s)=a;
        end
    catch ME
        error('%s read failed (%s): %s',p,f,ME.message);
    end
end
function i=nonemptyDim(x)
    if isempty(x),i=1;else,i=x;end
end
function [v,latV,lonV,timeV]=detectVars(I,p)
    names=string({I.Variables.Name}); low=lower(names);
    latV=pick(names,low,{'lat','latitude','y'});
    lonV=pick(names,low,{'lon','longitude','x'});
    timeV=pick(names,low,{'time','valid_time','date'});
    switch p
        case 'ERA5', wanted={'tp','total_precipitation','precipitation','precip','pr'};
        case 'SAOBS', wanted={'rr','precip','precipitation','pr'};
        case 'GPCP',wanted={'precip','precipitation','pr','rain'};
        otherwise,wanted={'precipitation','precip','pr','rainfall','rr'};
    end
    v=pick(names,low,wanted);
end
function x=pick(names,low,options)
    x='';
    for z=1:numel(options)
        k=find(low==options{z},1);
        if ~isempty(k),x=char(names(k));return;end
    end
end
function [d,ok]=getTime(f,v,C)
    d=datetime.empty(0,1);ok=false;
    if ~C.preferNetcdfTime||isempty(v),return;end
    try
        V=ncinfo(f,v);t=double(ncread(f,v));units=string(attrValue(V,'units'));
        tok=regexp(char(units),'(?i)(seconds?|minutes?|hours?|days?)\s+since\s+(\d{4}-\d{1,2}-\d{1,2})(?:[ T](\d\d:\d\d:\d\d))?','tokens','once');
        if isempty(tok),return;end
        base=datetime(tok{2},'InputFormat','yyyy-MM-dd');
        if numel(tok)>=3&&~isempty(tok{3}),base=base+timeofday(datetime(tok{3},'InputFormat','HH:mm:ss'));end
        u=lower(tok{1});
        if startsWith(u,'sec'),d=base+seconds(t(:));
        elseif startsWith(u,'min'),d=base+minutes(t(:));
        elseif startsWith(u,'hour'),d=base+hours(t(:));
        else,d=base+days(t(:));end
        ok=true;
    catch,ok=false;end
end
function d=filenameDates(f,p,T)
    [~,nm]=fileparts(f);d=datetime.empty(0,1);
    if strcmp(p,'CHIRPS')
        t=regexp(nm,'(?:^|\.)(19\d\d|20\d\d)\.days','tokens','once');
        if ~isempty(t),y=str2double(t{1});d=(datetime(y,1,1):caldays(1):datetime(y,12,31))';end
    end
    if strcmp(p,'MSWEP')
        t=regexp(nm,'(?<!\d)((?:19|20)\d\d)(\d{3})(?!\d)','tokens','once');
        if ~isempty(t),d=datetime(str2double(t{1}),1,1)+caldays(str2double(t{2})-1);end
        if isempty(d)
            t=regexp(nm,'(?<!\d)((?:19|20)\d\d)[-_]?([01]\d)[-_]?([0-3]\d)(?!\d)','tokens','once');
            if ~isempty(t),d=datetime(str2double(t{1}),str2double(t{2}),str2double(t{3}));end
        end
    end
    if T>1&&numel(d)==1,d=d+caldays(0:T-1)';end
    if T>1&&numel(d)~=T,d=datetime.empty(0,1);end
end
function val=attrValue(V,name)
    val='';
    for k=1:numel(V.Attributes)
        if strcmpi(V.Attributes(k).Name,name),val=V.Attributes(k).Value;return;end
    end
end
function f=extractAttrs(V)
    f=[];
    for k=1:numel(V.Attributes)
        if strcmpi(V.Attributes(k).Name,'_FillValue')||strcmpi(V.Attributes(k).Name,'missing_value')
            f=[f double(V.Attributes(k).Value(:)')]; %#ok<AGROW>
        end
    end
end
function a=cleanGridValues(a,fill,p,V)
    for k=1:numel(fill),a(a==fill(k))=NaN;end
    a(abs(a)>1e15)=NaN;
    % Scale/offset are applied by ncread in modern MATLAB. Avoid double scaling.
    if strcmp(p,'ERA5')
        units=lower(string(attrValue(V,'units')));
        if contains(units,'kg')&&contains(units,'s-1')
            warning('ERA5 rain is a rate; aggregate to daily totals BEFORE this analysis.');
        end
    end
    a(a<0)=NaN;
end
function v=skill(o,x,wet)
    m=isfinite(o)&isfinite(x);o=o(m);x=x(m);
    v=struct('R',NaN,'KGE',NaN,'RMSE',NaN,'Bias',NaN,'POD',NaN,'FAR',NaN,'CSI',NaN,'RX1dayBias',NaN);
    if numel(o)<5,return;end
    r=pearsonLocal(o,x);v.R=r;
    sdo=std(o);sdx=std(x);mo=mean(o);mx=mean(x);
    if sdo>0&&sdx>0&&mo>0
        v.KGE=1-sqrt((r-1)^2+(sdx/sdo-1)^2+(mx/mo-1)^2);
    end
    v.RMSE=sqrt(mean((x-o).^2));v.Bias=mean(x-o);
    H=sum(o>=wet & x>=wet);M=sum(o>=wet & x<wet);F=sum(o<wet & x>=wet);
    if H+M>0,v.POD=H/(H+M);end
    if H+F>0,v.FAR=F/(H+F);end
    if H+M+F>0,v.CSI=H/(H+M+F);end
    v.RX1dayBias=max(x)-max(o);
end
function r=pearsonLocal(a,b)
    m=isfinite(a)&isfinite(b);a=a(m);b=b(m);r=NaN;
    if numel(a)<4,return;end
    a=a-mean(a);b=b-mean(b);
    den=sqrt(sum(a.^2)*sum(b.^2));if den>0,r=sum(a.*b)/den;end
end
function q=quantileLocal(x,prob)
    x=sort(x(isfinite(x)));if isempty(x),q=NaN;return;end
    i=1+(numel(x)-1)*prob;k=floor(i);h=i-k;
    if k>=numel(x),q=x(end);else,q=(1-h)*x(k)+h*x(k+1);end
end

function med=annualRX1Bias(o,x,dates,mask)
% Median of annual RX1day(product)-RX1day(gauge) with paired daily dates.
    y=year(dates);ys=unique(y(mask));d=NaN(numel(ys),1);
    for i=1:numel(ys)
        z=mask & y==ys(i) & isfinite(o) & isfinite(x);
        if sum(z)<30,continue;end
        d(i)=max(x(z))-max(o(z));
    end
    med=median(d,'omitnan');
end

function T=normalizeTableCells(T)
% Different MATLAB releases preserve numeric cell columns differently.
    for c=1:width(T)
        z=T.(c);
        if iscell(z) && all(cellfun(@(x) isnumeric(x)&&isscalar(x),z))
            T.(c)=cell2mat(z);
        end
    end
end


function tf=zipContainsEvaluationYear(z,d1,d2)
% Inspect ZIP entry filenames without extraction. A ZIP with no target-year
% NetCDF files is skipped; no assumption about ZIP's own download timestamp.
    tf=false;
    h=java.util.zip.ZipFile(java.io.File(z));
    c=onCleanup(@() h.close());
    entries=h.entries();
    while entries.hasMoreElements()
        entry=entries.nextElement(); nm=char(entry.getName());
        if ~(endsWith(lower(nm),'.nc') || endsWith(lower(nm),'.nc4')),continue;end
        d=filenameDates(nm,'MSWEP',0);
        if ~isempty(d) && d>=d1 && d<=d2,tf=true;return;end
    end
end

function showZipContents(z)
% Java ZIP listing: no decompression of multi-GB archive for preflight.
try
    zis=java.util.zip.ZipFile(java.io.File(z));
    e=zis.entries(); k=0; total=0;
    while e.hasMoreElements()
        a=e.nextElement();total=total+1;
        if k<20, fprintf('    %s\n',char(a.getName())); k=k+1;end
    end
    zis.close();
    fprintf('  Entries in ZIP: %d (first %d shown)\n',total,k);
catch ME
    warning('Could not list ZIP; inspect manually with 7-Zip: %s',ME.message);
end
end
