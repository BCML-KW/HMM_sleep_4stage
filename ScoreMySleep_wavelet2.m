% calls in a continuous data EEGLAB dataset and scores each 30 sec stretch
% 
%
%
% resolution -- [time step in sec] number of ms to advance with each analysis window
% 

% function [pwr1,freqs,deeppeaks,spindlemax,rempeaks,mndpvsltpwr,perclt,percdp,percrem,percwake,mvmnt] = ScoreMySleep(filename,fullpath,resolution,frqlim,nfreqs,cycles,xfac,yticks,mvthresh,origtsmooth,chan,lfrqs);
function nhyppwrbands = ScoreMySleep_wavelet2(EEG,resolution,frqlim,nfreqs,cycles,xfac,yticks,mvthresh,origtsmooth,chan,lfrqs,bands)

if ~exist('lfrqs') % 'lfrqs' 함수가 없다 => ~exist('lfrqs') = 1 => 실행
  lfrqs{1} = []; % do not remove any noise
end;
 
adv = round(EEG.srate*resolution); % convert to number of frames advance

if ~exist('chan')
  if size(EEG.data,1) > 1
    chan = 2;
  else
    chan = 1; % in case chan of interest is not 1
  end;
elseif isempty(chan) % chan 배열이 비어있는지 확인 <chan 함수가 있는데 비어있는가>
  if size(EEG.data,1) > 1  % EEG.data의 1번째 차원(행)의 갯수가 1보다 큰가
    chan = 2;
  else
    chan = 1;
  end
end;

%  spectral variables:
if isempty(frqlim)
  frqlim = [1 250]; 
end;
if isempty(nfreqs)
  nfreqs = 500; 
end;
if isempty(cycles)
  cycles = [2 60];   % number of cycles at lowest and highest freqs 
end;
if isempty(xfac)
  xwidth = round(nfreqs/16); % for smoothing across freqs
else
  xwidth = round(nfreqs/xfac); 
end;
if isempty(yticks)
  yticks = [2,5,10,20,30,50,90,150];
end;
if isempty(mvthresh)
  mvthresh = 9;
end;
winsize=EEG.srate*cycles(1)+1;
preticks = num2cell(yticks);
for x=1:length(preticks)
  preticks{x} = num2str(preticks{x});
end;

%EEG.data = rmbase(EEG.data);

% step through data doing time/freq decomp:
tt = 1;clear allpwr alldat times % 5 second start/end buffer
for t = round(winsize/2)+2*EEG.srate : adv : size(EEG.data,2)-2*EEG.srate
  alldat(:,tt) = EEG.data(chan, t-round(winsize/2)+1 : t+floor(winsize/2) )'; 
  times(1,tt) = t/EEG.srate;tt=tt+1;
end;
alldat = alldat - repmat(mean(alldat,1),[size(alldat,1) 1]);% in case no high pass
% repmat :  mean(alldat,1) 행렬을 size(alldat,1)행 1열 사이즈에 맞춰서 복사
% tag windows that are x stds off
a=std(alldat,[],1); % 가중치 0 에 2번째차원(열)을 따라 표준편차 반환. 즉, 표준편차로 이루어진 열벡터
thresh = [mean(a) - mvthresh*std(a)  mean(a) + mvthresh*std(a)];
mvmnt = find(a < thresh(1) | a > thresh(2)); % find: 0이 아닌 값의 인덱스 출력 => 조건에 부합하는 인덱스 출력


contigtms = [];
for ep = 1:size(alldat,2) % alldat의 2번째 차원 벡터(열벡터) 크기
  if length(find(alldat(:,ep) == alldat(1,ep))) == size(alldat,1) % all values the same
    contigtms = [contigtms,ep,ep-1,ep+1]; % remove epochs before and after
  end;
end;
contigtms = unique(contigtms);

% find and save crazy time points
%a = mean(alldat,1);
%thresh = [mean(a) - mvthresh*std(a)  mean(a) + mvthresh*std(a)];
%mvmnt = find(a < thresh(1) | a > thresh(2));
sleepdata.raw = EEG.data;
EEG.data = []; clear a
% calculates time points to be in hours:------------------------------
times = (times/60)/60;  % in hours
tstep = times(2) - times(1); % time step in hours
tsmooth = origtsmooth;
tsmooth = (tsmooth/60)/60;
tsmooth = round(tsmooth/tstep); % how many samples to smooth across
bigsmooth = round(.0083/tstep); % long smooth for secondary preferred freq (30 sec)
% Do spectral decomposition:-----------------------------------
% Rey's wavelets transform method:----------------------
%freqs = linspace(frqlim(1), frqlim(2), nfreqs);  % linear freqs 
freqs = linspace(log(frqlim(1)), log(frqlim(end)), nfreqs);freqs = exp(freqs);% log
wavelets = computepsifamilyQodd(freqs,1/EEG.srate,cycles(1),cycles(2),winsize); %

pwr1 = wavelets * alldat; 
pwr1 = pwr1.*conj(pwr1);% pwr1 is now non-complex
pwr1(find(pwr1==0)) = 10^-30;
pwr1 = 10*log10(pwr1);% convert to log
tmppwr = pwr1; 
 
sleepdata.rawpower = pwr1;

%clear alldat wavelets

% remove line noise:
if ~isempty(lfrqs)
  if ~isempty(lfrqs{1})
    for r=1:length(lfrqs)
      linefr = find(freqs > lfrqs{r}(1) & freqs< lfrqs{r}(2)); % eliminate line noise freqs
      if linefr(end)+3<=size(tmppwr,1)
           surrfr = [linefr(1)-3:linefr(1)-1, linefr(end)+1:linefr(end)+3]; 
      elseif linefr(end)+3>size(tmppwr,1)
           surrfr = [linefr(1)-3:linefr(1)-1, size(tmppwr,1)-2:size(tmppwr,1)]; 
      end
      surrfr(find(ismember(surrfr,linefr))) = [];
      tmppwr(linefr,:) = repmat(mean(tmppwr(surrfr,:)),[length(linefr) 1]);
      %tmppwr(linefr(1):end,:) = 0; % for cases of extreme line noise
    end;
  end;
end;
% smooth across times and freqs:--------
if xwidth > 1
  fprintf('\nSmoothing across %s frequency bins (out of %s total)\n',int2str(xwidth),int2str(size(tmppwr,1)));
  [outdata,outx] = movav(tmppwr',[],xwidth,1); 
  freqs = freqs(round(outx));% smooth
else
  outdata = tmppwr';
end;

% big smooth for secondary preferred freq
%[bigsmoothed,bigx] = movav(outdata',[],bigsmooth,1); 
%bigtimes = times(round(bigx));

fprintf('\nSmoothing across %s time bins (%s seconds)\n',int2str(tsmooth),int2str(origtsmooth));
[outdata,outx] = movav(outdata',[],tsmooth,1); 
pwr1 = outdata;  
outx = floor(outx); %outx(length(times):end) = [];
times = times(outx);

% find movement times that no longer exist:-----
mvmnt = intersect(mvmnt,outx);
% intersect: 교집합
mvmnt(find(mvmnt>length(times))) = [];

contigtms = intersect(contigtms,outx);
contigtms(find(contigtms>length(times))) = [];

% for plotting movement frames:-------
mtimes = times(mvmnt); % find movement frames in hours
clear outdata outx 
x=find(pwr1==-Inf); % find when electrodes not attached
bp = pwr1; bp(x) = 0;
if ~isempty(contigtms) || ~isempty(mvmnt)
  notcontig = [1:size(bp,2)]; notcontig([contigtms,mvmnt]) = [];
  fprintf('\nSetting amp saturation periods to mean power for display\n');
  ncm = mean(bp(:,notcontig),2);
  for c = 1:length(contigtms)
    bp(:,contigtms(c)) = ncm;
  end;
  for c = 1:length(mvmnt)
    bp(:,mvmnt(c)) = ncm;
  end;
  ctimes = times([contigtms,mvmnt]); % find movement frames in hours
end;

 % remove baseline
pwr1 = pwr1-repmat(median(bp,2),[1 size(pwr1,2)]);
% median : 각 열의 중앙값 반환

if ~isempty(contigtms) % now get rid off the messy bits from pwr1
  pwr1(:,contigtms) = 0;
end;
%if ~isempty(mvmnt) % now get rid off the messy bits from pwr1
%  pwr1(:,mvmnt) = 0;
%end;


% change time in sec
timessec = times*60*60;
%% Calculate magnitude in each band
for ib = 1:size(bands,1)
    freqIndex{ib} = find(freqs>bands(ib,1) & freqs<bands(ib,2));
end
k = 1;
for it = 1:30:timessec(end)
    timeIndex = find(timessec>it & timessec<it+29);
    for ifq = 1:length(freqIndex)
        hyppwrbands(ifq,k) = mean(mean(pwr1(freqIndex{ifq},timeIndex)));
    end
    k = k+1;
end
hyppwrbands(:,end+1) = hyppwrbands(:,end);
nhyppwrbands = reshape(zscore(hyppwrbands(:)),size(hyppwrbands,1),size(hyppwrbands,2));
