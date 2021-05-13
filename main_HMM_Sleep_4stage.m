clc;clear;close all;

% fullpath = 'C:\Users\Cheolsoo\Copy\Sleep_PHYSIO\PHYSIONET\';
% fullpath = 'C:\Users\ipsl\Copy\Sleep_PHYSIO\PHYSIONET\';
% fullpath = 'C:\Users\BMCL\Desktop\HMM_HyunSoo\from_Prof\HMM_USB\HMM Training-20190521T025158Z-001\HMM Training\HMMcodes_CP\DATA\';
fileLocation ='F:\Sleep_data\Biopac\EEG_500Hz\EEG_500Hz\';
X = dir('F:\Sleep_data\Biopac\EEG_500Hz\EEG_500Hz\*.mat');
% filename_FpzCz ={'sc4002e0_FpzCz.ascii','sc4012e0_FpzCz.ascii','sc4102e0_FpzCz.ascii','sc4112e0_FpzCz.ascii','st7022j0_FpzCz.ascii','st7052j0_FpzCz.ascii','st7121j0_FpzCz.ascii','st7132j0_FpzCz.ascii'}; %
% filename_PzOz ={'sc4002e0_PzOz.ascii','sc4012e0_PzOz.ascii','sc4102e0_PzOz.ascii','sc4112e0_PzOz.ascii','st7022j0_PzOz.ascii','st7052j0_PzOz.ascii','st7121j0_PzOz.ascii','st7132j0_PzOz.ascii'}; %
% filename_EOG ={'sc4002e0_EOG.ascii','sc4012e0_EOG.ascii','sc4102e0_EOG.ascii','sc4112e0_EOG.ascii','st7022j0_EOG.ascii','st7052j0_EOG.ascii','st7121j0_EOG.ascii','st7132j0_EOG.ascii'}; %
% filename_EMG ={'sc4002e0_EMG_1Hz.ascii','sc4012e0_EMG_1Hz.ascii','sc4102e0_EMG_1Hz.ascii','sc4112e0_EMG_1Hz.ascii','st7022j0_EMG_1Hz.ascii','st7052j0_EMG_1Hz.ascii','st7121j0_EMG_1Hz.ascii','st7132j0_EMG_1Hz.ascii'}; %
% filenameHyp={'sc4002e0hyp.ascii','sc4012e0hyp.ascii','sc4102e0hyp.ascii','sc4112e0hyp.ascii','st7022j0hyp.ascii','st7052j0hyp.ascii','st7121j0hyp.ascii','st7132j0hyp.ascii'};  %
% load([fullpath,'TimeOfInterest.mat']);
% bands = [12 40;12 16;8 12;2 7;0.5 4];      % Gamma; Beta; Sigma; Theta; Delta
% bands = [35 50; 20 30; 10.15 15.75; 1 3]; %Gamma; Beta; Sigma; Delta
bands = [28.5 50; 20 28; 10.15 15.75; 1 3];

fs=500; % sampling frequency, in Hz
fs_EMG = 1;
% EM parameter setting
random_initial = 0;     % choose random initial parameters 1: yes, 0: no
all_combination = 1;    % try all different combinations of output X (1)
single_combination = 0; % try only a matrix

Wt = 100*[0    0.5000    0.7500    1.0000;
    -0.5000    0.2500    0.3750    0.5000;
    -0.7500    -0.3750    0.5625    0.7500;
    -1.0000    -0.5000    -0.7500    1.0000];
Wt2 = 100*[0    0.5000    0.7500    1.0000;
    -0.5000    0.2500    0.3750    0.5000;
    -0.7500    0.3750    0.5625    0.7500;
    -1.0000    0.5000    0.7500    1.0000];



Tnn = [];
PER = [];
Precision = [];
n=0;

for sub = 1%:length(X)
    sub
    A = load([fileLocation X(sub).name]);
    fileName = X(sub).name;
%     if isfield(EEG1,'EEG1a')% EEG1구조체 안에 각 변수가 있는지 확인 for Gtec
%         xx = EEG1.EEG1a;
%     else
%         xx = EEG1.EEG1b;
%     end
     
     xx=A.EEG(:,1)';

%     xx = xx((3600-(str2num(fileName(end-8:end-7))*60+str2num(fileName(end-5:end-4))))*fs:end);
    
    x=((xx+32768)*0.00581368734264134-192); % Calibration in uV
    
    % The sleep stages W, 1, 2, 3, 4, R, M and 'unscored' are coded in the
    % file as binaries 0, 1, 2, 3, 4, 5, 6 and 9
    %     hyp_tmp = ones(size(hyp)).*9;
    %     hyp_tmp(find(hyp==0)) = 1;     % Wake
    %     hyp_tmp(find(hyp==5)) = 2;     % REM
    %     hyp_tmp(find(hyp==1)) = 3;     % Light
    %     hyp_tmp(find(hyp==2)) = 3;     % Light
    %     hyp_tmp(find(hyp==3)) = 4;     % Deep
    %     hyp_tmp(find(hyp==4)) = 4;     % Deep
    %     hyp_tmp = hyp;
    %     hypnogram_simple_tmp = hyp_tmp;
    
    EEG.srate = fs;
    EEG.data = x;
    EEG.subject = X(sub).name;
    
    resolution = .5; % in seconds (resolution: time resolution)
    frqlim = [1 fs/2]; % in Hz [1 250]
    nfreqs = 500; % number of freq bins, log spaced (nfreqs: frequency resolution)
    cycles = [3 round(frqlim(2)*.1)]; % [3 25] number of cycles at lowest and highest freqs (number of cycles at lowest freq defines the length of window of wavelet, ex) lowest freq = 1Hz and cycle(1) = 3 => window size of wavelet function is 3s, -1.5s~1.5s))
    xfac = round(nfreqs/1.5); % [667 freq bins] vertical smoothing factor across frequencies (larger=less smoothing)
    mvthresh = 9; % stds from mean of raw EEG
    yticks = [3,5,10,20,30,50,90,150,225]; % for plotting, labels on y axis
    tsmooth = 6; % sec (smooth across this number of seconds)
    chan = []; % default to 2 if > 1 chan, [] takes chan 1
    lfrqs = {[48 76]}; % line noise ranges in Hz, can be multiple or []
    dcfrqs = {[0 0.5]};
    
    
    pwr30 = ScoreMySleep_wavelet2(EEG,resolution,frqlim,nfreqs,cycles,xfac,yticks,mvthresh,tsmooth,chan,lfrqs,bands);
    B{sub} = pwr30;
end


Tn = [];
Precision = [];
for sub = 1%:length(X)
    fileName = X(sub).name;
   A= load([fileLocation X(sub).name]);
%     if isfield(EEG1,'EEG1a')% EEG1구조체 안에 각 변수가 있는지 확인 for Gtec
%         xx = EEG1.EEG1a;
%     else
%         xx = EEG1.EEG1b;
%     end
     xx=A.EEG(:,1)';
   
    x=((xx+32768)*0.00581368734264134-192)'; % Calibration in uV
   
%          fileLocation_act = 'F:\Sleep_data\Actigraphy\60sec_timestamp_with_Label\';
%      Fol = dir(['F:\Sleep_data\Actigraphy\60sec_timestamp_with_Label\' fileName(1:4) ' ' fileName(11:13) '*.csv']);
%      airlinedata = readtable([fileLocation_act Fol(1).name]);
%      file_Table = airlinedata(:,{'Var1'});
%       Act = table2array(file_Table(5:end,1))';
%       Act_ex = [];%Actgraphy's epoch : 60sec
%       for ii = 1:length(Act)
%           Act_ex = [Act_ex Act(ii) Act(ii)];
%       end
%       Act_ex = string(Act_ex);

    % The sleep stages W, 1, 2, 3, 4, R, M and 'unscored' are coded in the
    % file as binaries 0, 1, 2, 3, 4, 5, 6 and 9
    %     hyp_tmp = ones(size(hyp)).*9;
    %     hyp_tmp(find(hyp==0)) = 1;     % Wake
    %     hyp_tmp(find(hyp==5)) = 2;     % REM
    %     hyp_tmp(find(hyp==1)) = 3;     % Light
    %     hyp_tmp(find(hyp==2)) = 3;     % Light
    %     hyp_tmp(find(hyp==3)) = 4;     % Deep
    %     hyp_tmp(find(hyp==4)) = 4;     % Deep
    %     hyp_tmp = hyp;
    %     hypnogram_simple_tmp = hyp_tmp;
    
    EEG.srate = fs;
    EEG.data = x;
    EEG.subject = X(sub).name;
    
    pwr30 = B{sub};
    % EEG Fpz-Cz
    pwr30_FpzCz = pwr30(:,:);
    %     hypnogram_simple = hypnogram_simple_tmp(:);
    npwr30_FpzCz = reshape(zscore(pwr30_FpzCz(:)),size(pwr30_FpzCz,1),size(pwr30_FpzCz,2));
    clear pwr30_FpzCz pwr30;% hypnogram_simple_tmp;
    
    observations = [npwr30_FpzCz];
    %     observations = [pwr30_FpzCz;pwr30_PzOz];
    std_obs = std(observations');
    for i = 1:size(observations,1)
        observations(i,find(abs(observations(i,:))>10*std_obs(i))) = -5*std_obs(i);
    end
    
    
    % TOI2 = 1:size(observations,2)-1
    %     clear npwr30_FpzCz npwr30_PzOz;
    M = 1;      % number of Gaussian mixture
    nLatentStates = 4;      % number of latent states, wake, REM, light, deep
    % index of sleep stage, 1:Wake, 2: REM, 3:Light, 4:Deep
    a = [1 2 3 4];
    
    b = observations;
    b = reshape(zscore(b(:)),size(b));
    LimY = 5*std(b');
    for iii=2:length(b)
        if (-1*LimY(1,1) > b(1,iii))||(LimY(1,1)<b(1,iii))
            b(:,iii) = b(:,iii-1);
        elseif (-1*LimY(1,2) > b(2,iii))||(LimY(1,2)<b(2,iii))
            b(:,iii) = b(:,iii-1);
        elseif (-1*LimY(1,3) > b(3,iii))||(LimY(1,3)<b(3,iii))
            b(:,iii) = b(:,iii-1);
        elseif (-1*LimY(1,4) > b(4,iii))||(LimY(1,4)<b(4,iii))
            b(:,iii) = b(:,iii-1);
        end
    end
    
    [nO,nTime] = size(observations);
    %% HMM Parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %initial guess of parameters
    if random_initial
        % initial state prob
        %     pi0 = zeros(1,nLatentStates);
        %     pi0(1) = 1;
        pi0 = rand(1,nLatentStates);
        pi0 = pi0./sum(pi0);
        
        % state transition matrix
        Q = rand(nLatentStates);
        for j = 1:nLatentStates
            Q(j,:) = normalise(Q(j,:));
        end
        
        
        % observation probability density
        for j = 1:nLatentStates
            for k = 1:M
                mu(j,:,k) = randn(nO,1);
                sigma(j,:,:,k) = eye(nO);
            end
        end
        R = R_matrix(b,mu,sigma);
    else
        % Theta (parameter) values = {pi,Q,R}
        % pi: Initialization Matrix P(X1=a1) where X takes the values a
        pi0 = [1 0 0 0];
        % Q: Transition Matrix P(Xi|X_i-1)
        % X = W   R  L   D
        Q = [.75 .01 .24 0;    % X = Wake
            .05 .88 .07 0;    % X = REM
            .18 .11 .55 .16;  % X = Light
            .02 0 .14 .84];   % X = Deep
        % Mu Matrix relates to R matrix = P(Yi|Xi)
        % Standard magnitude in all bands corresponding to different sleep stages
        % muMat: S X B (Sleep Stages by Frequency Bands)
        % X =     W  R  L  D
        muMat = [.7 .1 .1 .1;    % Y = Gamma [35 Hz]
            .1 .7 .1 .1;    % Y = Beta  [20-30 Hz]
            .1 .1 .7 .1;    % Y = Sigma [10.25-15.75 Hz]
            .1 .1 .1 .7];   % Y = Delta [1-3 Hz]
        
        sigma = [1 0 0 0;
            0 1 0 0;
            0 0 1 0;
            0 0 0 1];
        
        % Normalize GTM
        normMuMat = reshape(zscore(muMat(:)),size(muMat));
        % Building Mean and Variance matrices for M Gaussian mixture
        for j=1:nLatentStates
            for k=1:M
                mu(:,j,k)=normMuMat(:,j);
                sigma(:,:,j,k)=eye(nLatentStates);
            end
        end
        % Normalize Y (observation) matrix before finding R
        R = mixgauss_prob(b, mu, sigma); %nLatentState X TimeStep
        
    end
    
    % R matrix
    %
    Preci = 0;
    for nn = 1::30
        nn
        %% Test using fixed parameters
        
        % Test using MAP
        Xmap = MAPSEQ_Estimation(b,a,pi0,Q,R);
        
        % figure;
        % set(gca,'fontsize',14);
        % plot(Xmap,'b','linewidth',2);hold;
        % plot(hypnogram_simple','r--','linewidth',2);
        % legend('MAP','True');
        % title('MAP Sequence');
        % axis([0 length(hypnogram_simple) 0 5]);
        %% Calculate Accuracy
        %         hypnogram_ROI = hypnogram_simple;
        
        Xmap_ROI = Xmap;
        
        % Before EM
        %         PER_MAP = 100*length(find(hypnogram_ROI(find(hypnogram_ROI(1:end-1)~=9))'==Xmap_ROI(find(hypnogram_ROI(1:end-1)~=9))))./length(Xmap_ROI);
        
        %% EM algorithm
        iterN = 30;
        Rall = R;
        pi0update_tmp = zeros(1,nLatentStates);
        Qupdate_tmp = zeros(nLatentStates,nLatentStates);
        Rupdate_tmp.mu = zeros(nLatentStates,nO);
        Rupdate_tmp.sigma = zeros(nLatentStates,nO,nO);
        iterN_EM = 1;
        for i = 1:iterN_EM
            [pi0update_,Qupdate_,Rupdate_] = EM_CONT(pi0,b,Rall,R,Q,iterN);
            pi0update_tmp = pi0update_tmp +pi0update_;
            Qupdate_tmp = Qupdate_tmp+Qupdate_;
            Rupdate_tmp.mu = Rupdate_tmp.mu+Rupdate_.mu;
            Rupdate_tmp.sigma = Rupdate_tmp.sigma+Rupdate_.sigma;
            clear pi0update_ Qupdate_ Rupdate_;
        end
        pi0update = pi0update_tmp./iterN_EM;
        Qupdate = Qupdate_tmp./iterN_EM;
        Rupdate.mu = Rupdate_tmp.mu./iterN_EM;
        Rupdate.sigma = Rupdate_tmp.sigma./iterN_EM;
        
        [~,bb] = max(squeeze(Rupdate.mu(:,:,1)),[],2);
        [aaa,bbb]= hist(bb,1:4);
        
        if single_combination == 1
            v = [1 2 3 4];
            
        elseif all_combination || any(aaa>=3) || length(find(aaa==2))>=2   % if the maximum of band power is not well sorted
            
            v = perms(a);
            disp('=========== Try all combinations! =========== ');
        else
            [~,v] = sort(bb);
            v = v';
            disp(['=========== Try only ',num2str(v),' ===========' ]);
        end
        
        for vi = 1:size(v,1)
            k = 1;
            tmp = v(vi,:);
            for i = tmp
                mu(k,:) = Rupdate.mu(i,:);
                sigma(k,:,:) = squeeze(Rupdate.sigma(i,:,:));
                pi0(1,k) = pi0update(i);
                k2 = 1;
                for i2 = tmp
                    Q(k2,k) = Qupdate(i2,i);
                    k2 = k2+1;
                end
                
                k = k+1;
            end
            
            % R matrix update using updated parameters
            %         R = mixgauss_prob(b, permute(mu,[2 1 3]), squeeze(permute(sigma,[2 3 1 4])));
            R = R_matrix(b,mu,sigma);
            
            %% Test using updated parameters
            
            % Test using MAP
            Xmap_EM = MAPSEQ_Estimation(b,a,pi0,Q,R);
            
            %% Calculate Accuracy after updating parameters
            Xmap_EM_ROI = Xmap_EM;
%             PER_MAP_EM(vi) = 100*length(find(hypnogram_ROI(find(hypnogram_ROI(1:end-1)~=9))'==Xmap_EM_ROI(find(hypnogram_ROI(1:end-1)~=9))))./length(Xmap_EM_ROI);
            para.pi0{vi} = pi0;
            para.Q{vi} = Q;
            para.mu{vi} = mu;
            para.sigma{vi} = sigma;
%             clear Q R pi0;
            
        end
%         
%         [PER_MAX,bb] = max(PER_MAP_EM);
%         PER_MAX
%         best_para.pi0 = (para.pi0{bb});
%         best_para.Q = para.Q{bb};
%         best_para.mu = para.mu{bb};
%         best_para.sigma = para.sigma{bb};
%         
%         disp(['The best combination is ', num2str(v(bb,:))]);
%         disp('Q matrix is ');
%         Qbar = best_para.Q
%         disp('mu matrix is ');
%         best_para.mu;
%         
%         k = 1;
%         tmp = v(bb,:);
%         for i = tmp
%             mu(k,:) = best_para.mu(i,:);
%             sigma(k,:,:) = best_para.sigma(i,:,:);
%             pi0(k) = best_para.pi0(i);
%             k2 = 1;
%             for i2 = tmp
%                 Q(k2,k) = best_para.Q(i2,i);
%                 k2 = k2+1;
%             end
%             
%             k = k+1;
%         end
        
        % R matrix update using updated parameters
        % R = mixgauss_prob(b, permute(mu,[2 1 3]), squeeze(permute(sigma,[2 3 1 4])));
        R = R_matrix(b,mu,sigma);
        
        Xmap_EM = MAPSEQ_Estimation(b,a,pi0,Q,R);
    Qbar = Q;

        Siz = length(Xmap_EM);
        times = Siz*30*fs;

        %% Estimated scoring
        perclt = 100*(length(find(Xmap_EM==3))/length(Xmap_EM));
        percdp = 100*(length(find(Xmap_EM==4))/length(Xmap_EM));
        percrem = 100*(length(find(Xmap_EM==2))/length(Xmap_EM));
        percwake = 100*(length(find(Xmap_EM==1))/length(Xmap_EM));
        count_r = 0;
        for jjj = 1:length(Xmap_EM)-1
            if ((Xmap_EM(jjj+1)==1)&&(Xmap_EM(jjj)>1))
                count_r = count_r+1;
            end
        end
        
        SQ_a = ( percrem*0.5 + perclt*0.75 + percdp)*100/(percrem + perclt+ percdp);
        ZQ_a = ( times/fs/3600*(1-percwake/100) + (times/fs/3600*percrem/100*0.5 + times/fs/3600*percdp/100*1.5) - (times/fs/3600*percwake/100*0.5 + count_r/15))*8.5;
        
        PQ2_a = ( (Qbar(2,2)*0.5+Qbar(2,3)*0.75+Qbar(2,4))*0.5 + (Qbar(3,2)*0.5+Qbar(3,3)*0.75+Qbar(3,4))*0.75 + (Qbar(4,2)*0.5+Qbar(4,3)*0.75+Qbar(4,4)) )*100/(Qbar(2,2)+Qbar(2,3)+Qbar(2,4)+Qbar(3,2)+Qbar(3,3)+Qbar(3,4)+Qbar(4,2)+Qbar(4,3)+Qbar(4,4));
        PQ7_a = ( (Qbar(2,2)*0.5+Qbar(2,3)*0.75+Qbar(2,4))*0.5 + ((-1)*Qbar(3,2)*0.5+Qbar(3,3)*0.75+Qbar(3,4))*0.75 + ((-1)*Qbar(4,2)*0.5-Qbar(4,3)*0.75+Qbar(4,4)) )*100/(Qbar(2,2)+Qbar(2,3)+Qbar(2,4)+Qbar(3,2)+Qbar(3,3)+Qbar(3,4)+Qbar(4,2)+Qbar(4,3)+Qbar(4,4));
        TQ_a = sum(sum(Qbar.*Wt));
        TQ2_a = sum(sum(Qbar.*Wt2));
        
        %         Precc = 0;
        %         for cnt=1:Siz
        %             if Xmap_EM(cnt)==hypnogram_simple(cnt)
        %                 Precc = Precc+1;
        %
        %             end
        %         end
        %         Precc = Precc/Siz*100;
        Xmap_EM2 = Xmap_EM;
%         for jj = 1:length(Act_ex)
%             if( Act_ex(jj) == "W"); Xmap_EM2(jj) = 1;end
%             if( Act_ex(jj) == "S")&&(Xmap_EM2(jj) ==1); Xmap_EM2(jj) = 0;end
%         end
        
        perclt = 100*(length(find(Xmap_EM2==3))/length(Xmap_EM2));
        percdp = 100*(length(find(Xmap_EM2==4))/length(Xmap_EM2));
        percrem = 100*(length(find(Xmap_EM2==2))/length(Xmap_EM2));
        percwake = 100*(length(find(Xmap_EM2==1))/length(Xmap_EM2));
        percunknown = 100*(length(find(Xmap_EM2==0))/length(Xmap_EM2));
        count = 0;
        for jjj = 1:length(Xmap_EM2)-1
            if ((Xmap_EM2(jjj+1)==1)&&(Xmap_EM2(jjj)>1))
                count = count+1;
            end
        end
        

        %         end
    end

    
    close all;
    figure;
    set(gca,'fontsize',14);
    plot(Xmap_EM,'b','linewidth',2);hold;
%     plot(Xmap_EM2','r--','linewidth',2);
    legend('MAP','True');
    title('MAP Sequence');
    axis([1 Siz 0 5]);
    set(gca,'YTick',1:4);
    set(gca,'YTickLabel',{ 'WAKE','REM','LIGHT','DEEP'});
    set(gca,'YDir','reverse');
       saveas(gcf,['F:\Sleep_data\Hypnogram\after_band\EEG1_' fileName '.jpg'])

    %     [PER_MAX(sub) best_para] = physio_EM(observations,hypnogram_simple,random_initial, 1:size(observations,2),all_combination,single_combination);
    %     fileID = fopen('C:\Users\ipsl\Copy\Sleep_PHYSIO\PHYSIONET\RESULTS\resultffteog_band2_gamma.txt','w');
    %     fprintf(fileID,num2str(PER_MAX));
    %     fclose(fileID);
end





%     writetable(Tn, 'Qmatrix_public_sleep_DB3.xls');
save('SNU_experiment_final_after_band_500Hz_EEG1.mat', 'Tnn');

% save('TQ_PQ_stack.mat', 'stack2');




