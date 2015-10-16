% fmrib_pas_EEGActivations() - Remove pulse artifacts from BCG-related selected 
%     ICs by PROJIC. A choice of methods is available. This program is merely an 
%     adaptation of the original code fmrib_pas(), as part of the EEGLAB toolbox, 
%     within the FMRIB plug-in.
%     
% Usage:
%    >> data = fmrib_pas_EEGActivations(data, chs, QRSevents, method, npc, avg_window, fs, delay)
%
% Inputs:
%   data       - IC time-courses returned by a given ICA decomposition
%                 algorithm
%   chs        - ICs for which the artefact occurrencs will be corrected
%                 for
%   QRSevents  - vector containing the R peak annotations
%   method: 'obs' for artifact principal components.  You need to specify
%           'npc', which is the number of PC to use in fitting the artifacts.  If unsure, use 4.
%           'mean' for simple mean averaging.
%           'gmean' for Gaussian weighted mean.
%           'median' for median filter.
%   npc        - number of principal components (set n_win=[] when using the 'obs' method)
%   avg_window - number of windows (set npc=[] when using the 'mean', 
%                 'gmean' or 'median' correction methods).
%   fs         - sampling rate
%   delay_QRS  - time delay between the R peak annotations and the BCG
%
% Outputs:
%   data - IC time-courses corrected for the artefact occurrences using the
%          specified method
%
%   Author:  Rami K. Niazy (adapted by Rodolfo Abreu, ISR/IST, Universidade de Lisboa, 2015)

function data = fmrib_pas_EEGActivations(data, chs, QRSevents, method, npc, avg_window, fs, delay)

nargchk(3,4,nargin);

data = data(chs, :);

switch method
    
    case 'obs'
        
        if nargin < 4
            error('Incorrect number of input arguments');
        end
        %init
        %-----
        [channels samples]=size(data);
        pcs=npc-1;
        
        %standard delay between QRS peak and artifact (allen,1998)
%         delay=round(0.21*fs);
        
        Gwindow=20;
        GHW=floor(Gwindow/2);
        rcount=0;
        firstplot=1;
        
        % memory allocation
        %------------------
        avgdata=zeros(Gwindow,round(fs*2));
        drift=zeros(1,round(fs*2)*Gwindow);
        fitted_art=zeros(channels,samples);
        mPA=zeros(1,round(2*fs));
        peakplot=zeros(1,samples);
        mEye=eye(channels);
        
        %Extract QRS events
        %------------------
        
        peakplot(QRSevents)=1;
        sh=zeros(1,delay);
        np=length(peakplot);
        peakplot=[sh peakplot(1:np-delay)];
        peakI=find(peakplot>0);
        peakcount=length(peakI);
        pa=peakcount;
        
        %make filter
        %-----------
        a=[0 0 1 1];
        f=[0 0.4/(fs/2) 0.9/(fs/2) 1];
        ord=round(3*fs/0.5);
        fwts=firls(ord,f,a);
        
        
        
        %Artifact Subtraction
        %---------------------
        
        % for the first window/2 points use arthemitic mean for averageing.
        % findg mean QRS peak-to-peak (R-to-R) interval
        for ch=1:channels
%             if ch==1
%                 %init waitbar
%                 %------------
%                 barth=5;
%                 barth_step=barth;
%                 Flag25=0;
%                 Flag50=0;
%                 Flag75=0;
%                 fprintf('\nPulse artifact subtraction in progress...Please wait!\n');
%             end
            
            RR=diff(peakI);
            mRR=median(RR);
            sRR=std(RR);
            %PArange=round(1.25*(mRR+sRR)/2); % half RR
            PArange=round(1.5*mRR/2); % half RR
            midP=PArange+1;
            
            %         if ch==1
            %             while ((QRSevents(pa-1)+2*PArange-1) > samples)
            %                 pa=pa-1;
            %             end
            %             steps=channels*pa;
            %         end
            
            if ch==1
                while (peakI(pa)+PArange > samples)
                    pa=pa-1;
                end
                steps=channels*pa;
                peakcount=pa;
            end
            
            eegchan=filtfilt(fwts,1,double(data(ch,:)));
            pcamat=zeros(pa-1,2*PArange+1);
            dpcamat=pcamat;
            for p=2:pa
                pcamat(p-1,:)=eegchan(peakI(p)-PArange:peakI(p)+PArange);
            end
            pcamat=detrend(pcamat','constant')';
            meaneffect=mean(pcamat);
            dpcamat=detrend(pcamat,'constant');
            [apc,ascore,asvar]=pca_calc(dpcamat');
            
            papc=[meaneffect' ascore(:,1:pcs)];
            
            
            try
                pad_fit=double(papc)*(double(papc)\...
                    double(detrend(data(ch,peakI(1)-PArange:...
                    peakI(1)+PArange)','constant')));
                
                fitted_art(ch,peakI(1)-PArange:peakI(1)+PArange)=...
                    pad_fit';
            catch
            end
            
            
            
            
            
            %-----------------------------------------------------
            for p=2:GHW+1
                
                pad_fit=double(papc)*(double(papc)\...
                    double(detrend(data(ch,peakI(p)-PArange:...
                    peakI(p)+PArange)','constant')));
                
                fitted_art(ch,peakI(p)-PArange:peakI(p)+PArange)=...
                    pad_fit';
                
                %update bar
                %----------
%                 percentdone=floor(((ch-1)*pa+p)*100/steps);
%                 if floor(percentdone)>=barth
%                     if percentdone>=25 & Flag25==0
%                         fprintf('25%% ')
%                         Flag25=1;
%                     elseif percentdone>=50 & Flag50==0
%                         fprintf('50%% ')
%                         Flag50=1;
%                     elseif percentdone>=75 & Flag75==0
%                         fprintf('75%% ')
%                         Flag75=1;
%                     elseif percentdone==100
%                         fprintf('100%%\n')
%                     else
%                         fprintf('.')
%                     end
%                     
%                     while barth<=percentdone
%                         barth=barth+barth_step;
%                     end
%                     if barth>100
%                         barth=100;
%                     end
%                 end
            end
            
            %---------------- Processing of central data ---------------------
            %cycle through peak artifacts identified by peakplot
            rcount=GHW;
            for p=GHW+2:peakcount-GHW-1
                PreP=ceil((peakI(p)-peakI(p-1))/2);
                PostP=ceil((peakI(p+1)-peakI(p))/2);
                if PreP > PArange
                    PreP=PArange;
                end
                if PostP > PArange
                    PostP=PArange;
                end
                
                pad_fit=double(papc(midP-PArange:midP+PArange,:))*...
                    (double(papc(midP-PArange:midP+PArange,:))\...
                    double(detrend(data(ch,peakI(p)-PArange:...
                    peakI(p)+PArange)','constant')));
                
                fitted_art(ch,peakI(p)-PreP:peakI(p)+PostP)=...
                    pad_fit(midP-PreP:midP+PostP)';
                
                
                
                %update bar
                %----------
%                 percentdone=floor(((ch-1)*pa+p)*100/steps);
%                 if floor(percentdone)>=barth
%                     if percentdone>=25 & Flag25==0
%                         fprintf('25%% ')
%                         Flag25=1;
%                     elseif percentdone>=50 & Flag50==0
%                         fprintf('50%% ')
%                         Flag50=1;
%                     elseif percentdone>=75 & Flag75==0
%                         fprintf('75%% ')
%                         Flag75=1;
%                     elseif percentdone==100
%                         fprintf('100%%\n')
%                     else
%                         fprintf('.')
%                     end
%                     
%                     while barth<=percentdone
%                         barth=barth+barth_step;
%                     end
%                     if barth>100
%                         barth=100;
%                     end
%                 end
            end
            
            %---------------- Processing of last GHW peaks------------------
            sectionPoints=samples-(peakI(peakcount-GHW)-PreP)+1;
            
            
            for p=peakcount-GHW:peakcount
                try
                    
                    pad_fit=double(papc(midP-PArange:midP+PArange,:))*...
                        (double(papc(midP-PArange:midP+PArange,:))\...
                        double(detrend(data(ch,peakI(p)-PArange:...
                        peakI(p)+PArange)','constant')));
                    
                    fitted_art(ch,peakI(p)-PreP:peakI(p)+PostP)=...
                        pad_fit(midP-PreP:midP+PostP)';
                    
                    
                catch
                end
                
                
                %update bar
                %----------
%                 percentdone=floor(((ch-1)*pa+p)*100/steps);
%                 if floor(percentdone)>=barth
%                     if percentdone>=25 & Flag25==0
%                         fprintf('25%% ')
%                         Flag25=1;
%                     elseif percentdone>=50 & Flag50==0
%                         fprintf('50%% ')
%                         Flag50=1;
%                     elseif percentdone>=75 & Flag75==0
%                         fprintf('75%% ')
%                         Flag75=1;
%                     elseif percentdone==100
%                         fprintf('100%%\n')
%                     else
%                         fprintf('.')
%                     end
%                     
%                     while barth<=percentdone
%                         barth=barth+barth_step;
%                     end
%                     if barth>100
%                         barth=100;
%                     end
%                 end
            end
            
        end
        
        data=data-fitted_art;
        
        
    otherwise
        
        %init
        %-----
        [channels samples]=size(data);
        
        %standard delay between QRS peak and artifact (allen,1998)
%         delay=round(0.21*fs);
        
        if isempty(avg_window), Gwindow=20; else Gwindow = avg_window; end;
        GHW=floor(Gwindow/2);
        rcount=0;
        firstplot=1;
        
        % memory allocation
        %------------------
        avgdata=zeros(channels,round(fs*2),Gwindow);
        drift=zeros(channels,round(fs*2)*Gwindow);
        clean=zeros(channels,samples);
        mPA=zeros(channels,round(2*fs));
        mMag=zeros(channels,Gwindow);
        peakplot=zeros(1,samples);
        mEye=eye(channels);
        
        %Extract QRS events
        %------------------
        
        peakplot(QRSevents)=1;
        sh=zeros(1,delay);
        np=length(peakplot);
        peakplot=[sh peakplot(1:np-delay)];
        peakI=find(peakplot>0);
        peakcount=length(peakI);
        
        %Artifact Subtraction
        %---------------------
        
        % for the first window/2 points use arthemitic mean for averageing.
        % findg mean QRS peak-to-peak (R-to-R) interval
        for p=2:GHW+1
            RR(p-1)=peakI(p)-peakI(p-1);
        end
        mRR=mean(RR);
        sRR=std(RR);
        PArange=round((mRR+sRR)/2); % half RR
        
        t=(0:peakI(GHW+1)+PArange-1)/fs;
        
        % subtraction of low freq trend; could use detrend.m function instead
        for n=1:channels
            pdrift=polyfit(t,data(n,1:peakI(GHW+1)+PArange),1);
            drift(n,1:peakI(GHW+1)+PArange)=polyval(pdrift,t);
        end
        
        data(:,1:peakI(GHW+1)+PArange)=...
            data(:,1:peakI(GHW+1)+PArange)-...
            drift(:,1:peakI(GHW+1)+PArange);
        
        
        for p=2:GHW+1
            avgdata(:,1:2*PArange+1,p-1)=...
                data(:,peakI(p)-PArange:peakI(p)+PArange);
        end
        
        for chan =1:channels
            avgdata(chan,:,:)=detrend(squeeze(avgdata(chan,:,:)),'constant');
        end
        
        mPA(:,1:2*PArange+1)=mean(avgdata(:,1:1+2*PArange,1:GHW),3);
        
        if peakI(1) > PArange
            
            alphaV=...
                sum(detrend(data(:,peakI(1)-PArange:peakI(1)+PArange)','constant')'...
                .*mPA(:,1:2*PArange+1),2)./...
                sum(mPA(:,1:2*PArange+1).*mPA(:,1:2*PArange+1),2);
            
            alphaD=diag(alphaV);
            
        else
            
            alphaV=sum(detrend(data(:,1:peakI(1)+PArange)','constant')'...
                .*mPA(:,PArange-peakI(1)+2:2*PArange+1),2)./...
                sum(mPA(:,PArange-peakI(1)+2:2*PArange+1).*...
                mPA(:,PArange-peakI(1)+2:2*PArange+1),2);
            
            alphaD=diag(alphaV);
            
        end
        
        
        data(:,1:peakI(GHW+1)+PArange)=...
            data(:,1:peakI(GHW+1)+PArange)+...
            drift(:,1:peakI(GHW+1)+PArange);
        
        
        %init waitbar
        %------------
%         barth=5;
%         barth_step=barth;
%         Flag25=0;
%         Flag50=0;
%         Flag75=0;
%         fprintf('\nPulse artifact subtraction in progress...Please wait!\n');
        
        
        if peakI(1) > PArange
            
            
            clean(:,peakI(1)-PArange:peakI(1)+PArange)=...
                data(:,peakI(1)-PArange:peakI(1)+PArange)-alphaD*mPA(:,1:2*PArange+1);
            startpoints=peakI(1)-PArange-1;
            clean(:,1:startpoints)=data(:,1:startpoints);
        else
            
            clean(:,1:peakI(1)+PArange)=...
                data(:,1:peakI(1)+PArange)-alphaD*mPA(:,PArange-peakI(1)+2:2*PArange+1);
        end
        
        %update bar
        %----------
%         percentdone=floor(1*100/peakcount);
%         if floor(percentdone)>=barth
%             if percentdone>=25 & Flag25==0
%                 fprintf('25%% ')
%                 Flag25=1;
%             elseif percentdone>=50 & Flag50==0
%                 fprintf('50%% ')
%                 Flag50=1;
%             elseif percentdone>=75 & Flag75==0
%                 fprintf('75%% ')
%                 Flag75=1;
%             elseif percentdone==100
%                 fprintf('100%%\n')
%             else
%                 fprintf('.')
%             end
%             
%             while barth<=percentdone
%                 barth=barth+barth_step;
%             end
%             if barth>100
%                 barth=100;
%             end
%         end
        
        
        %-----------------------------------------------------
        for p=2:GHW+1
            
            alphaV=sum(data(:,peakI(p)-PArange:peakI(p)+PArange).*...
                mPA(:,1:2*PArange+1),2)./...
                sum(mPA(:,1:2*PArange+1).*mPA(:,1:2*PArange+1),2);
            
            alphaD=diag(alphaV);
            
            clean(:,peakI(p)-PArange:peakI(p)+PArange)=...
                data(:,peakI(p)-PArange:peakI(p)+PArange)-...
                alphaD*mPA(:,1:2*PArange+1);
            
            %update bar
            %----------
%             percentdone=floor(p*100/(peakcount-1));
%             if floor(percentdone)>=barth
%                 if percentdone>=25 & Flag25==0
%                     fprintf('25%% ')
%                     Flag25=1;
%                 elseif percentdone>=50 & Flag50==0
%                     fprintf('50%% ')
%                     Flag50=1;
%                 elseif percentdone>=75 & Flag75==0
%                     fprintf('75%% ')
%                     Flag75=1;
%                 elseif percentdone==100
%                     fprintf('100%%\n')
%                 else
%                     fprintf('.')
%                 end
%                 
%                 while barth<=percentdone
%                     barth=barth+barth_step;
%                 end
%                 if barth>100
%                     barth=100;
%                 end
%             end
        end
        
        %---------------- Processing of central data ---------------------
        %cycle through peak artifacts identified by peakplot
        for p=GHW+2:peakcount-GHW-1
            PreP=ceil((peakI(p)-peakI(p-1))/2); % interval before peak artifact
            PostP=ceil((peakI(p+1)-peakI(p))/2);% interval after peak artifact
            sectionPoints=peakI(p+GHW)+PostP-(peakI(p-GHW)-PreP-1);
            
            % subtract drift
            t=(0:sectionPoints-1)/fs;
            for n=1:channels
                pdrift=...
                    polyfit(t,data(n,peakI(p-GHW)-PreP:peakI(p+GHW)+PostP),1);
                drift(n,1:sectionPoints)=polyval(pdrift,t);
            end
            
            data(:,peakI(p-GHW)-PreP:peakI(p+GHW)+PostP)=...
                data(:,peakI(p-GHW)-PreP:peakI(p+GHW)+PostP)-...
                drift(:,1:sectionPoints);
            
            rr=1;
            
            for r=p-GHW:p+GHW
                mMag(:,rr)=mean(abs(data(:,peakI(r)-PreP:peakI(r)+PostP)),2);
                rr=rr+1;
            end
            
            minMag=min(mMag,[],2);
            
            %only count those EEG data with no motion or line artifact by rejecting
            %excessively large variations.  Use channel 1 as ref (Cz).
            oldavgdata=avgdata;
            oldrcount=rcount;
            rcount=0;
            for r=p-GHW:p+GHW
                if (mean(abs(data(1,peakI(r)-PreP:peakI(r)+PostP))) ...
                        < 4*minMag(1) ) ...
                        & (max(abs(data(1,peakI(r)-PreP:peakI(r)+PostP)))...
                        < 2* 175)
                    avgdata(:,1:PreP+PostP+1,rcount+1)=...
                        data(:,peakI(r)-PreP:peakI(r)+PostP);
                    rcount=rcount+1;
                end
            end
            
            for chan =1:channels
                avgdata(chan,:,:)=detrend(squeeze(avgdata(chan,:,:)),'constant');
            end
            
            %     rcount
            %if more than 8 good data sections found then do averaging, else use
            %previous.
            if rcount > 5
                switch method
                    case 'gmean'
                        if p==peakcount-GHW-1
                            mPA=mean(avgdata(:,:,1:rcount),3);
                        elseif rcount==(2*GHW+1)
                            G=gausswin(rcount);
                            G=G/sum(G);
                            for r=1:rcount
                                avgdata(:,:,r)=avgdata(:,:,r)*G(r);
                            end
                            mPA=sum(avgdata(:,:,1:rcount),3);
                        else
                            mPA=median(avgdata(:,:,1:rcount),3);
                        end
                    case 'mean'
                        mPA=mean(avgdata(:,:,1:rcount),3);
                    case 'median'
                        mPA=median(avgdata(:,:,1:rcount),3);
                end
            else
                switch method
                    case 'mean'
                        mPA=mean(oldavgdata(:,:,1:oldrcount),3);
                    case 'median'
                        mPA=median(oldavgdata(:,:,1:oldrcount),3);
                    case 'gmean'
                        mPA=mean(oldavgdata(:,:,1:oldrcount),3);
                end
            end
            
            alphaV=sum(detrend(data(:,peakI(p)-PreP:peakI(p)+PostP)','constant')'.*...
                mPA(:,1:PreP+PostP+1),2)./...
                sum(mPA(:,1:PreP+PostP+1).*mPA(:,1:PreP+PostP+1),2);
            
            alphaD=diag(alphaV);
            
            data(:,peakI(p-GHW)-PreP:peakI(p+GHW)+PostP)=...
                data(:,peakI(p-GHW)-PreP:peakI(p+GHW)+PostP)+...
                drift(:,1:sectionPoints);
            
            
            clean(:,peakI(p)-PreP:peakI(p)+PostP)=...
                data(:,peakI(p)-PreP:peakI(p)+PostP)-...
                alphaD*mPA(:,1:PreP+PostP+1);
            
            %update bar
            %----------
%             percentdone=floor(p*100/(peakcount-1));
%             if floor(percentdone)>=barth
%                 if percentdone>=25 & Flag25==0
%                     fprintf('25%% ')
%                     Flag25=1;
%                 elseif percentdone>=50 & Flag50==0
%                     fprintf('50%% ')
%                     Flag50=1;
%                 elseif percentdone>=75 & Flag75==0
%                     fprintf('75%% ')
%                     Flag75=1;
%                 elseif percentdone==100
%                     fprintf('100%%\n')
%                 else
%                     fprintf('.')
%                 end
%                 
%                 while barth<=percentdone
%                     barth=barth+barth_step;
%                 end
%                 if barth>100
%                     barth=100;
%                 end
%             end
        end
        
        %---------------- Processing of last GHW peaks------------------
        sectionPoints=samples-(peakI(peakcount-GHW)-PreP)+1;
        
        if samples-peakI(peakcount) >= PostP
            
            alphaV=sum(detrend(data(:,peakI(peakcount)-...
                PreP:peakI(peakcount)+PostP)','constant')'.*...
                mPA(:,1:PreP+PostP+1),2)./...
                sum(mPA(:,1:PreP+PostP+1).*mPA(:,1:PreP+PostP+1),2);
            
            alphaD=diag(alphaV);
            
            clean(:,peakI(peakcount)-PreP:peakI(peakcount)+PostP)=...
                data(:,peakI(peakcount)-PreP:peakI(peakcount)+PostP)- ...
                alphaD*mPA(:,1:PreP+PostP+1);
            
            endpoints=samples-peakI(peakcount)-PostP;
            
            clean(:,samples-endpoints+1:samples)=...
                data(:,samples-endpoints+1:samples);
            
        else
            
            alphaV=sum(detrend(data(:,peakI(peakcount)-PreP:samples)','constant')'.*...
                mPA(:,1:samples-(peakI(peakcount)-PreP)+1),2)./...
                sum(mPA(:,1:samples-(peakI(peakcount)-PreP)+1).*...
                mPA(:,1:samples-(peakI(peakcount)-PreP)+1),2);
            
            alphaD=diag(alphaV);
            
            clean(:,peakI(peakcount)-PreP:samples)=...
                data(:,peakI(peakcount)-PreP:samples)-...
                alphaD*mPA(:,1:samples-(peakI(peakcount)-PreP)+1);
        end
        
        
        
        for p=peakcount-GHW:peakcount-1
            
            alphaV=sum(detrend(data(:,peakI(p)-PreP:peakI(p)+PostP)','constant')'.*...
                mPA(:,1:PreP+PostP+1),2)./...
                sum(mPA(:,1:PreP+PostP+1).*mPA(:,1:PreP+PostP+1),2);
            
            alphaD=diag(alphaV);
            
            clean(:,peakI(p)-PreP:peakI(p)+PostP)=...
                data(:,peakI(p)-PreP:peakI(p)+PostP)-...
                alphaD*mPA(:,1:PreP+PostP+1);
            
            %update bar
            %----------
%             percentdone=floor(p*100/(peakcount-1));
%             if floor(percentdone)>=barth
%                 if percentdone>=25 & Flag25==0
%                     fprintf('25%% ')
%                     Flag25=1;
%                 elseif percentdone>=50 & Flag50==0
%                     fprintf('50%% ')
%                     Flag50=1;
%                 elseif percentdone>=75 & Flag75==0
%                     fprintf('75%% ')
%                     Flag75=1;
%                 elseif percentdone==100
%                     fprintf('100%%\n')
%                 else
%                     fprintf('.')
%                 end
%                 
%                 while barth<=percentdone
%                     barth=barth+barth_step;
%                 end
%                 if barth>100
%                     barth=100;
%                 end
%             end
        end
        data=clean;
        
end
return;
