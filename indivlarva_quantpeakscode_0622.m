%% Quantitative Peaks Code - linear time

%Find and plot peak ratiometric values for muscles of interest in each segment on
%3D averaged, ratiometric data

%then compare to activity from Aref's crawling functional muscle data

%Update to orig: peakscodedistlineslr.m & peakscoderationov.m
%PC 8/2021; 11/2021; 3/2022; 5/2022; finalized 6/2022

%% load background data
clear all
close all

larvanames = {'larva2_run2','larva3_run1','larva4_run2_eL2','larva10_run1_L2','larva10_run4_L2'};
larvae = struct;

for larva = 1:length(larvanames)
    %roll line structures
    lines = nan(3,2,60,2); %seg,L & R muscles
    
    %load in single larvae and their muscle mats
    clear peak_time
    clear ratio
    clear musclename
    clear frames
    
    cd 'C:\Users\coone\Desktop\Patricia\newdataset-3D';
    larvafolder = [pwd,'\',larvanames{larva}];
    larvaf = dir([larvafolder '/**/','*reapprois-0522']);
    cd(strcat(larvafolder,'\',larvaf.name));
    
    matfilenames = dir([pwd,'/*recalcreapp','*_avgdata.mat']);
    test = load(matfilenames(1).name);
    frames = length(test.avgdata.newtraces);
    
    q = 1;
    for m = 1:length(matfilenames)
        loaded = load(matfilenames(m).name,'avgdata');
        newtraces = loaded.avgdata.newtraces; 
        ratio = (newtraces(:,4));
        musclename = matfilenames(m).name(13:17);
        expername = loaded.avgdata.expername;

        %make the arrays for quick plotting
        full_names{q} = musclename;
        [~, peak_time(q)] = max(ratio);
        q = q+1;
    end
end
    %% get ready to plot
    %offset for text
    dx = 0.1;
    dy = -0.025;
    
    %%
    %add offset values so we can actually read the data
    %data_count_indiv = zeros(4,3); %for dlv lines
    data_count_indiv = zeros(3,1,2);

    %calculate offsets for the lines - preallocate structs
    segy = nan(1,size(peak_time,2));
    c = nan(size(peak_time,2),3);
    cmap_indiv = colormap(parula(60));
    
    %lists for the indexing below
    %sidelist = {'l','r'};

    %calc xy coords for data points, single larva
    for i = 1:size(peak_time,2)
        %set the y-axis based on segment, side, and muscgrp
        %seg
        segnum = str2double(extractBetween(full_names(i),2,2));
        %side
        if contains(full_names(i),'l')
            lr = 1;
        elseif contains(full_names(i),'r')
            lr = 2;
        end
        mnum = str2double(extractBetween(full_names(i),4,5));
            
            if ~isempty(mnum) && mnum~=0
                segy(i) = segnum + data_count_indiv(segnum-1,1,lr);
                lines(segnum-1,mnum,lr,:) = [peak_time(i), segy(i)]; %store as line

                data_count_indiv(segnum-1,1,lr) = data_count_indiv(segnum-1,1,lr) + 0.1; %update for jitter

                %plot the peaks for each muscle across segments
                c(i,:) = cmap_indiv(mnum,:);
                figure(1)
                scatter(peak_time(i),segy(i),[],c(i,1:3),'filled')
                text(peak_time(i)+dx, segy(i)+dy, full_names(i))
                ylim([1.5 5.5])
                hold on 
            else
            end
    end

    lines(lines==0)=nan; %3, 30, 2, 2 %seg, muscle, side, tr
    
    %finish lines for indiv muscles
    for l = 1:size(lines,2)
        for s = 1:size(lines,3)
            for se = 1:size(lines,1)-1
                figure(1)
                plot(lines(se:se+1,l,s,1),lines(se:se+1,l,s,2),'k');
                hold on
            end
        end
    end

    figure(1)
    title(strcat('Highest Captured Ratio Per Muscle - ',expername));
    xlabel('Frames');
	ylabel('Segments (a2 - a4)'); 
    
    saveas(gcf,strcat(expername,'_','indiv_peaks'),'jpeg')
    saveas(gcf,strcat(expername,'_','indiv_peaks'),'svg')
    saveas(gcf,strcat(expername,'_','indiv_peaks.fig'))

    larvae(larva).peak_time = peak_time;
    larvae(larva).lines = lines;
    larvae(larva).frames = frames;
    close all
%save the nice struct
save('alllarvae_lines.mat','larvae');

%%
%keep going to dists comparisons?
clear all
kg = input('Compare distances now? 1 = yes');

%load back in roll
load('alllarvae_lines.mat')

%load in the crawl data
if kg == 1
    crawldata = load('redoratio_crawl-887_allrois_wlines.mat');
    clines = crawldata.allrois.lines;
    for m = 1:30
        for n = 1:2
            crawldists(m,n) = clines(n+1,m,1) - clines(n,m,1); %calculate diff between same musc in a2-a3 and a3-a4
        end
    end
    avgcrawldists = nanmean(crawldists,2);
    avgcrawldistsperc = avgcrawldists./crawldata.allrois.frames; %normtime
    avgcrawldiststimes = avgcrawldists.*0.2; %convert to seconds, 5 fps
end 

alltimes = nan(30,5*2);
allpercs = nan(30,5*2);
d = 1;
f = d+1;
        
for lar = 1:size(larvae,2)
    dists = nan(60,2); %1-30 = left, 31-60 = right, 1 = a2-a3, 2 = a3-a4

    %reshape to combine l & r into single dim
    relines = reshape(larvae(lar).lines,[3,60,2]);
    frames = larvae(lar).frames;

    %Find the distances between the peak times of each muscle between segs
    for m = 1:60
        for n = 1:2
            dists(m,n) = relines(n+1,m,1) - relines(n,m,1);
        end
    end

    %convert to % roll cycle & to time, plot both because need to compare to crawl
    %roll%
    distsperc = dists./frames;
    %time
    diststimes = dists.*0.1; % 10 fps = 100ms per frame --> seconds

    %find avg dists between muscle peak times, all segments
    avgdistsperc = nanmean(distsperc,2);
    avgdiststimes = nanmean(diststimes,2);

    %% stats and meta analysis bar plot on avgdists for crawl vs roll
    pdists_perc = nan(1,2);
    pdists_times = nan(1,2);

    for s = 1:2
        if s == 1
            pdists_perc(lar,s) = ranksum(avgcrawldistsperc,avgdistsperc(1:30));  
            pdists_times(lar,s) = ranksum(avgcrawldiststimes,avgdiststimes(1:30));
        elseif s == 2
            pdists_perc(lar,s) = ranksum(avgcrawldistsperc,avgdistsperc(31:60));
            pdists_times(lar,s) = ranksum(avgcrawldiststimes,avgdiststimes(31:60));
        end
    end

    %compare to normal distribution with mu = 0 sigma = 1
    norm = normrnd(0,1,[1,30]);
    pnorm = nan(3,2);
    for r = 1:3
        if r == 1
            pnorm(r,1) = ranksum(norm,avgcrawldistsperc);
            pnorm(r,2) = ranksum(norm,avgcrawldiststimes);
        elseif r == 2
            pnorm(r,1) = ranksum(norm,avgdistsperc(1:30));
            pnorm(r,2) = ranksum(norm,avgdiststimes(1:30));
        elseif r == 3
            pnorm(r,1) = ranksum(norm,avgdistsperc(31:60));
            pnorm(r,2) = ranksum(norm,avgdiststimes(31:60));
        end
    end
    larvaroll(lar).times = [avgdiststimes(1:30)',  avgdiststimes(31:60)'];
    larvaroll(lar).perc = [avgdistsperc(1:30)',  avgdistsperc(31:60)'];
    larvaroll(lar).pvalsnorm = pnorm;
    larvaroll(lar).pvalscomp = [pdists_perc; pdists_times];

    %make summ arrays
    alltimes(1:30,d) = larvaroll(lar).times(1:30)';
    alltimes(1:30,f) = larvaroll(lar).times(31:60)';
    allpercs(1:30,d) = larvaroll(lar).perc(1:30)';
    allpercs(1:30,f) = larvaroll(lar).perc(31:60)';
    d = d + 2;
    f = d + 1;
end

%summary stats
alltimesmean = nanmean(alltimes,2);
allpercsmean = nanmean(allpercs,2);
allpnormpercs = ranksum(norm,alltimesmean);
allpnormtimes = ranksum(norm,allpercsmean);
allpcomppercs = ranksum(avgcrawldistsperc,allpercsmean);
allpcomptimes = ranksum(avgcrawldiststimes,alltimesmean);
    
figure(1)
boxplot([avgcrawldiststimes,larvaroll(1).times(1:30)',larvaroll(1).times(31:60)',larvaroll(2).times(1:30)',larvaroll(2).times(31:60)',larvaroll(3).times(1:30)',larvaroll(3).times(31:60)',larvaroll(4).times(1:30)',larvaroll(4).times(31:60)',larvaroll(5).times(1:30)',larvaroll(5).times(31:60)'],{'Crawl','Roll L - To 1','Roll R - To 1','Roll L - To 2','Roll R - To 2','Roll L - Away 1','Roll R - Away 1','Roll L - Away 2','Roll R - Away 2','Roll L - Away 3','Roll R - Away 3'}); 
ylabel('Average Time Diff (s)');
title('Crawl & Roll Intersegmental Raw Time Differences');
saveas(gcf,'comp_crawl-roll_peaktimes-sep_redocrawlratio','svg')
saveas(gcf,'comp_crawl-roll_peaktimes-sep_redocrawlratio','jpeg')
saveas(gcf,'comp_crawl-roll_peaktimes-sep_redocrawlratio.fig')

figure(2)
boxplot([avgcrawldiststimes,alltimesmean],{'Crawl','Roll'}); 
ylabel('Average Time Diff (s)');
title('Crawl & Roll Intersegmental Raw Time Differences');
saveas(gcf,'comp_crawl-roll_peaktimes-altog_redocrawlratio','svg')
saveas(gcf,'comp_crawl-roll_peaktimes-altog_redocrawlratio','jpeg')
saveas(gcf,'comp_crawl-roll_peaktimes-altog_redocrawlratio.fig')

figure(3)
boxplot([avgcrawldistsperc,larvaroll(1).perc(1:30)',larvaroll(1).perc(31:60)',larvaroll(2).perc(1:30)',larvaroll(2).perc(31:60)',larvaroll(3).perc(1:30)',larvaroll(3).perc(31:60)',larvaroll(4).perc(1:30)',larvaroll(4).perc(31:60)',larvaroll(5).perc(1:30)',larvaroll(5).perc(31:60)'],{'Crawl','Roll L - To 1','Roll R - To 1','Roll L - To 2','Roll R - To 2','Roll L - Away 1','Roll R - Away 1','Roll L - Away 2','Roll R - Away 2','Roll L - Away 3','Roll R - Away 3'});
title('Crawl & Roll Normalized Time Differences');
ylabel('Average Normalized Time Diff (%)');
saveas(gcf,'overall_crawl-roll_normpeakperc-sep_redocrawlratio','svg')
saveas(gcf,'overall_crawl-roll_normpeakperc-sep_redocrawlratio','jpeg')
saveas(gcf,'overall_crawl-roll_normpeakperc-sep_redocrawlratio.fig')

figure(4)
boxplot([avgcrawldistsperc,allpercsmean],{'Crawl','Roll'});
title('Crawl & Roll Normalized Time Differences');
ylabel('Average Normalized Time Diff (%)');
saveas(gcf,'overall_crawl-roll_rawpeaktimes-altog_redocrawlratio','svg')
saveas(gcf,'overall_crawl-roll_rawpeaktimes-altog_redocrawlratio','jpeg')
saveas(gcf,'overall_crawl-roll_rawpeaktimes-altog_redocrawlratio.fig')

save('allcomp_peakdists_pvals_redocrawlratio.mat','larvaroll');
save('comb_pvals_redocrawlratio.mat','alltimes','allpercs','allpercsmean','alltimesmean','allpnormpercs','allpnormtimes','allpcomppercs','allpcomptimes');