%plot muscle motion, unity lines, and rotational velocity

%quantitative comparisons
%1. loads in data from segmented muscles of each larva
%2. finds centroid of each ROI over time
%3. plot bend side vs. curve side activity:
    %b. categorize: muscle peak <=50% muscle present, >50% muscle present
        %i. divide by bend-roll type:
            %RCW bend, to <=50% = bend
            %RCCW bend, away <=50% = bend
            %LCW bend, away <=50% = stretch
            %LCCW bend, to <=50% = stretch %b/c all of these vals are set
            %based on y position NOT time because diff rotational
            %velocities per different muscles
    %c. take mean of ratiometric val's for bend and stretch frames
    %d. plot unity line plot of the mean ratiometric val's for each muscle:
        %i. once for each indiv larva, all muscles
        %ii. once for all combined
%4. plot rotational velocity of different muscles on average per larva
    %where x is each muscle or muscle group across segments
    %in a single video and y is the avg combined velocity of the muscle
    %across segments/muscle groups (periodicity to rolling?)
        
%Patricia Cooney, 04/2022, 5/2022, 6/2022 
%Grueber Lab
%Columbia University

%% UPDATE: 5/22
%do raw ratio vals instead of drr:
%comment out velocity/position for now
%also tried separating values on each side further by calculating
%exponential of ratio values and by taking the top 50% and top 25% of bend inds and top
%50% and top 25% of stretch inds, but did not separate much further

%% 6/22
%divide into 4 chunks of the traces (stretch, intermed 1, intermed 2, and
%bend positions) and just plot stretch vs. bend
 
%% load ROIs per larva
clear all
close all

larvanames = {'larva2_run2','larva3_run1','larva4_run2_eL2','larva10_run1_L2','larva10_run4_L2'};
color = {'r','b','m','g','k'};
dx = -0.001;
dy = -0.001;

for larva = 1:length(larvanames)
    clear cents
    clear centsy
    clear centsymois
    clear medianratiobend
    clear medianratiostretch
    clear meanratiobend
    clear meanratiostretch
    clear mid
    clear polyin
    cd 'C:\Users\coone\Desktop\Patricia\newdataset-3D';
    larvafolder = [pwd,'\',larvanames{larva}];
    larvaf = dir([larvafolder '/**/','*reapprois-0522']);
    cd(strcat(larvafolder,'\',larvaf.name));
    
    matfilenames = dir([pwd,'/*recalcreapp','*_avgdata.mat']);
    test = load(matfilenames(1).name);
    frames = length(test.avgdata.newtraces);

    dorsalmuscs = {'01','09','10'};
    lateralmuscs = {'11','19','20','04','05'};
    ventralmuscs = {'12','27','15','16','17'};
    allmuscs = [dorsalmuscs, lateralmuscs, ventralmuscs]; %indices here are in order of DV, then go by left right

    ratiomois = nan(length(allmuscs),2,2,3,length(larvanames)); %muscle,medians,lr,seg,larva
    
    a = 1;
    b = 1;
    for m = 1 : length(matfilenames)
        clear loaded
        if contains(matfilenames(m).name,'D')||contains(matfilenames(m).name,'L')||contains(matfilenames(m).name,'V')
            %do nothing
        else
            loaded = load(matfilenames(m).name,'avgdata');
            newtraces = loaded.avgdata.newtraces; 
            ratio = (newtraces(:,4));
            
            if larva == 4 || larva == 5
                bend = 1;
            elseif larva == 1 || larva == 2 || larva == 3
                bend = 2;
            end
            
            musclename = matfilenames(m).name(13:17);
            expername = loaded.avgdata.expername;
            rois = loaded.avgdata.newrois;
            
            if max(contains(allmuscs(:),musclename(4:5))>0)
                moi = 1;
                muscnum = find(contains(allmuscs(:),musclename(4:5)));
                if contains(musclename,'l')
                    lr = 1;
                else
                    lr = 2;
                end
                seg = str2double(musclename(2))-1;
            end
            
            for f = 1:frames
                if ~isempty(rois{f})
                    %find centroid positions at each frame
                    polyman = rois{f};
                    cents(f,:) = mean(polyman,1,'omitnan');%find centroids by taking mean of x and y vals
                    centsy(f) = cents(f,2);
                else
                    cents(f,:) = [nan, nan];
                    centsy(f) = nan;
                end
            end
            
            if moi == 1
                centsymois(muscnum,:,seg,lr) = centsy(:);
                b = b+1;
            end
            
%             %midway point of centroid for unity plots
%             %top of image position yvals < bottom of image position yvals
%             %organize ascending so that all oriented top to bottom
%             actcentsy = sort(centsy(~isnan(centsy)),'ascend'); %organize such that all traces are organized by position
%             if length(actcentsy) >= 4
%                 mid(a) = median(actcentsy,'omitnan'); %store midpoint y vals
%                 closestIndex = find(actcentsy == mid(a));
%                 if rem(size(actcentsy,2),2)==0
%                     midrpt = repmat(mid(a),[1 length(actcentsy)]);
%                     [minValue,closestIndex] = min(abs(midrpt-actcentsy));
%                     mid(a) = actcentsy(closestIndex);
%                 end
%                 firstq = actcentsy(1:closestIndex);
%                 quart(a) = median(firstq,'omitnan');
%                 qrpt = repmat(quart(a),[1 length(actcentsy)]);
%                 [minValq,closestIndq] = min(abs(qrpt-actcentsy));
%                 trueq(a) = actcentsy(closestIndq);
% 
%                 thirdq = actcentsy(closestIndex:length(actcentsy));
%                 tert(a) = median(thirdq,'omitnan');
%                 trpt = repmat(tert(a),[1 length(actcentsy)]);
%                 [minValt,closestIndt] = min(abs(trpt-actcentsy));
%                 truet(a) = actcentsy(closestIndt);         	
% 
%                 ratio(ratio==0) = nan;
%                 %find mean ratio for the stretch vs. bend side
%                 %larva2 = RCCW; <mid = bend; >mid = stretch
%                 %larva3 = RCCW; <mid = bend; >mid = stretch
%                 %larva4 = RCW; <mid = bend; >mid = stretch
%                 %larva10-1 = LCW; >mid = bend; <mid = stretch;
%                 %larva10-4 = LCCW; >mid = bend; <mid = stretch;
%                 if bend == 1 %L to and away
%                     indsbend = find(centsy >= truet(a));
%                     indsstretch = find(centsy <= trueq(a));
%                 elseif bend == 2 %R to and away
%                     indsbend = find(centsy <= trueq(a));
%                     indsstretch = find(centsy >= truet(a));
%                 end
%                 meanratiobend(a) = mean(ratio(indsbend),'omitnan');
%                 meanratiostretch(a) = mean(ratio(indsstretch),'omitnan');
%                 %% figures
%                 %plot unity line + ratiocomp for indiv larva
%                 figure(larva)
%                 scatter(meanratiobend(a),meanratiostretch(a),color{larva})
%                 text(meanratiobend(a)+dx, meanratiostretch(a)+dy, extractBetween(musclename,3,5))
%                 hold on
% 
%                 %plot unity line + ratiocomp for all larvae
%                 figure(8)
%                 scatter(meanratiobend(a),meanratiostretch(a),color{larva})
%                 text(meanratiobend(a)+dx, meanratiostretch(a)+dy, extractBetween(musclename,3,5))
%                 hold on
% 
%                 %store for the average to or average away plots
%                 if moi == 1
%                     ratiomeanmois(muscnum,:,lr,seg,larva) = [meanratiobend(a),meanratiostretch(a)];
%                 end
             end
%             a = a+1;
          %end
      end
%     %finish out the figure details
%     figure(larva); %all for single larva
%     %make unity line
%     xvals = [0, 3];
%     yvals = [0, 3];
%     plot(xvals,yvals,'--')
%     title(strcat('Ratio Values - ',expername))
%     xlabel('Mean Ratio for Frames in Bend')
%     ylabel('Mean Ratio for Frames in Stretch')
%     hold off
%     saveas(gcf,strcat(expername,'mean_ratio_unityplots_quarts'),'jpeg')
%     saveas(gcf,strcat(expername,'mean_ratio_unityplots_quarts','.fig'))
%     saveas(gcf,strcat(expername,'mean_ratio_unityplots_quarts'),'svg')
%end
% 
% %finish out the figure details
% figure(8) %all
% plot(xvals,yvals,'--')
% title('Mean Values for All Larvae')
% xlabel('Mean Ratio for Frames in Bend')
% ylabel('Mean Ratio for Frames in Stretch')
% saveas(gcf,'all_mean_ratio_unityplots_quarts','jpeg')
% saveas(gcf,'all_mean_ratio_unityplots_quarts.fig')
% saveas(gcf,'all_mean_ratio_unityplots_quarts','svg')
% %%
% %plot averages of all larvae
% %ratiomois(muscnum,:,lr,seg,larva) = [medianratiobend(a),medianratiostretch(a)]; %muscle,medians,lr,seg,larva
% meanratiomois = mean(ratiomeanmois,5,'omitnan');
% meancombratiomois = mean(meanratiomois,4,'omitnan'); %remove grp by segs
% meansinglevalratiomois = mean(meancombratiomois,3,'omitnan'); %remove LR
% 
% figure(9)
% colorseg = {'b','g','m'};
% side = {'l','r'};
% for se = 1:3
%     for si = 1:2
%         for mn = 1:size(meanratiomois,1)
%             scatter(meanratiomois(mn,1,si,se),meanratiomois(mn,2,si,se),colorseg{se});
%             text(meanratiomois(mn,1,si,se)+dx, meanratiomois(mn,2,si,se)+dy,strcat(num2str(se+1),side{si},allmuscs{mn}));
%             hold on
%         end
%     end
% end
% 
% figure(9) %avg all w/ seg and side
% plot(xvals,yvals,'--')
% title('Combined Mean Ratio Values for All Larvae')
% xlabel('Average of Mean Ratio for Frames in Bend')
% ylabel('Average of Mean Ratio for Frames in Stretch')
% saveas(gcf,'mean_ratio_unityplots_avgallbyseg_redo_quarters','jpeg')
% saveas(gcf,'mean_ratio_unityplots_avgallbyseg_redo_quarters.fig')
% saveas(gcf,'mean_ratio_unityplots_avgallbyseg_redo_quarters','svg')
% 
% %average all (combined seg and side)
% figure(10)
% for mn = 1:size(meansinglevalratiomois,1)
%     scatter(meansinglevalratiomois(mn,1),meansinglevalratiomois(mn,2));
%     text(meansinglevalratiomois(mn,1)+dx, meansinglevalratiomois(mn,2)+dy,allmuscs{mn});
%     hold on
% end
% 
% figure(10) %avg all
% plot(xvals,yvals,'--')
% title('Combined Mean Ratio Values for All Larvae')
% xlabel('Average of Mean Ratio for Frames in Bend')
% ylabel('Average of Mean Ratio for Frames in Stretch')
% saveas(gcf,'mean_ratio_unityplots_avgallcombined_redo_quarters','jpeg')
% saveas(gcf,'mean_ratio_unityplots_avgallcombined_redo_quarters.fig')
% saveas(gcf,'mean_ratio_unityplots_avgallcombined_redo_quarters','svg')
% 
% %% post hoc regression to show diff b/t unity and actual activity
% mdl = fitlm(meansinglevalratiomois(:,1),meansinglevalratiomois(:,2));
% figure(9)
% plot(mdl)
% 
% figure(10)
% plot(mdl)
% save('forregression_unityinfo_quarters_allcomb.mat','mdl');

%% taking multiple muscles across the same larva, plot ex positions - seg sync
%plot position raw, w/ time trace color for each mois muscle on sep fig
xin = [1, 3, 5];
xp = [0, 1];

for posplo = 1:size(centsymois,1)
    for se = 1:3
        for si = 1:2
            figure(posplo)
            %L2,R2,L3,R3,L4,R4 plot order, plus color for timing
            y = centsymois(posplo,:,se,si);
            x = repmat(xin(se)+xp(si),size(y));
            z = zeros(size(y)); %check what this does..
            col = 1:length(y); %change this to be the index of y (aka the frame num)
            surface([x;x],[y;y],[z;z],[col;col],...
                'facecol','no',...
                'edgecol','interp',...
                'linew',2);
            hold on
        end
    end
    %find index in allmuscs to title
    tmusc = allmuscs(posplo);
    title(strcat('Position across segments: ',tmusc))
    hold on
    xlim = [0,xin(end)+xp(end)+1];
    xlabel('Segments');
    ylim = max(max(max(centsymois)));
    ylabel('Position through Time')
    colorbar
    
    saveas(gcf,strcat('musclepositons_',tmusc{1}),'jpeg')
    saveas(gcf,strcat('musclepositions_',tmusc{1},'.fig'))
    close all
    pause(1)
end

%% with same ex muscles, make beeplots of the velocities all on one
%calculate average velocity for the muscles across all 3 segs:
velmusctime = nan(size(centsymois,1),3,2);
velmusc = nan(size(centsymois,1),2);

xml = 1:2:26;
xmr = 2:2:26;
ybee = nan(3,length(allmuscs)*2);
xbee = nan(size(ybee));

figure
for mus = 1:size(centsymois,1)
    for side = 1:2
        for se = 1:3
            velmusc = abs(diff(centsymois(mus,:,se,side)));%find velocity from subtracting y position vals
            velmusctime(mus,se,side) = nanmean(velmusc,2);%find avg vel trace for muscle in each seg
        end
        if side == 1
            ybee(:,mus) = velmusctime(mus,:,side)';
            xbee(:,mus) = repmat(xml(mus),[3 1])';
        else
            ybee(:,mus+13) = velmusctime(mus,:,side)';
            xbee(:,mus+13) = repmat(xmr(mus),[3 1])';
        end
    end
end
% xbee(isnan(ybee))=0;
% ybee(isnan(ybee))=0;
xbeer = reshape(xbee,(size(xbee,1)*size(xbee,2)),1);
ybeer = reshape(ybee,(size(ybee,1)*size(ybee,2)),1);
% beeswarm(xbeer,ybeer);
% xlabel('Muscle')
% ylabel('Average Velocity during Imaging (px/sec)')
% title('Rotational Velocity of Different Spatial Muscle Groups')

% try boxplot instead %% this is for lr together but the matrix ybee is
% 1-13 left and 14 - 26 right

for add = 1:length(allmuscs)*2
    if add <= 13
        newx = xml(add);
        xgrps(newx) = strcat('l',allmuscs(add));
    elseif add > 13
        newx = xmr(add - 13);
        xgrps(newx) = strcat('r',allmuscs(add-13));
    end
end

xgrpsredo = [xgrps(xml),xgrps(xmr)];
%% boxplot orig order LR separate
figure
boxplot(ybee,xgrpsredo)
xlabel('Muscle')
ylabel('Average Velocity during Imaging (px/sec)')
title(strcat('Rotational Velocity of Different Spatial Muscle Groups: ',expername))
saveas(gcf,strcat(expername,'_rotvel_spatial_lrsep'),'jpeg')
saveas(gcf,strcat(expername,'_rotvel_spatial_lrsep.fig'))
%% boxplot orig order LR together
% for add = 1:length(allmuscs)*2
%     if add <= 13
%         newx = xml(add);
%         xgrps(newx) = strcat('l',allmuscs(add));
%     elseif add > 13
%         newx = xmr(add - 13);
%         xgrps(newx) = strcat('r',allmuscs(add-13));
%     end
% end

for add = 1:length(allmuscs)*2
    if add <= 13
        newy = xml(add);
        ybeelrnext(:,newy) = ybee(:,add);
    elseif add > 13
        newy = xmr(add - 13);
        ybeelrnext(:,newy) = ybee(:,add);
    end
end

figure
boxplot(ybeelrnext,xgrps)
xlabel('Muscle')
ylabel('Average Velocity during Imaging (px/sec)')
title(strcat('Rotational Velocity of Different Spatial Muscle Groups: ',expername))
saveas(gcf,strcat(expername,'_rotvel_spatial_lrtog'),'jpeg')
saveas(gcf,strcat(expername,'_rotvel_spatial_lrtog.fig'))

% % circum order boxplot
% %old order: dorsalmuscs = {'01', '09', '10'}; lateralmuscs = {'04', '05', '11', '19', '20'};
% ventralmuscs = {'12', '15', '16', '17', '27'};
% dorsalmuscs = {'01', '09', '10'};
% lateralmuscs = {'11', '19', '20', '04', '05'};
% ventralmuscs = {'12', '27', '17', '16', '15'};
% allmuscsleft = [dorsalmuscs, lateralmuscs, ventralmuscs];
% allmuscsright = flip(allmuscsleft);
% allmuscscirc = [allmuscsleft, allmuscsright];
% ltransmat = [1:3, [6, 7, 8, 4, 5], [9, 13, 12, 11, 10]];
% rtransmat = flip(ltransmat);
% fulltransmat = [ltransmat, rtransmat];
% 
% for yu = 1:size(ybee,2)
%     if yu <= 13
%         newbox(:,yu) = ybee(:,fulltransmat(yu));
%     elseif yu > 13 
%         newbox(:,yu) = ybee(:,fulltransmat(yu)+13);
%     end
% end
% 
% for add = 1:length(allmuscs)*2
%     if add <= 13
%         newx = xml(add);
%         xgrps(newx) = strcat('l',allmuscs(add));
%     elseif add > 13
%         newx = xmr(add - 13);
%         xgrps(newx) = strcat('r',allmuscs(add-13));
%     end
% end
% figure
% boxplot(ybee,xgrps)
% xlim([0 27])
% xlabel('Muscle')
% ylabel('Average Velocity during Imaging (px/sec)')
% title(strcat('Circumferential Order Rotational Velocity of Different Spatial Muscle Groups: ',expername))
end