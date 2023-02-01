% Boxplots comparing sides for bend, initiate, roll

%Patricia Cooney, 6/2022
%Grueber Lab
%Columbia University

%%
clear all
close all

%load mat files for each larva
larvanames = {'larva2_run2','larva3_run1','larva10_run1_L2'};

behavs = {'bend','init','roll'};
sides = {'DLs','VOs'};
bs = {'bs','ss'};

%array for behavs,bendvstrethc,D or V muscs,larva
ratiovals = nan(1,length(behavs)*length(bs),length(sides),length(larvanames));
allratiovals = nan(length(larvanames),length(behavs)*length(bs),length(sides));

for larva = 1:length(larvanames)
    %load in single larvae and their muscle mats
    cd 'C:\Users\coone\Desktop\Patricia\newdataset-3D';
    larvafolder = [pwd,'\',larvanames{larva}];
    larvaf = dir([larvafolder '/**/','*reapprois-0522']);
    cd(strcat(larvafolder,'\',larvaf.name));
    
    matfilenames = dir([pwd,'/*recalcreapp','*_avgdata_wozeros.mat']);
    test = load(matfilenames(1).name);

    for m = 1:length(matfilenames)
        %check if init, roll , or bend and sort
        if max(contains(behavs,matfilenames(m).name(21:24)))==1
            km = find(contains(behavs,matfilenames(m).name(21:24))==1);
            if km == 1
                km = 0;
            elseif km == 3
                km = 4;
            end  
            bt = find(contains(bs,matfilenames(m).name(18:19))==1);
            si = find(contains(sides,matfilenames(m).name(14:16))==1);
            %load in the ratio trace and store in appropriate index
            loaded = load(matfilenames(m).name,'avgdata');
            newtraces = loaded.avgdata.newtraces; 
            ratiovals(1:3,km+bt,si,larva) = newtraces(:,4);
        end
    end 
    if larva == 1
        startlarva = 1;
    elseif larva == 2
        startlarva = 4;
    elseif larva == 3
        startlarva = 7;
    end
    allratiovals(startlarva:startlarva+2,:,:) = (ratiovals(:,:,:,larva)); %row = larva obs, cols = behav b s, behav bs, ...; 3rd dim = dv
end

allratiovals(allratiovals == 0) = nan;
%% for all larvae, dorsal and ventral sep then dv combo
%run stats dorsal and ventral sep
for dv = 1:2
    %make boxplots
    figure(dv)
    for k = 1:size(allratiovals,2) 
        j = k;
        while j <= 6
            [~,pdv(k,j,dv)] = ttest2(allratiovals(:,k,dv),allratiovals(:,j,dv));
            j = j + 1;
        end
    end
    boxplot(allratiovals(:,:,dv));
    title(strcat('Bend vs. Stretch Side',sides(dv)))
    saveas(gcf,strcat('compbehavmoments_',sides{dv}),'jpeg')
    saveas(gcf,strcat('compbehavmoments_',sides{dv}),'svg')
    saveas(gcf,strcat('compbehavmoments_',sides{dv},'.fig'))
end

%combined ranksum and boxplots
combdvratiovals = mean(allratiovals,3,'omitnan');
figure(dv + 1)
for k = 1:size(combdvratiovals,2)
    j = k;
    while j <= 6
        [~,pcomb(k,j)] = ttest2(combdvratiovals(:,k),combdvratiovals(:,j));
        j = j + 1;
    end
end
boxplot(combdvratiovals);
title(strcat('Bend vs. Stretch Side - DV combined'))
saveas(gcf,strcat('compbehavmoments_combined'),'jpeg')
saveas(gcf,strcat('compbehavmoments_combined'),'svg')
saveas(gcf,strcat('compbehavmoments_combined.fig'))

%save organized data and pvals
save('compbehavmoments.mat','ratiovals','allratiovals','combdvratiovals','pcomb','pdv');