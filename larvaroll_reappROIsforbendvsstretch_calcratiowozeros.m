%% reapply rois for side calculations comparing behaviors

%Patricia Cooney, 6/2022
%Grueber Lab
%Columbia University

%% Locate the tiffs of interest
clear all
close all
    
larvanames = {'larva2_run2','larva3_run1','larva10_run1_L2'};

larva2frames = [nan, nan, nan; nan, 283, 325];
larva3frames = [1, 4, 42; 278, 325, 512];
larva10frames = [354, 397, 419; nan, 16, 36];

allframes = cat(3,larva2frames,larva3frames,larva10frames);

behavs = {'bend','init','roll'};
sides = {'DLs','VOs'};
bs = {'bs','ss'};

for larva = 1:length(larvanames)
    %clear previous variables
    clear matfilenames
    clear filenames
    clear i_tempgreen
    clear i_tempred
    
    %load in single larvae and their muscle mats
    cd 'C:\Users\coone\Desktop\Patricia\newdataset-3D';
    larvafolder = [pwd,'\',larvanames{larva}];
    larvaf = dir([larvafolder '/**/','*reapprois-0522']);
    cd(strcat(larvafolder,'\',larvaf.name));
    
    filenames = dir('*.tif');
    file = filenames(1);
    expername = extractBetween(convertCharsToStrings(file.name),"gcamp_",".tif");

    %% Load in the separate channels
    %load
    for fi = 1:length(filenames)
        if contains(filenames(fi).name,'fullrun') && contains(filenames(fi).name,'gcamp')
            tiffnamegreen = filenames(fi).name;
            newinfo = imfinfo(tiffnamegreen);
            for a = 1 : 600
                if larva == 4
                    i_tempgreen(:,:,a) = imread(tiffnamegreen, a)./40; %adjustment line for when gcamp is far brighter than mcherry (larva10_run1 = ./40; larva10_run4 = ./40)
                else
                    i_tempgreen(:,:,a) = imread(tiffnamegreen, a);
                end
            end
        elseif contains(filenames(fi).name,'fullrun') && contains(filenames(fi).name,'mcherry')
            tiffnamered = filenames(fi).name;
            newinfo = imfinfo(tiffnamered);
            for a = 1:600
                if larva == 3 
                    i_tempred(:,:,a) = imread(tiffnamered, a)./50; %adjust divide by 50 for larva4 expression
                else
                    i_tempred(:,:,a) = imread(tiffnamered, a);
                end
            end
        else
        end
    end
    
    %load those mat files and check if already done and in spreadsheet 
    matfilenames = dir([pwd,'/**/','*recalcreapp_','*_avgdata.mat']);

    for m = 1:length(matfilenames)
        musclename = matfilenames(m).name(13:24);

        if max(contains(behavs,matfilenames(m).name(21:24)))==1
            clear loaded
            loaded = load(matfilenames(m).name,'avgdata');
            newrois = loaded.avgdata.newrois;
            
            km = find(contains(behavs,matfilenames(m).name(21:24))==1);
            si = find(contains(sides,matfilenames(m).name(14:16))==1);
            
            frames = allframes(si,km,larva):allframes(si,km,larva)+2;

            %just reapply all the old ROIs and check trace
            [newtraces] = runseg(i_tempred,i_tempgreen,newrois,frames);
            %plot those vals
            plottraces(newtraces,musclename);
            
            %save your values
            close all
            disp('Ok, saving all')
            avgdata.expername = expername;
            avgdata.musclename = musclename;
            avgdata.newtraces = newtraces;
            avgdata.newrois = newrois;

            save(strcat('recalcreapp_',musclename,'_avgdata_wozeros.mat'),'avgdata','-v7.3')
            disp('Done')
        end   
    end
end
%%
%%
%%
%% local functions:
function [newtraces] = runseg(i_tempred,i_tempgreen,newrois,frames)
    
    for foi = 1:size(newrois,2)
        fr = frames(foi);
        i_segred = i_tempred(:,:,fr);
        i_seggreen = i_tempgreen(:,:,fr);

        vertsarray = cell2mat(newrois(foi));          
        newmred = poly2mask(vertsarray(:,1),vertsarray(:,2),size(i_segred,1),size(i_segred,2));
        newmgreen = poly2mask(vertsarray(:,1),vertsarray(:,2),size(i_seggreen,1),size(i_seggreen,2));

        %remove mcherry puncta and calculate mean values
        %all channels
        withpunctared_90 = double(i_segred(newmred));
        withpunctagreen_90 = double(i_seggreen(newmgreen));
        [~,reminds] = rmoutliers(withpunctared_90,'percentiles',[0,90]); %remove any fluor val's that are stat outlier for ROI in each frame
        wopunctared_90 = withpunctared_90(~reminds);
        wopunctagreen_90 = withpunctagreen_90(~reminds);

        %UPDATE: set zeros to nans 
        wopunctared_90(wopunctared_90==0) = nan;
        wopunctagreen_90(wopunctagreen_90==0) = nan;

        newtraces(foi,1) = nanmean(wopunctagreen_90);
        newtraces(foi,2) = nanmean(i_segred(newmred));
        newtraces(foi,3) = nanmean(wopunctared_90); 

        %calculate the ratio of the ROI, then take the mean
        roiratio = (wopunctagreen_90./wopunctared_90);
        newtraces(foi,4) = nanmean(roiratio);
    end
end
   
%%
function plottraces(newtraces,musclename) 
%isolate just the data from when muscles in the FOV
if find(newtraces(:,1)<=0,1,'first') > find(newtraces(:,1)>0,1,'first')
    spaceind = find(newtraces(:,1)<=0,1,'first');
    nonzerotraces = reshape(newtraces(newtraces>0),[],4);
    nonzerotraces = [nonzerotraces(1:spaceind-1,1:4); zeros(5,4); nonzerotraces(spaceind:end,1:4)];
else
    nonzerotraces = reshape(newtraces(newtraces>0),[],4);
end

    figure
    colororder({'k','b'})

    yyaxis left
    plot(newtraces(:,1),'g');
    hold on
    plot(newtraces(:,3),'r');
    hold on
    ylabel('Raw Traces')

    yyaxis right
    plot(newtraces(:,4),'b');
    hold on
    ylabel('Ratio Raw')

    xlabel('Frames')
    title(strcat(musclename,' Roll Activity'))

    saveas(gcf,strcat(musclename,'_behavcomp_wozeros'),'jpeg')
    saveas(gcf,strcat(musclename,'_behavcomp_wozeros.fig'))
end