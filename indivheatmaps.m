%% Individual Heatmaps for each larva

%heatmaps for all ratio traces in most consistently measured muscles

%Patricia Cooney, 6/2022
%Grueber Lab
%Columbia University

%% load ROIs per larva
clear all
close all

larvanames = {'larva2_run2','larva3_run1','larva4_run2_eL2','larva10_run1_L2','larva10_run4_L2'};
larvae = struct;

%pull out the muscles of interest (most consistently measured)
dorsalmuscs = {'01','09','10'};
lateralmuscs = {'11','19','20','04','05'};
ventralmuscs = {'12','27','15','16','17'};
allmuscs = [dorsalmuscs, lateralmuscs, ventralmuscs];

for larva = 1:length(larvanames)    
    clear frames
    cd 'C:\Users\coone\Desktop\Patricia\newdataset-3D';
    larvafolder = [pwd,'\',larvanames{larva}];
    larvaf = dir([larvafolder '/**/','*reapprois-0522']);
    cd(strcat(larvafolder,'\',larvaf.name));

    matfilenames = dir([pwd,'/*recalcreapp','*_avgdata.mat']);
    test = load(matfilenames(20).name);
    frames = length(test.avgdata.newtraces);
    
    allratios = nan(13,frames,2,3); %musclenum,maxtracelength,lr,seg

    for m = 1:length(matfilenames)
        clear ratio
        clear musclename
        musclename = matfilenames(m).name(13:17);
        if max(contains(allmuscs(:),musclename(4:5))>0)
            loaded = load(matfilenames(m).name,'avgdata');
            newtraces = loaded.avgdata.newtraces; 
            ratio = (newtraces(:,4)'); %make it a row

            %sort out details for indexing
            seg = str2double(musclename(2))-1;
            muscnum = find(contains(allmuscs(:),musclename(4:5)));
            if contains(musclename,'l')
                lr = 1;
            else
                lr = 2;
            end
            allratios(muscnum,:,lr,seg) = ratio;
        end
    end
    
    %heatmap for indiv larvae: all muscs, sep segs, sep sides
    figure
    
    pp = parula;
    pp(1,:) = 0; %set frames where muscle is out of field to black
    
    for lr = 1:2
        for seg = 1:3
            if lr == 2
                subplot(3,2,lr*seg)
            else
                subplot(3,2,(lr+1)*seg-1)
            end
            imagesc(allratios(:,:,lr,seg))
            colormap(pp)
            yticks(1:13)
            yticklabels(allmuscs);
            if lr == 1 && seg == 1
                title('Left Side')
            elseif lr == 2 && seg == 1
                title('Right Side')
            elseif seg == 2
                ylabel('Muscle Ratio Traces')
            elseif seg == 3
                xlabel('Frames')
            end
            hold on
        end
        hold on 
    end
    sgtitle(strcat('Heat Map of Muscles across Segments - ',larvanames{larva}));
    
    %SAVE
    saveas(gcf,strcat(larvanames{larva},'_indivheatmaps_black'),'jpeg')
    saveas(gcf,strcat(larvanames{larva},'_indivheatmaps_black'),'svg')
    saveas(gcf,strcat(larvanames{larva},'_indivheatmaps_black','.fig'))
end
save('alllarvae_allratios.mat','larvae');