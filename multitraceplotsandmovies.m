%% final multitrace plots -- movies and figures

%Patricia Cooney - 6/2022
%Grueber Lab
%Columbia University
%% load key ROIs for the multitrace movies and figures
clear all
close all

%specify
larva = input('which larva? (1 - 5)');
muscs = input('which muscle? (e.g. 05, 09, 11)');

for mu = 1:length(muscs)
    moi = muscs{mu};
    vidt = strcat(moi,'LR');

    larvanames = {'larva2_run2','larva3_run1','larva4_run2_eL2','larva10_run1_L2','larva10_run4_L2'};

%% load the ROI + trace data
cd 'C:\Users\coone\Desktop\Patricia\newdataset-3D';
larvafolder = [pwd,'\',larvanames{larva}];
larvaf = dir([larvafolder '/**/','*reapprois-0522']);
cd(strcat(larvafolder,'\',larvaf.name));
matfilenames = dir([pwd,'/*recalcreapp','*_avgdata.mat']);
test = load(matfilenames(20).name,'avgdata');
frames = size(test.avgdata.newtraces,1);
muscnum = 1;

for m = 1 : length(matfilenames)
    clear loaded
    musclename = matfilenames(m).name(13:17);
    if max(contains(moi,musclename(4:5))>0)
        loaded = load(matfilenames(m).name,'avgdata');
        newtracescomb(:,muscnum) = loaded.avgdata.newtraces(:,4);
        roiscomb{muscnum,:} = loaded.avgdata.newrois;
        mois{muscnum} = musclename;
        muscnum = muscnum + 1;
    end
end

%set 0s to nans so not plotted
newtracescomb(newtracescomb == 0) = nan;

%% make the composite traces plot
plottraces(newtracescomb,mois,vidt,frames);

%% let's make a movie
%load ratio image mat data
close all
clear ratio
clear M
cd 'C:\Users\coone\Desktop\Patricia\newdataset-3D';
larvafolder = [pwd,'\',larvanames{larva},'\final3davg-reapprois-0522'];
cd(larvafolder);
ratiof = load(strcat('fullrun_',larvanames{larva},'_ratio.mat'));
ratio = ratiof.ratio;

rollframes = [[286,330];[6,44];[28,62];[400,420];[79,123]];

expername = larvanames{larva};
roll = rollframes(larva,1):rollframes(larva,2);

%image contrast params, customized across frames per larva
clow = [80,80,100,80,100]; %1 = 80; 2 = 80; 3 = 100; 4 = 80; 5 = 100
chigh = [4400,4200,600,380,880];%1 = 4400; 2 = 4200; 3 = 600; 4 = 380; 5 = 880
cax = [clow(larva),chigh(larva)];

%scale bar params for ums
pxumconv = 1.056;
ums = 100;
xsb = 1124;
ysb = 600;
col1 = round(xsb-(ums*pxumconv));
col2 = round(xsb+(ums*pxumconv));

%run it
makemovie(ratio,roll,cax,xsb,ysb,col1,col2,newtracescomb,mois,roiscomb,expername,frames,vidt);

end

%%%
%%%
%%%

%%% local functions
function plottraces(newtracescomb,mois,vidt,frames) 
    %plot all the traces - segments on different plots
    figure
    colororder({'k','b'})
    cols = {'green','magenta'};
    
    for m = 1:length(mois)
        if contains(mois{m},'a2')
            mi = 1;
        elseif contains(mois{m},'a3')
            mi = 2;
        elseif contains(mois{m},'a4')
            mi = 3;
        end
        if rem(m,2) == 0 
            colp = 2;
        else
            colp = 1;
        end
        
        subplot(3,1,mi);
        plot(newtracescomb(:,m),cols{colp});
        hold on
        ylabel('Ratio Traces')
        ylim([min(newtracescomb(:,m)), max(newtracescomb(:,m))*1.1])
        xlabel('Frames')
        xlim([1,frames])
        title(mois{m}(1:2))
    end
    %overall title
    sgtitle(vidt)
    %SAVE
    saveas(gcf,strcat(vidt,'_combtraces'),'jpeg')
    saveas(gcf,strcat(vidt,'_combtraces','.fig'))
end

%%
function makemovie(ratio,roll,cax,xsb,ysb,col1,col2,newtracescomb,mois,roiscomb,expername,frames,vidt)
    ratio(ratio < cax(1)) = cax(1);
    ratio(ratio > cax(2)) = cax(2);    

    ind = 1;
    for fr = roll(1):roll(end)
        fig = figure;
        moicols = {'green','cyan','magenta'};
        %make subplot for SCAPE data
        subplot(7,2,[1:8]);
        ratioframe = ratio(:,:,fr);
        imagesc(imrotate(ratioframe,180));
        axis off
        hold on
        for m = 1:length(mois)
            if contains(mois{m},'a2')
                mo = 1;
            elseif contains(mois{m},'a3')
                mo = 2;
            elseif contains(mois{m},'a4')
                mo = 3;
            end
            rois = roiscomb{m};
            if isempty(rois{ind}) == 0
                color = moicols{mo};
                %colored ROI based on saved verts
                vertsarray = cell2mat(rois(ind));
                %drawpolygon('Position',vertsarray);
                do = polyshape(vertsarray); 
                h = plot(do);
                h.EdgeColor = color; 
                h.FaceColor = 'none';
            end 
        end
        
        %add colorbar, scalebar, set caxis, etc
        colormap('gray')
        axes = gca;
        axes.FontSize = 14;
        caxis(cax);
        colorbar('FontSize',14);
        hold on
        
        %add in scale bar and time ticker
        line([col1, col2], [ysb, ysb], 'Color','w','LineWidth',3);
        hold on
        
        disptime = fr/10;
        text(xsb-100,ysb-20,[num2str(disptime),' s'],'Color', 'w','FontSize',14);
        hold off
        
        %now plot the traces below
        %make subplot for each set of trace data; add vertical bar for each
        %frame value, each plot
        for m = 1:length(mois) 
            if contains(mois{m},'a2')
                ma = 9;
                sub = 'a2';
            elseif contains(mois{m},'a3')
                ma = 11;
                sub = 'a3';
            elseif contains(mois{m},'a4')
                ma = 13;
                sub = 'a4';
            end
            if rem(m,2) == 0
                ma = ma + 1; %10, 12, 14
                subtitle = strcat(sub,' Right');
                colm = 'm';
            else
                subtitle = strcat(sub,' Left');
                colm = 'g';
            end
            subplot(7,2,ma);
            plot(newtracescomb(:,m),colm);
            
            xlabel('Frames')
            xlim([1,frames]);
            ylabel('Ratio')
            ylim([min(newtracescomb(:,m)), max(newtracescomb(:,m))*1.1]);
            title(subtitle)

            %time bar
            xline(ind,'m','LineWidth',2);
            hold on
        end

        %final touches
        %make full frame size and improve text res
        myfigsize = [1,41,1920,963];
        fig.Position = myfigsize;
        set(gcf,'PaperPosition',myfigsize);
        
        %capture each frame for movie
        M(ind) = getframe(gcf);
        close(fig)
        pause(2)
        ind = ind + 1;
    end
    
    %save it all
    vidname5 = strcat(expername,'_',vidt,'_ratiomovie_5fps_gray.avi');
    vid5 = VideoWriter(vidname5,'Uncompressed AVI');
    set(vid5,'FrameRate',5);
    open(vid5);
    writeVideo(vid5,M)
    close(vid5);
end