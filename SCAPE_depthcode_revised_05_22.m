%% SCAPE pre-processing code

%Elizabeth Hillman, 8/2021,
%with revisions by 
%Patricia Cooney 1/2022, 
%Wenze Li - 04/2022,
%PC - 05/2022

%this code processes 3D tiffs for R and G channels of SCAPE data

%loads in pre-allocated frames of interest for each video
%loads tiffs of interest
%smooths voxels, bkgd subtract
%choose top-view or side-view
%calculates 3D avg for G & R
%also makes a ratio image from 3D avgs
%saves outputs as tiffs for green and red 16-bit (2D projections of 3D averaged data)
%save output as mat for ratio (floating point numbers)

%%
%load in the excel with frame numbers of interest (col1-name, col2-start,
%col3-stop)

%eois = readtable('framesofinterest-decemberdualcolor_update.xlsx');

%choose top or sideview
%perspective = input('Top or Side? (top = 1, side = 2)');
perspective = 1;

%loop all vids of interest within dir
% for la = 1:height(eois)
%     larvaname = eois.Var1{la};
%     frames = eois.Var2(la):eois.Var3(la); 
% end

%run in batch:
larvanames = {'larva2_run2','larva3_run1','larva4_run2_eL2','larva10_run1_L2','larva10_run4_L2'};
for l = 1:length(larvanames)
    clear ratio
    larvaname = larvanames{l};

    cd 'C:\Users\coone\Desktop\Patricia\newdataset-3D';
    larvafolder = [pwd,'\',larvaname];
    larvaf = dir([larvafolder '/**/','*allsecs']);
    cd(strcat(larvafolder,'\',larvaf.name));

    %troubleshooting options for 3D-Gaussian smoothing of orig data (sf), pixel fluor
    %threshold before including in average (pt), and window size (w)
    sf = 5;
    pt = 0.5;
    w = 15;
    % sf = [3, 5, 7];
    % pt = [0.25, 0.5, 0.75];
    % w = [5, 10, 15];

    frames = 1:601;

    for j = 1:length(frames)
        fnameG = strcat('G_',larvaname,'_',num2str(frames(j)),'.tiff');
        fnameR = strcat('R_',larvaname,'_',num2str(frames(j)),'.tiff');
        Gdata = tiffLoad(fnameG);
        Rdata = tiffLoad(fnameR);

        if perspective == 1 %top-down view
            GDatap = permute(Gdata(:,:,1:end-1),[1, 3, 2]);
            RDatap = permute(Rdata(:,:,1:end-1),[1, 3, 2]);
        elseif perspective == 2
            %keep same dims - sideview
        end

        ss = size(RDatap);

        % smoothing
        GDatasm = smooth3(GDatap,'gaussian',[sf,sf,sf]);
        RDatasm = smooth3(RDatap,'gaussian',[sf,sf,sf]);

        % finding a threshold value over which to keep data - better than
        % processing loads of zeros for the more complicated maximum calcs
        tmpR = max(RDatasm,[],2);
        [hR, vR] = hist(tmpR(:),100);
        threshR2 = vR(2); %takes the max value across pixels and only look at what's above the lowest 2%

        keepR = find(RDatasm>threshR2);
        tmpR = uint8(zeros(size(RDatasm)));
        tmpR(keepR)= 1; % 3D mask of 1's above threshold

        keepyR = find(max(tmpR,[],2)>0); % convert to 2D mask

        % reshape smoothed data for faster processing
        GDatasmL = reshape(permute(GDatasm,[1 3 2]),[ss(1)*ss(3),ss(2)]);
        RDatasmL = reshape(permute(RDatasm,[1 3 2]),[ss(1)*ss(3),ss(2)]);

        %reshape the raw data for the actual averages
        GDatapL = reshape(permute(GDatap,[1 3 2]),[ss(1)*ss(3),ss(2)]);
        RDatapL = reshape(permute(RDatap,[1 3 2]),[ss(1)*ss(3),ss(2)]);

        % prepare variables
        frontsurface = zeros([ss(1),ss(3)]);
        backsurface = zeros([ss(1),ss(3)]);
        thickness = zeros([ss(1),ss(3)]);
        gcampave2 = zeros([ss(1),ss(3)]);
        redave2 = zeros([ss(1),ss(3)]);
        maxpos = zeros([ss(1),ss(3)]);
        L = size(RDatasm,2);
        pixthresh = pt; % average adjacent pixels that are > % of the max.
        win = w; % using a window of +/- w around the max (within which we take the top 50% pixels)

        %loops through the super-threshold pixels and makes new 2D
        %images with dynamic z-window size
        for i = 1:length(keepyR)
            mm = min(RDatasmL(keepyR(i),:));
            [p pp] = max(RDatasmL(keepyR(i),:)-mm,[],2);
            maxpos(keepyR(i)) = pp; %finding the 3D max pixel location in 2D space
            winnys = max([1 pp-win]); %so comparing window with 1 and L allows you to deal with borders of the image
            winnye = min([pp+win, L]);
            use2 = winnys-1+find(RDatasmL(keepyR(i),winnys:winnye)-mm>pixthresh*p);
            thickness(keepyR(i)) = length(use2);
            frontsurface(keepyR(i)) = winnys;
            backsurface(keepyR(i)) = winnye;
            gcampave2(keepyR(i)) = mean(GDatapL(keepyR(i),use2)); %here, taking mean across x and z values
            redave2(keepyR(i)) = mean(RDatapL(keepyR(i),use2));
        end


        pause(0.1);


        %% resize and rotate images so larva is proportional
        imsz = size(gcampave2);
        gcresz = imresize(gcampave2,[imsz(1)*1.0999,imsz(2)*4]);
        gcreszrot = rot90(gcresz,-1);
        rcresz = imresize(redave2,[imsz(1)*1.0999,imsz(2)*4]);
        rcreszrot = rot90(rcresz,-1);

        ratio(:,:,j) = ((gcreszrot)-80./rcreszrot);
        
        %%
        %save gcamp depth-averaged 2D
        gcampsave = uint16(gcreszrot);
        filename = strcat('fullrun_resize_rotate_pxthresh',num2str(pt),'_smooth',num2str(sf),'_win',num2str(w),'_imageunsmooth_top_R_based_depth_gcamp_',larvaname);
        imwrite(gcampsave, strcat(filename, '.tiff'), 'Tiff', 'Compression', 'none', 'WriteMode', 'Append');
        %save mcherry depth-averaged 2D
        redsave = uint16(rcreszrot);
        filename = strcat('fullrun_resize_rotate_pxthresh',num2str(pt),'_smooth',num2str(sf),'_win',num2str(w),'_imageunsmooth_top_R_based_depth_mcherry_',larvaname);
        imwrite(redsave, strcat(filename, '.tiff'), 'Tiff', 'Compression', 'none', 'WriteMode', 'Append');

    end
    %save ratio
    ratio(ratio==Inf) = 0;
    ratio(ratio==-Inf) = 0;
    save(strcat('fullrun_',larvaname,'_ratio.mat'),'ratio','-v7.3');
end
