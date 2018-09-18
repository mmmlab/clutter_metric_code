%% CLUTTER METRIC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Below is a MATLAB implementation of the clutter metric described in Semizer and Michel (2018). 
% This metric is a modified version of the clutter metric described in Bravo and Farid (2008).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% STEPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 1. Define all the parameters
% 2. Compute number of regions in each image in the data set 
% 3. Compute the distribution of exponents over the image data set
% 4. Compute the average log slope
% 5. Compute the clutter distribution of images
% 6. Compute clutter ranks
% 7. Save the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [nr_regions,clutter_distro,clutter_ranks,seg_images,fileNames] = getClutter(image_dir)
%% STEP 1: Define all the parameters here
dirData = dir(sprintf('%s/*.png',image_dir));           %data from the corresponding image directory
fileNames = {dirData.name};                             %names of the image files
pl_crits = 2.^(6.5:0.5:12);                             %scale of segmentation
nr_crit = length(pl_crits);                             %number of scales of segmentation
nr_images = numel(fileNames);                           %number of images
nr_regions = zeros(nr_crit,nr_images);                  %array to contain number of regions in segmented images
seg_images = cell(1,nr_images);                         %cell to contain segmented images
raw_clutter_distro = zeros(1,nr_images);                %array to contain distribution of raw clutter values
clutter_distro = zeros(1,nr_images);                    %array to contain distribution of clutter values
exponent_distro = zeros(1,nr_images);                   %array to contain distribution of exponent values
nr_samples = 8;                                         %number of random samples to generate sections from each image
img_samp_size = 1024;                                   %size of each random sample
min_pix_size = 2^floor(log2(1e-4*pi*512^2));            %minimum size of a region within an image

%% STEP 2: Compute number of regions in each image in the data set 
for i = 1:nr_images
    i
    tic
    file_name = sprintf('%s/%s',image_dir,fileNames{i}) 
    full_img = imread(file_name); %get image
    [ht,wd,dim] = size(full_img); %get size of the image
    if(dim~=3) %check the dim of the image
        full_img = cat(3,full_img,full_img,full_img); %create three channels
        sprintf('...You have a grayscale image...')
    end
    
    %Sample the image
    samp_imgs = cell(1,nr_samples);
    for samp = 1:nr_samples
            h_start = randi([0,max(ht-img_samp_size,0)])+1;
            w_start = randi([0,max(wd-img_samp_size,0)])+1;
            h_end = min(h_start+img_samp_size-1,ht);
            w_end = min(w_start+img_samp_size-1,wd);
            samp_imgs{samp} = full_img(h_start:h_end,w_start:w_end,:);
    end
    
    %Compute segmentation of the image for each sample 
    for l = 1:nr_crit
        pl_crit = pl_crits(l);
        sample_regions=zeros(1,nr_samples);
        for samp = 1:nr_samples
            img = samp_imgs{samp};
            [~, noRegions]= segmentFelzenszwalb(img, 0.5, pl_crit, min_pix_size, false, 0 ,false);
            sample_regions(samp) = max(noRegions,1);
        end
        %Compute the (geometric) mean of the sample results
        nr_regions(l,i) = prod(sample_regions)^(1/nr_samples);
    end
    %Compute a segmentation of the entire image (for evaluation purposes)
    [segResult, ~]= segmentFelzenszwalb(full_img, 0.5, 2048, min_pix_size, true, 0 , true);
    seg_images{i} = segResult;
    toc
end

%% STEP 3: Compute the distribution of exponents over the image data set
for i = 1:nr_images
    y = log(nr_regions(:,i));
    k = log(pl_crits');
    % find the least-squares log-linear fit
    coeff = polyfit(k,y,1);
    log_slope = coeff(1);
    log_intercept = coeff(2);
    raw_clutter_distro(i) = exp(log_intercept);
    exponent_distro(i) = exp(log_slope);
end

%% STEP 4: Compute the average log slope
avg_log_slope = mean(log(exponent_distro));

%% STEP 5: Compute the clutter distribution of images
for i = 1:nr_images
    y = log(nr_regions(:,i));
    k = log(pl_crits');
    obj_func = @(c) sum((k.*avg_log_slope+c-y).^2);
    log_clutter = fminsearch(obj_func,log(raw_clutter_distro(i)));
    clutter_distro(i) = round(exp(log_clutter));
end

%% STEP 6: Compute clutter ranks
indices = 1:nr_images;
temp = sortrows([indices(:),clutter_distro(:)],2);
clutter_ranks = temp(:,1);

%% STEP 7: Save the results
save('clutter_metric_data.mat','nr_regions','clutter_distro','clutter_ranks','seg_images','fileNames');

end