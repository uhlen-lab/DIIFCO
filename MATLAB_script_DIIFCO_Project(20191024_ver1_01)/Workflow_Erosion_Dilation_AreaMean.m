function Workflow_Erosion_Dilation_AreaMean(nucleus_dir, Immuno_dir, InSitu_dir, out_dir)

% Scripts for DIIFCO(for diagnosing in situ immunofluorescence-labelled
% cleared oncosamples) project
%
% by Per Uhlen lab 
% Ver 1.01
%
% These scripts are written by Shigeaki Kanatani 
% 
% Contact: Per Uhlén, per.uhlen@ki.se
%          Shigeaki Kanatani, shigeaki.kanatani@ki.se
% 
%          Karolinska Institutet, Stockholm, Sweden 

%% Heading

% This is workflow script for 3d cell counting using erosion-dilation technique
% nucleus data is segmented nucleus directory
% Immuno_dir and InSitu_dir is the directories of intensity data
% out_dir is save directory
% functions used in this script is E:Erosion_Dilation_Cell_Intensity folder


%% Define structure element

	SE_shape = strel('diamond',3); % This structure element was used in the paper
	se = SE_shape.Neighborhood;

%% Erosion here

    % In the paper, we did erosion in Amira. 
    [BW] = ImageSeries_openBinary(nucleus_dir);
    CC_erosion = bwconncomp(BW,26);
    save([out_dir '\' 'CC_erosion.mat'],'CC_erosion'); clear BW;
     
    % The matlab script for erosion is below.    
	% CC_erosion = erosion(nucleus_dir, out_dir);
	% after erosion, get CC_erosion and save image
	
    
%% Dilation here
	
	[Dilate_Pixel_Object] = AreaIndDilation(nucleus_dir, out_dir, CC_erosion, se);
	% Dilate_Pixel_Object is a cell array which has dilated cell pixels in each object
	% after dilation anad save image(for checking)
	
%% Get mean intensity value from dilated area

	Immuno_mean = Area2Intensity_Mean(Immuno_dir, Dilate_Pixel_Object, out_dir);
	InSitu_mean = Area2Intensity_Mean(InSitu_dir, Dilate_Pixel_Object, out_dir);

	% Get mean intensity data for Immuno and InSitu data sets
	% Dilate_Pixel_Object is dilated pixel are from CC_erosion 

    % You can use parfor for this script, but large memory will be required.
    % Immuno_mean = Area2Intensity_Mean_parfor(Immuno_dir, Dilate_Pixel_Object, out_dir);
    % InSitu_mean = Area2Intensity_Mean_parfor(InSitu_dir, Dilate_Pixel_Object, out_dir);
    
%% Calculate centroids(pixel) and change relolution from pixel to real um value and save data

	Intensity2csv(InSitu_mean, Immuno_mean, CC_erosion, out_dir);
	% resolution is [0.585 0.585 5] for the COLM system
	% The save data is 'Merge_mean.csv', which has [ID, ceintroids[X Y Z], InSitu_mean, Immuno_mean]

end


