function Workflow_DIFCO_cell_vessel_distance (nucleus_data_dir, Int_data_dir, threshold, area_dir, vessel_dir, out_dir)

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

%% Summary
% nucleus_data is binary tif file directory which is going to be counted
% Int_data_dir is 16 bit tif file directory (usually in situ hybridization data)
% threshold is the value you will cut-off in Int_data_dir
% area_dir is the binary tif file directory which is segmented for whole sample shape.
% out_dir is the directory you will have data

%% Counting cells and basic statistic 
out_dir1 = [out_dir '\' 'Cell_Segment_Intensity'];
mkdir (out_dir1);
% Save data for counting data and summary

F_XYZID = Hmaxima2PointCloud_3D(nucleus_data_dir, out_dir1); 
% Count segmented cells (centroids) using bwconncomp and regionprops3

SE_shape= strel('disk', 8,0);
se=SE_shape.Neighborhood; 
%se is 2 dimension only. You can change if you want.
%This parameter is for the COLM system, 10X lens.

F_XYZIDint = Point3Dto2D_Dilate_MeanInt_Par (Int_data_dir, F_XYZID, se, out_dir1);
% Generate circle area defined with se (above) in each cell and measure mean intensity in Int_data_dir.

F_XYZIDint_Space = Resolution_adjustment(F_XYZIDint, out_dir1);
% Adjustment of pixel to actual volume
% [0.585 0.585 5] is resolution specific to the COLM light sheet system 

[F_XYZIDint_Space_ThUpper, F_XYZIDint_Space_ThLower] = Thresholding_PointCloud(F_XYZIDint_Space, threshold, out_dir1);
% This function outputs the point population which is more than threshold and less than threshold.
 
area_image = ImageSeries_openBinary(area_dir);
Volume3D = length(find(area_image == 1)) * 0.585 * 0.585 * 5; 
%[0.585 0.585 5] is resolution specific to the COLM light sheet system 
%Volume of space used for analysis is calculated.

Basic_summary(F_XYZIDint_Space, F_XYZIDint_Space_ThUpper, F_XYZIDint_Space_ThLower, Volume3D, threshold, out_dir1);
% This save the summary in out_dir. 
% This script generate summary csv file of 
% ["Mean intensity"; "Median intensity"; "All cell number"; ...
% "Positive cell number"; "Negative cell number"; ...
% "positive cell ratio"; "total density of positive cell(mm3)"; ...
% "Volume(um3)"; "Threshold"];

%% Cell to vessel distance calculation

out_dir2 = [out_dir '\' 'Cell_to_Vessel_Analysis']; % Save directory
mkdir (out_dir2);
image_resolution = [0.585 0.585 5]; % The resolution of the COLM system

PointList_space = F_XYZIDint_Space_ThUpper(:, 1:3); % positive cell list in real space axis

Points_to_Object_distance_ptCloud(vessel_dir, image_resolution, PointList_space, out_dir2);

% vessel_dir should have binary tif images (segmented vessel images)
% PointList_space is X, Y, Z of real space array, not pixel.
% image resolution is for the COLM system  [0.585 0.585 5]
% out_dir2 is to save data in the directory


end