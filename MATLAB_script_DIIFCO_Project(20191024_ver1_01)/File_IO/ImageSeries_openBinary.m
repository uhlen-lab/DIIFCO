function [BW] = ImageSeries_openBinary(data_dir)

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
% Karolinska Institutet, Stockholm, Sweden 

%% Summary
% This script is to read tiff file from directory specified and store as
% logical array.

 Im=datastore(data_dir, 'FileExtensions', '.tif','Type', 'image');
 list = dir([data_dir '\' '*.tif']); 
 Size_image=logical(readimage(Im,1));
 [y, x] = size(Size_image);
 
 BW = false(y, x, length(list));
 
 poolobj = gcp('nocreate');
 delete(poolobj);
 no_of_workers = 12;
 parpool ('local',no_of_workers);  
 
 parfor d=1:length(list)
 
     T_image=logical(readimage(Im,d));
     BW(:,:,d)= T_image;
 
 end
 
 poolobj = gcp('nocreate');
 delete(poolobj);
 
end