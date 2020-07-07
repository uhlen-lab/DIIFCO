function PointList = Binary_image_to_PointList(data_dir)

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

%% Open binary image

BW = ImageSeries_openBinary(data_dir); %Open binary image
BW2 = bwperim(BW,26);
image_size= size(BW2);
Voxel_list_ind = find(BW2 == 1);

[Y, X, Z] = ind2sub(image_size,Voxel_list_ind);

PointList = cat(2, X, Y, Z); %Convert the CC.PixelList to Array(PointList)

end