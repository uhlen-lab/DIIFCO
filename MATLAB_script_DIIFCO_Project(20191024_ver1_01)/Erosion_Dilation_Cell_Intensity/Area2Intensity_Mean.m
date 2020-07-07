function Intensity_mean = Area2Intensity_Mean(data_dir, PixelArea, out_dir)

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

 % data_dir is direct path to folder which have tif series of intensity data
 % out_dir is direct path to save output
 % PixelArea is cell array which have object pixel area list in each cell. Cell(number, 1)

%% Get image size
	
	Im=datastore(data_dir, 'FileExtensions', '.tif','Type', 'image');%image data
	list = dir([data_dir '\' '*.tif']); % get image file names and assign to variable 'list'
	First_image=logical(readimage(Im,1));
	sz=size(First_image);
	x=sz(2);
	y=sz(1);
	z=length(list);
	Image_size=[y,x,z];

%% read intensity images

	Intensity_Image = zeros(y, x, z);

	poolobj = gcp('nocreate');
	delete(poolobj);
	no_of_workers = 12;
	parpool ('local',no_of_workers);  
 
	parfor d=1:z
 
     		T_image=readimage(Im,d);   
     		Intensity_Image(:,:,d)= T_image;
  
	end
 
	poolobj = gcp('nocreate');
	delete(poolobj);

%% Read object area
    % Use PixelArea as it is if you have file, or use datastore
    % ds= datastore(PixelArea_dir);
    % PixelArea=tall(ds);
	
	Object_Number = length(PixelArea);

%% Preallocation

	Intensity_mean = zeros(length(PixelArea),1);    % To store mean intensity value for each object

%% Mean calculation

%   Sometimes parfor does not work. In this case use for.

% 	poolobj = gcp('nocreate'); %parallel pool for parfor
% 	delete(poolobj);
%  	no_of_workers = 12;
%  	parpool ('local',no_of_workers);  
%   parfor obj=1:Object_Number   
	
    for obj=1:Object_Number   
    
		disp(['Object number=' num2str(obj)]);
    
            % Object pixel retrieval
    		Pixel_object=PixelArea(obj);                    % Get list of pixel indext of objects
    		Matrix3D=cell2mat(Pixel_object);                % Conver the index value of pixels to array
    		Area_Intensity = Intensity_Image(Matrix3D);     % Get pixel value of Intensity_Images
    		Intensity_mean(obj,1) = mean2(Area_Intensity);
    
    end     
   
% 	poolobj = gcp('nocreate');
%  	delete(poolobj);
 
    dir_name = inputname(1);
	save([out_dir '\' 'Intensity_mean_' dir_name '.mat'],'Intensity_mean');
 
end 
 