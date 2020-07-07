function [CC] = erosion(data_dir, out_dir)

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
 % This script aim to perform erosion operation for big data. 
 % data_dir has segmented tif image (nuclear staining).
 % erosion shape used in the paper is strel('diamond',3), but you can change
 % processed data(CC = bwconncomp(BWd,26)) is saved in out_dir

%% Open image data
Im=datastore(data_dir, 'FileExtensions', '.tif','Type', 'image');
list = dir([data_dir '\' '*.tif']); % get image file names and assign to variable 'list'
First_image=logical(readimage(Im,1));
sz=size(First_image);
x=sz(2);
y=sz(1);
z=length(list);

BW = false(y, x, z);  
 
poolobj = gcp('nocreate');
delete(poolobj);
no_of_workers = 12;
parpool ('local',no_of_workers);  
 
parfor d=1:length(list)
 
     T_image=logical(readimage(Im,d));
     BW(:,:,d)= T_image;
 
end
 
%% Erosion preallocation

    shape_name = 'diamond';
    se_size = 3;

    SE_shape = strel(shape_name, se_size);
    se = SE_shape.Neighborhood;
    BWd = false(y, x, z);


%% Erosion
parfor d=1:length(list)
 
     T_image=BW(:,:,d);
     T_image = imerode(T_image,se);
     BWd(:,:,d)= T_image;
     disp(['Slice' num2str(d) ' Processed']);

end

%% Connected component analysis
    disp(['Analyzing connected components....']); 
    CC = bwconncomp(BWd,26);

%% Save data

out_dir1 = [out_dir '\' 'erosion'];
mkdir (out_dir1);

parfor d=1:length(list)
    
    disp(['Saving ' num2str(d) '...']); 
    image_tif=BWd(:,:,d);
    fname_out =[out_dir1 '\erosion_' shape_name num2str(se_size) '_' num2str(d,'%04i') '.tif'];
    imwrite(image_tif, fname_out);
    
end
    disp(['Saving CC....']); 
    save([out_dir '\' 'CC_erosion_' shape_name num2str(se_size) '.mat'],'CC');
    disp('Done!!!');   

end