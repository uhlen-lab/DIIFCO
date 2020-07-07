function [Dilate_Pixel_Object] = AreaIndDilation(data_dir, out_dir, CC, se)

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

    % This script dilate objects defined with CC using se structure element.
    % The structure element used in the paper is below.
	% SE_shape = strel('diamond',3);
	% se=SE_shape.Neighborhood;

    % This complex scripting is for using imdilate for processing
    % I think entire script can be more simplified.

%% Get image size
    Im=datastore(data_dir, 'FileExtensions', '.tif','Type', 'image');%image data for chennel 1
    list = dir([data_dir '\' '*.tif']); % get image file names and assign to variable 'list'
    First_image=logical(readimage(Im,1));
    sz=size(First_image);
    x=sz(2);
    y=sz(1);
    z=length(list);
    Image_size=[y,x,z];

%% Get object number from CC file (from binary image)
    
    Pixel_object_All = CC.PixelIdxList;
    Object_Number = CC.NumObjects;


%% Preallocation
    
    Dilate_Pixel_Object = cell(CC.NumObjects,1);
    %Dilate_Pixel_Object = tall(cell(CC.NumObjects,1));

%% Make diluted objects
    
    poolobj = gcp('nocreate'); %parallel pool for parfor
    delete(poolobj);
    no_of_workers = 12;
    parpool ('local',no_of_workers);  
 
    parfor obj=1:Object_Number   
    
        disp(['Object number=' num2str(obj)]);
    
        % Object pixel retrieval
        Pixel_object = Pixel_object_All(obj);               % Get list of pixel indext of object(obj)
        Matrix3D = cell2mat(Pixel_object);                  % Convert the pixel_object to an array
        [indY, indX, indZ] = ind2sub(Image_size, Matrix3D); % Convert index to y, x, z matrix
    
    %% Make temporary image box to do imerode
        
        % Make cell object box for dilation
        ImageBoxY = max(indY)-min(indY)+1 + 6;  % Y-axis size. +6 is for 'diamond' 3 for up and down.
        ImageBoxX = max(indX)-min(indX)+1 + 6;  % X-axis size. +6 is for 'diamond' 3 for left and right.
        ImageBoxZ = max(indZ)-min(indZ)+1;  % Z-axis size. No 3D dilation.
        ImageBoxSize = [ImageBoxY,ImageBoxX,ImageBoxZ];     % Box size
        T_array = zeros(ImageBoxY, ImageBoxX, ImageBoxZ);   % Temporary imaze array only contains one cell
    
        % Make the each cell [Y X Z] minimum and fit the Imagebox
        T_array_Y = indY-min(indY)+1 +3;    % Make the Y,X,Z coordinate minimum, but shifted +3 right and down
        T_array_X = indX-min(indX)+1 +3;    % For dilation by 'diamond', 3, extra 3 dots are added
        T_array_Z = indZ-min(indZ)+1;       % No 3D dilation
    
        linearInd_T_array = sub2ind(ImageBoxSize, T_array_Y, T_array_X, T_array_Z); % Convert y,x,z to linear index
        T_array(linearInd_T_array) = 1 ;    % Linear index value is 1.
    
        T_array = imdilate(T_array, se);    % This is what we want do here. Dilate temporal each cell image
                                            % This complex scripting is for using imdilate for processing
                                            % I think entire script can be more simplified.
    
        Dilate_ind = find(T_array == 1);    % Get dilated linear index of pixel =1
        [DindY_T, DindX_T, DindZ_T] = ind2sub(ImageBoxSize, Dilate_ind);    % Convert linear index to matrix
    
        DindY = DindY_T + min(indY)-4;      % Convert temporal box to real value by adding min(ind) and -4
        DindX = DindX_T + min(indX)-4;
        DindZ = DindZ_T + min(indZ)-1;
    
        %% Check the value is <0 or >max size
        
        % pixels located edge of images will be removed
        DindY_C = (DindY<= y) & (DindY >= 1); % width check (Get logical array) (select the pixel =>1 and <= image max)
        DindX_C = (DindX<= x) & (DindX >= 1); % height check 
        Dind_XY= DindY_C & DindX_C; % Get logical array which has 0 in Y and X
        DindY(Dind_XY==0) = []; % remove out of range pixels. If the value == 0, then remove from matrix.
        DindX(Dind_XY==0) = []; 
        DindZ(Dind_XY==0) = []; 
       
        Dilate_Pixel_Object{obj,1} = sub2ind(Image_size, DindY, DindX, DindZ);    % Store the dilated cell pixels
         
    end
    
 %% Saving data
 
    disp('Saving...');   
    save([out_dir '\' 'Dilate_Pixel_Object.mat'],'Dilate_Pixel_Object');
    
    disp('Done!!!');   

    poolobj = gcp('nocreate');
    delete(poolobj);
    
end