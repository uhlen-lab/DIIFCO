function Merge_distance = points_on_distanceMap (CC_erosion, BW_vessel_dir, out_dir)

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
% BW_cell_dir is a folder of binary images of segmented cells which was eroded already.
% BW_vessel_dir is folder of segmented vessel images for distance map(tiff)
% Cell and vessel binary images should have the equivalent resolution for each dimension
% (this case 5x5x5um is assumed)

%% Getting points cloud and locate them in the image cordinate.

% Component analysis and get centroid information
SS = regionprops(CC_erosion,'Centroid');
centroids = cat(1, SS.Centroid);
 
% Resolution adjustment
T_centroidsY = centroids(:,1).* 0.585;
T_centroidsX = centroids(:,2).* 0.585;
T_centroidsZ = centroids(:,3).* 5;

centroids_res = cat(2, T_centroidsY, T_centroidsX, T_centroidsZ); % Store resolution adjusted coordinates

% Preallocation
Temp1=false(length(centroids),1); %make preallocation
ID=find(not(Temp1)); %get ID list using find zero

% Merge arrays
centroids_res_ID = cat(2, ID, centroids_res);

%% Open vessel image data

[vessel_BW] = ImageSeries_openBinary(BW_vessel_dir);

%% Generation of distance map

D1 = bwdist(vessel_BW,'euclidean');
D1_res5 = D1*5;

Distance = zeros(length(ID),1);

parfor d=1:length(ID)
    
    point = centroids_res(d,:)/5; % get point coordinate
    current_pixel = floor(point); %Round down to get pixel coordinate
    nearest_pixel = ceil(point); % Round up to get nearest pixel to the point
    
    if nearest_pixel(1) > size(D1,1) nearest_pixel(1)=current_pixel(1); end   %check if value of nearest pixel is out of range
    if nearest_pixel(2) > size(D1,2) nearest_pixel(2)=current_pixel(2); end
    if nearest_pixel(3) > size(D1,3) nearest_pixel(3)=current_pixel(3); end
    if current_pixel(1) == 0         current_pixel(1)=nearest_pixel(1); end   %check if value of current is 0
    if current_pixel(2) == 0         current_pixel(2)=nearest_pixel(2); end
    if current_pixel(3) == 0         current_pixel(3)=nearest_pixel(3); end
    
    current_pixel_ind = sub2ind(size(D1), current_pixel(2), current_pixel(1), current_pixel(3)); % index transform
    nearest_pixel_ind = sub2ind(size(D1), nearest_pixel(2), nearest_pixel(1), nearest_pixel(3)); % index transform
    
    dist_C = D1_res5(current_pixel_ind); % get actual distance value of current pixel from distance map*5um file
    dist_N = D1_res5(nearest_pixel_ind); % get actual distance value of nearest pixel from distance map*5um file
        
    current_nearest = [current_pixel ; nearest_pixel]; % make matrix to calculate ratio of current/nearest pixel
    cn_distance = pdist(current_nearest); % calculate the distance between current/nearest pixel
    current_point = (point-current_pixel).*(nearest_pixel-current_pixel);           % get vector which have value only along the axis of current and nearest pixel 
    cp_distance = pdist([0 0 0 ; current_point]);    % calculate the distance of point along the axis of current and nearest_pixel
    
    if current_pixel == nearest_pixel     % check if current and nearest pixel is the same (this case cn_distance will be inf, then error)
        dist_p = dist_C;
    else
        
    ratio = cp_distance/cn_distance;      % calculate length ratio of current to point and current to nearest
    dist_p = ratio .* dist_N + (1-ratio) .* dist_C;  % calculate the distance of the point using ratio
    end
    
    Distance(d, 1) = dist_p;
    disp(['Object number=' num2str(d) 'is done!']);
    
end

%% Save data

% Merge arrays
Merge_distance = cat(2, ID, centroids_res, Distance);
mkdir([out_dir '\distance_map_tif']);

parfor f=1:size(D1_res5,3)
    
    disp(['Saving distance map image' num2str(f) '...']); 
    image_tif=uint16(D1_res5(:,:,f));
    fname_out =[out_dir '\distance_map_tif\distance_map_5um' num2str(f,'%04i') '.tif'];
    imwrite(image_tif, fname_out);
    
end

disp(['Saving distance result....']);
   
    save([out_dir '\' 'Merge_distance.mat'],'Merge_distance');
    save([out_dir '\' 'centroids.mat'],'centroids');
    save([out_dir '\' 'centroids_res.mat'],'centroids_res');
    save([out_dir '\' 'distance_map.mat'],'D1');
    save([out_dir '\' 'distance_map_5um.mat'],'D1_res5');

%csvwrite([out_dir '\' 'Merge_mean.csv'],Merge_mean); % This lose precision
dlmwrite([out_dir '\' 'Merge_distance.csv'],Merge_distance, 'precision', 8, 'delimiter', ',');
    
header_string ='ID,X,Y,Z,Distance';
filename=[out_dir '\' 'Merge_distance_table.csv'];   %'Merge_mean_table.csv'
fid = fopen(filename,'w');
fprintf(fid,'%s\r\n',header_string);
fclose(fid);
dlmwrite(filename, Merge_distance,'-append','delimiter',',', 'precision', 8);
        
    disp('Done!!!');   

end
