function Points_to_Object_distance_ptCloud(data_dir, image_resolution, PointList_space, out_dir)
 
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
% data_dir is Tiff binary data folder of the object (segmented vessel image)
% PointList has X, Y, Z coordinate (n x 3 array)
% Parallel loop used.
% image_resolution = [0.585 0.585 5];


%% Perimeter list acquisition

Vessel_PointList = Binary_image_to_PointList(data_dir);
x = image_resolution(1);
y = image_resolution(2); 
z = image_resolution(3);
Vessel_PointList_Space = [x * Vessel_PointList(:,1) y * Vessel_PointList(:,2) z * Vessel_PointList(:,3)];
% pointcloud is preconverted to [y x z] because it use regionprops3, you don't need to change axis


%% Preallocation

ptCloud = pointCloud(Vessel_PointList_Space);

point_number = length(PointList_space);
%PC_PointList = pointCloud(PointList_space); % Here PointList should be only [x y z] array without intensity
ind_dists = zeros(point_number,2);


%% Cell to vessel distance calculation

poolobj = gcp('nocreate');
delete(poolobj);
no_of_workers = 12;
parpool ('local',no_of_workers);  

        
    parfor p=1:point_number
            
            [indices, dists] = findNearestNeighbors(ptCloud, PointList_space(p,:),1);
            ind_dists(p,:) = [double(indices), double(dists)];
            disp(['Processing point=' num2str(p)]);
            
    end
        
     
%% Save data

    save([out_dir '\' 'ind_dists.mat'],'ind_dists');
    
    mean_dis = mean(ind_dists(:,2));
    P2O_distance_summary = cat(2, mean_dis, point_number);
        
    Header ={'Data_Mean_distance', 'The_number_of_points'};    
    %csvwrite_with_headers([out_dir '\' 'P2O-distance_summary.csv'], P2O_distance_summary, Header);
    T = array2table(P2O_distance_summary, 'VariableNames', Header);
    writetable(T, [out_dir '\' 'P2O-distance_summary.csv']);
    
    save([out_dir '\' 'P2O_distance_summary.mat'],'P2O_distance_summary');
    
    F_XYZ_Space_NND = cat(2, PointList_space, ind_dists(:,2));
    Header = {'X', 'Y', 'Z', 'NN_distance'};    
    %csvwrite_with_headers([out_dir '\' 'F_XYZ_Space_NND.csv'], F_XYZ_Space_NND, Header);
    T = array2table(F_XYZ_Space_NND, 'VariableNames', Header);
    writetable(T, [out_dir '\' 'F_XYZ_Space_NND.csv']);
    
    csvwrite([out_dir '\' 'Vessel_PointList_Space.csv'], Vessel_PointList_Space)
    save([out_dir '\' 'Vessel_PointList_Space.mat'],'Vessel_PointList_Space');
    
poolobj = gcp('nocreate');
delete(poolobj);       
        
end











