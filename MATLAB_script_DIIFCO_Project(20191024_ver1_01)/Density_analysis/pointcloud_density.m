function [numberOfpoints] = pointcloud_density (pointData, radius, out_dir)

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

    % pointCloud file is a array of three column, X, Y, Z. Table format shows error.
    % radius value used in the publication is 100. Unit is um.
    % I used built-in function of pointCloud in MATLAB below. 
    % [indices,dists] = findNeighborsInRadius(ptCloud,point,radius);

    
%%
p_cloud = pointCloud(pointData);

indices = cell(p_cloud.Count,1);
dists = cell(p_cloud.Count,1);
numberOfpoints = zeros(p_cloud.Count,1);
   
poolobj = gcp('nocreate');
delete(poolobj);
no_of_workers = 12;
parpool ('local',no_of_workers);  
 
parfor d=1:p_cloud.Count
 
     point = pointData(d,:);
     [T_indices,T_dists] = findNeighborsInRadius(p_cloud, point, radius);
     
     indices{d} = T_indices;
     dists{d} = T_dists;
     numberOfpoints(d) = length(T_indices);
     
     disp(['Processing point=' num2str(d)]);

end

pointData_Density=cat(2, pointData, numberOfpoints);
save ([out_dir '\' 'pointData_Density' num2str(radius) 'um.mat'],'pointData_Density');
csvwrite([out_dir '\' 'pointData_Density' num2str(radius) 'um.csv'], pointData_Density);

mean_pointNum = mean(numberOfpoints);
median_pointNum = median(numberOfpoints);
positive_cell_number = length(pointData);
Label_Column = ["Mean density"; "Median density"; "positive cell number"; "radius"];
            
Data_Column = cat(1, mean_pointNum, median_pointNum, positive_cell_number, radius);

Density_data_summary = table(Data_Column,'RowNames',Label_Column);

save([out_dir '\' 'Density_data_summary.mat'],'Density_data_summary');
writetable(Density_data_summary, [out_dir '\' 'Density_data_summary.csv'], 'WriteRowNames',true);

save([out_dir '\' 'indices.mat'],'indices');
save([out_dir '\' 'dists.mat'],'dists');
save([out_dir '\' 'numberOfpoints.mat'],'numberOfpoints');

poolobj = gcp('nocreate');
delete(poolobj);

end