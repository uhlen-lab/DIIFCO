function mean_NN_MonteCarlo_3D_pdist(area_dir, F_XYZIDint, threshold, trial_number, out_dir)

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
 % This script perform MonteCarlo simulation for distribution of cells in the space defined by area_dir. 
 % The number of voxel of spaces are randomly selected and average nearest distance is obtained  
 % The number of trial (trial_number) is how many times you will do simulation and get average value.
 % Threshold is used to extract points more than (>) thresold intensity
 % Nearest neibor index(NNI) is to show how dense or sparse cells are distributed.
 % If NNI is less than 1, cells distribute densly more than simulation
 % If NNI is more than 1, cells distribute sparsely more than simulation
 
%% preallocation

area_image = ImageSeries_openBinary(area_dir);
indM = find(area_image == 1);

Thre = F_XYZIDint(:,5) > threshold;
F_XYZIDint_Threshold = F_XYZIDint(Thre,:);
point_number = length(F_XYZIDint_Threshold);

emNN_mean = zeros(trial_number,1);
sz=size(area_image); %Image_size=[y,x,z]; %pixel value

ind_dists_All = zeros(point_number,2, trial_number);
random_points_space_All = zeros(point_number, 3, trial_number);

%% Monte Carlo simulation

poolobj = gcp('nocreate');
delete(poolobj);
no_of_workers = 12;
parpool ('local',no_of_workers);  

for t=1:trial_number

        ind_dists = zeros(point_number,2);
        random_points = datasample(indM, point_number);
        [y, x, z] = ind2sub(sz, random_points);
        random_points_space = [0.585*x 0.585*y 5*z];
       
        parfor p=1:point_number
            
            [dists, indices] = pdist2(random_points_space, random_points_space(p,:), 'euclidean', 'Smallest' ,2);
            ind_dists(p,:) = [dists(2), indices(2)];
            %disp(['Processing point=' num2str(p)]);
            
        end
        
        emNN_mean(t) = mean(ind_dists(:,1)); 
        ind_dists_All(:,:,t) = ind_dists;
        random_points_space_All(:,:,t) = random_points_space;
        
        disp(['Simulation trial=' num2str(t) ' done.']);
end


%% Calculation of the mean NN value of actual data

        data_ind_dists = zeros(point_number,2);
        y = F_XYZIDint_Threshold(:,1);
        x = F_XYZIDint_Threshold(:,2);
        z = F_XYZIDint_Threshold(:,3);
        data_points_space = [0.585*x 0.585*y 5*z]; % Resolution is optimized for the COLM system [0.585*x 0.585*y 5*z]
     
        parfor p=1:point_number
            
            [dists, indices] = pdist2(data_points_space, data_points_space(p,:), 'euclidean', 'Smallest' ,2);
            data_ind_dists(p,:) = [dists(2), indices(2)];
            disp(['Processing point=' num2str(p)]);
            
        end
        
        data_mean = mean(data_ind_dists(:,1));

%% Calculation of NN index       
   

    emNN_random_mean = mean(emNN_mean);
    emNN_random_std = std(emNN_mean);
    %enNN_random_stE = emNN_random_std / sqrt(point_number);

    [h,p,ci,zval] = ztest(data_ind_dists(:,1), emNN_random_mean, emNN_random_std);
    
    NN_index = data_mean/ emNN_random_mean;
    
    
%% Save data

        save ([out_dir '\' 'data_ind_dists.mat'],'data_ind_dists');       
        NNA_summary = cat(2, data_mean, emNN_random_mean, emNN_random_std, ... 
                                  point_number, trial_number, NN_index, ...
                                  h, p, ci(1), ci(2), zval);
        Header ={'Data_Mean_distance', 'Simulation_Mean_distance', 'Simulation_Standard_deviation', ...
                 'Point_number', 'Trial_number', 'NN_index', ...
                 'Hypothesis', 'p-value', 'Confidence interval1', 'Confidence interval2','z-value'};
        %csvwrite_with_headers([out_dir '\' 'NNA_summary.csv'], NNA_summary, Header);
        T = array2table(NNA_summary, 'VariableNames',Header);
        writetable(T, [out_dir '\' 'NNA_summary.csv']);
               
        save([out_dir '\' 'NNA_summary.mat'],'NNA_summary');

        save ([out_dir '\' 'ind_dists_All.mat'],'ind_dists_All');
        save ([out_dir '\' 'emNN_mean.mat'],'emNN_mean');
        save ([out_dir '\' 'emNN_random_mean.mat'],'emNN_random_mean');
        save ([out_dir '\' 'emNN_random_std.mat'],'emNN_random_std');
        
        save([out_dir '\' 'random_points_space_All.mat'],'random_points_space_All');
        
        random_points_space1_forcheck = random_points_space_All(:,:,1);
        csvwrite([out_dir '\' 'random_points_space1_forcheck.csv'], random_points_space1_forcheck);
        
poolobj = gcp('nocreate');
delete(poolobj);

end