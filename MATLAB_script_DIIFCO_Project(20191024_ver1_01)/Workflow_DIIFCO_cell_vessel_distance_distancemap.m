function Workflow_DIFCO_cell_vessel_distance_distancemap (BW_cell_dir, BW_vessel_dir, BW_inSitu_dir, out_dir)

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
% BW_cell_dir is binary tif file directory which is segmented using nuclear staining
% BW_vessel_dir is binary tif file directory of vessel images which is segmented
% BW_inSitu_dir is a binary tif file directory of In Situ hybridization images which is segmented
% out_dir is the directory you will have data

%% Calculate Centroids

    % calculate centroids of nuclear staining
    out_dir1 = [out_dir '\' 'Nucleus_centroids'];
    mkdir (out_dir1);
    [CC_Nucleus_erosion] = erosion(BW_cell_dir, out_dir1);

    % calculate centroids of in situ staining
    out_dir2 = [out_dir '\' 'In_Situ_centroids'];
    mkdir (out_dir2);
    [CC_inSitu_erosion] = erosion(BW_inSitu_dir, out_dir2);

%% Distance calculation
    out_dir3 = [out_dir '\' 'Distance_analysis'];
    mkdir (out_dir3);
    Merge_distance = points_on_distanceMap (CC_Nucleus_erosion, BW_vessel_dir, out_dir3);
    % Merge_distance is an array of [cell ID, centroidX, centroidY, centroidZ, distance to vessel].

%% Colocalize check

    out_dir4 = [out_dir '\' 'Colocalize_check'];
    mkdir (out_dir4);

    colocarize_check(CC_Nucleus_erosion, CC_inSitu_erosion, Merge_distance, out_dir4)
    % Remove the InSitu positive area which does not overlap with Nuclear staining
    % Remove the InSitu positive signal which overlap with more than two Nuclear cells
    
    % InSitu_distance file is the output of the InSitu positive cells to vessel.
end

