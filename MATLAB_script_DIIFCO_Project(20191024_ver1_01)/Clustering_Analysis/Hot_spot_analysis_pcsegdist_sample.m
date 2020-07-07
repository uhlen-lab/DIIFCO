function Hot_spot_analysis_pcsegdist_sample(positive_cell_list, distance_num, cell_num, Volume3D, out_dir)

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
 
 % positive cell list should have X, Y, Z of real space, not pixel number.
 % distance_num is distance(micrometer) to do clustering(pcsegdist). 20um is used in the publication.
 % cell_num is used to filter out clusters which is small number of cells. 5 cells is used for the publication.
 % Volume_3D is real space volume of 3D tissue sample (um3). Conversion to mm3 is done to this scrpit.
 % out_dir is save directory of analysis

ptCloud = pointCloud(positive_cell_list(:,1:3));

%% pcsegdist
    d=distance_num;

    Cluster_cell_list = double(pcsegdist(ptCloud, d));
    [C, ia, ic] = unique(Cluster_cell_list);
    Cluster_list_count = accumarray(ic,1);
    value_counts = [C, Cluster_list_count];

%% cluster spatial data
    Column_Header = {'X', 'Y', 'Z', 'ID', 'Density', 'Cluster'};
 
    positive_cell_list_pcsedist = cat(2, positive_cell_list, Cluster_cell_list);
    save ([out_dir '\' 'positive_cell_list_pcsedist' num2str(distance_num) 'distance_1cell.mat'],'positive_cell_list_pcsedist');       
    %csvwrite_with_headers([out_dir '\' 'positive_cell_list_pcsedist ' num2str(distance_num) 'distance_1cell.csv'], positive_cell_list_pcsedist, Column_Header);
    T = array2table(positive_cell_list_pcsedist, 'VariableNames', Column_Header);
    writetable(T, [out_dir '\' 'positive_cell_list_pcsedist ' num2str(distance_num) 'distance_1cell.csv']);
    
%% Cell number thresholding
      value_counts_Cth = value_counts(:,2) >= cell_num;
      value_counts_Cth_list = value_counts(value_counts_Cth, :);
              
%% Spatial data of cell threshoulded cluster
      ic_Cth = ismember(ic, value_counts_Cth_list(:,1));
      positive_cell_list_pcsedist_ClasterCellThresholded = positive_cell_list_pcsedist(ic_Cth, :);
      
      Column_Header = {'X', 'Y', 'Z', 'ID', 'Density', 'Cluster'};
      save ([out_dir '\' 'positive_cell_list_pcsedist' num2str(distance_num) 'distance_' num2str(cell_num)  'cell.mat'],'positive_cell_list_pcsedist_ClasterCellThresholded');       
      %csvwrite_with_headers([out_dir '\' 'positive_cell_list_pcsedist ' num2str(distance_num) 'distance_' num2str(cell_num)  'cell.csv'], positive_cell_list_pcsedist_ClasterCellThresholded, Column_Header);
      T = array2table(positive_cell_list_pcsedist_ClasterCellThresholded, 'VariableNames', Column_Header);
      writetable(T, [out_dir '\' 'positive_cell_list_pcsedist ' num2str(distance_num) 'distance_' num2str(cell_num)  'cell.csv']);
    
%% Convert cluster number to 1,2,3 manner
    
      [unique_cluster_num, ia, ic] = unique(positive_cell_list_pcsedist_ClasterCellThresholded(:,6));
      
      positive_cell_list_pcsedist_ClasterCellThresholded_ID = cat(2, positive_cell_list_pcsedist_ClasterCellThresholded, ic);
      
      Column_Header = {'X', 'Y', 'Z', 'ID', 'Density', 'Cluster', 'Cluster_ID'};
      save ([out_dir '\' 'positive_cell_list_pcsedist' num2str(distance_num) 'distance_' num2str(cell_num)  'cell_withClusterID.mat'],'positive_cell_list_pcsedist_ClasterCellThresholded_ID');       
      %csvwrite_with_headers([out_dir '\' 'positive_cell_list_pcsedist ' num2str(distance_num) 'distance_' num2str(cell_num)  'cell_withClusterID.csv'], positive_cell_list_pcsedist_ClasterCellThresholded_ID, Column_Header);
      T = array2table(positive_cell_list_pcsedist_ClasterCellThresholded_ID, 'VariableNames', Column_Header);
      writetable(T, [out_dir '\' 'positive_cell_list_pcsedist ' num2str(distance_num) 'distance_' num2str(cell_num)  'cell_withClusterID.csv']);
    
%%  Density culculation    
      
      % Data of distance_num, 1 cell
      Num_Cluster = length(value_counts(:,1));
      Mean_Cell_Cluster = mean(value_counts(:,2));
      Num_Cell = size(positive_cell_list, 1);
        
      Cluster_density = Num_Cluster ./ Volume3D .* (10 ^ 9); %um3 to mm3
      Num_Cell_density = Num_Cell ./ Volume3D .* (10 ^ 9);
        
      % Data of distance_num, cell_num    
      Num_Cluster_Cth = length(value_counts_Cth_list(:,1));
      Mean_Cell_Cluster_Cth = mean(value_counts_Cth_list(:,2));
      Cluster_density_Cth = Num_Cluster_Cth ./ Volume3D .* (10 ^ 9); %um3 to mm3
      Num_Cell_Cth = size(positive_cell_list_pcsedist_ClasterCellThresholded, 1);
      Num_Cell_density_Cth = Num_Cell_Cth ./ Volume3D .* (10 ^ 9);      
      
           
      % Data summary save
      Label_Column = ["Number_of_cluster"; "Mean_cell_number"; "Cluster_density"; "Number_of_Cell"; "Cell density"; ...
                      "Number_of_cluster_Cell_thresholded"; "Mean_cell_number_Cell_thresholded"; ...
                      "Cluster_density_Cell_thresholded"; "Number_of_Cell_Cell_thresholded"; "Cell density_Cell_thresholded"; ...
                      "Distance_th"; "Cell_th"];            
      Data_Column = cat(1, Num_Cluster, Mean_Cell_Cluster, Cluster_density, Num_Cell, Num_Cell_density, ...
                    Num_Cluster_Cth, Mean_Cell_Cluster_Cth, Cluster_density_Cth, Num_Cell_Cth, Num_Cell_density_Cth, ...
                    distance_num, cell_num);

      Density_data_summary = table(Data_Column,'RowNames',Label_Column);
        
      save([out_dir '\' 'Density_data_summary.mat'],'Density_data_summary');
      writetable(Density_data_summary, [out_dir '\' 'Density_data_summary.csv'], 'WriteRowNames',true); 
 

end