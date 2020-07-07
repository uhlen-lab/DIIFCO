function  Basic_summary(F_XYZIDint_Space, F_XYZIDint_Space_ThUpper, F_XYZIDint_Space_ThLower, Volume3D, threshold, out_dir)

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
 
 % This script generate summary csv file of 
 % ["Mean intensity"; "Median intensity"; "All cell number"; ...
 % "Positive cell number"; "Negative cell number"; ...
 % "positive cell ratio"; "total density of positive cell(mm3)"; ...
 % "Volume(um3)"; "Threshold"];

mean_intensity = mean(F_XYZIDint_Space(:,5));
median_intensity = median(F_XYZIDint_Space(:,5));

all_cell_number = length(F_XYZIDint_Space);
positive_cell_number = length(F_XYZIDint_Space_ThUpper);
negative_cell_number = length(F_XYZIDint_Space_ThLower);
positve_cell_ratio = (positive_cell_number/ all_cell_number) * 100;
total_density_positive_cell = positive_cell_number ./ Volume3D * (10 ^ 9);

Label_Column = ["Mean intensity"; "Median intensity"; "All cell number"; ...
                "Positive cell number"; "Negative cell number"; ...
                "positive cell ratio"; "total density of positive cell(mm3)"; ...
                "Volume(um3)"; "Threshold"];
            
Data_Column = cat(1, mean_intensity, median_intensity, all_cell_number, ...
                positive_cell_number, negative_cell_number, positve_cell_ratio, ...
                total_density_positive_cell, Volume3D, threshold);

Basic_data_summary = table(Data_Column,'RowNames',Label_Column);

save([out_dir '\' 'Basic_data_summary.mat'],'Basic_data_summary');
writetable(Basic_data_summary, [out_dir '\' 'Basic_data_summary.csv'], 'WriteRowNames',true);

end
