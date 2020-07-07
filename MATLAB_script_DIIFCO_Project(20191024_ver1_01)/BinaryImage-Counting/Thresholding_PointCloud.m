function [F_XYZIDint_Space_ThUpper, F_XYZIDint_Space_ThLower] = Thresholding_PointCloud(F_XYZIDint_Space, threshold, out_dir)
 
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

 % F_XYZIDintSpace is arrary of X(adjusted), Y(adjusted), Z(adjusted), ID, Intensity.
 % Adjustment is used for the resolution of [0.585 0.585 5].
 % This script output the point population which is more than threshold and
 % less than threshold.

 %%
    sele_upper = F_XYZIDint_Space(:,5) >= threshold;
    sele_lower = F_XYZIDint_Space(:,5) < threshold;
    F_XYZIDint_Space_ThUpper= F_XYZIDint_Space(sele_upper,:);
    F_XYZIDint_Space_ThLower= F_XYZIDint_Space(sele_lower,:);
    
    save([out_dir '\' 'F_XYZIDint_Space_ThUpper' num2str(threshold) '.mat'], 'F_XYZIDint_Space_ThUpper');
    save([out_dir '\' 'F_XYZIDint_Space_ThLower' num2str(threshold) '.mat'], 'F_XYZIDint_Space_ThLower');
    
    Header ={'X', 'Y', 'Z', 'ID', 'Intensity'};
    %csvwrite_with_headers([out_dir '\' 'F_XYZIDint_Space_ThUpper' num2str(threshold) '.csv'], F_XYZIDint_Space_ThUpper, Header);
    T = array2table(F_XYZIDint_Space_ThUpper, 'VariableNames',Header);
    writetable(T, [out_dir '\' 'F_XYZIDint_Space_ThUpper' num2str(threshold) '.csv']);
    
    pointData = F_XYZIDint_Space_ThUpper(:,1:3);
    save([out_dir '\' 'pointData_Th' num2str(threshold) ' .mat'], 'pointData');
    
    Header ={'X', 'Y', 'Z'};
    %csvwrite_with_headers([out_dir '\' 'pointData_Th' num2str(threshold) '.csv'], pointData, Header);
    T = array2table(pointData, 'VariableNames',Header);
    writetable(T, [out_dir '\' 'pointData_Th' num2str(threshold) '.csv']);
end