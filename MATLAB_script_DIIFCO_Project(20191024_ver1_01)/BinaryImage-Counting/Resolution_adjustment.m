function F_XYZIDint_Space = Resolution_adjustment(F_XYZIDint, out_dir)

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
 % Adjustment of pixel to actual volume
 % [0.585 0.585 5] is resolution specific to the COLM light sheet system
 % F_XYZIDint is arrary of X, Y, Z, ID, Intensity.
 
    F_XYZIDint_Space = F_XYZIDint;
    F_XYZIDint_Space(:,3) = F_XYZIDint_Space(:,3)*5;
    F_XYZIDint_Space(:,2) = F_XYZIDint_Space(:,2)*0.585;
    F_XYZIDint_Space(:,1) = F_XYZIDint_Space(:,1)*0.585;
    
    save([out_dir '\' 'F_XYZIDint_Space.mat'],'F_XYZIDint_Space');
    
    Header ={'X', 'Y', 'Z', 'ID', 'Intensity'};    
    T = array2table(F_XYZIDint_Space, 'VariableNames',Header);
    writetable(T, [out_dir '\' 'F_XYZIDint_Space.csv']);
    
    pointData = F_XYZIDint_Space(:,1:3);
    save([out_dir '\' 'pointData.mat'],'pointData');
    
end