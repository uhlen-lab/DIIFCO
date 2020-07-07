function Intensity2csv(InSitu_mean, Immuno_mean, CC, out_dir)

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

%%

% Get centroids from CC

	s = regionprops(CC,'centroid');
	centroids = cat(1, s.Centroid);
  
% Resolution adjustment
	
	T_centroidsY = centroids(:,1).* 0.585;
	T_centroidsX = centroids(:,2).* 0.585;
	T_centroidsZ = centroids(:,3).* 5;

	centroids_res = cat(2, T_centroidsY, T_centroidsX, T_centroidsZ);

% Preallocation

	%F_XYZ= fix(centroids_res); %Remove floating value
	Temp1=false(length(centroids),1); %make preallocation
	ID=find(not(Temp1)); %get ID list using find zero

% Merge arrays
	
	Merge_mean = cat(2, ID, centroids_res, InSitu_mean, Immuno_mean);

% Save file

	%csvwrite([out_dir '\' 'Merge_mean.csv'],Merge_mean); % This lose precision
	dlmwrite([out_dir '\' 'Merge_mean.csv'],Merge_mean, 'precision', 8, 'delimiter', ',');

	%Make csv file with header
	%This code is from 'CSV with column headers'
	header_string ='ID,X,Y,Z,InSitu,Immuno';
	filename=[out_dir '\' 'Merge_mean_table.csv'];   %'Merge_mean_table.csv'
	fid = fopen(filename,'w');
	fprintf(fid,'%s\r\n',header_string);
	fclose(fid);
	dlmwrite(filename, Merge_mean,'-append','delimiter',',', 'precision', 8);


%Code below does not work
	% Write csv file
	% T = array2table(Merge_mean,...
	% 'VariableNames',{'ID','X','Y','Z', 'InSitu', 'Immuno'}); %Convert array to table

	%writetable(Merge_mean,[out_dir '\' 'Merge_mean_table.csv']);
	%writetable(T,[out_dir '\' 'Merge_mean_table_excel.xlsx']);

