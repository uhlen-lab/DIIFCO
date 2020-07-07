function colocarize_check(CC_Nucleus, CC_inSitu, Merge_distance, out_dir)

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

%% Generate matrix from cell array

    % This section make a list of pixel and cell ID for both CC1(Nucleus-1) and CC2 (In situ).

    pixel_Nucleus  = CC_Nucleus.PixelIdxList;  %Get pixel list of each cell (cell array format) from Nucleus-1
    pixel_InSitu = CC_inSitu.PixelIdxList; %Get pixel list of each cell (cell array format) from In Situ image

    EachCellPixel_Nucleus = cellfun('length',CC_Nucleus.PixelIdxList);      %Get array of the number of pixels for each cell in Nucleus
    EachCellPixel_InSitu = cellfun('length',CC_inSitu.PixelIdxList);    %Get array of the number of pixels for each cell in InSitu


    poolobj = gcp('nocreate');
    delete(poolobj);
    no_of_workers = 12;
    parpool ('local',no_of_workers);  

    Nucleus_Idx_cellID = {};                  %Preallocation
    parfor o=1:CC_Nucleus.NumObjects          %Make cell array with cell ID
        Nucleus_Idx_cellID{o} = [pixel_Nucleus{o} ones(length(pixel_Nucleus{o}), 1) .* o];
        disp(['Object number=' num2str(o) 'is done!']);
    end

    InSitu_Idx_cellID = {};                  %Preallocation
    parfor m=1:CC_inSitu.NumObjects          %Make cell array with cell ID
        InSitu_Idx_cellID{m} = [pixel_InSitu{m} ones(length(pixel_InSitu{m}), 1) .* m];
        disp(['Object number=' num2str(m) 'is done!']);
    end

    Nucleus_Idx_mat = cell2mat(Nucleus_Idx_cellID');    %Convert cell array to matrix with cell ID
    InSitu_Idx_mat = cell2mat(InSitu_Idx_cellID');

    disp(['Saving matrix result....']);
        save([out_dir '\' 'Nucleus_Idx_mat.mat'],'Nucleus_Idx_mat');
        save([out_dir '\' 'InSitu_Idx_mat.mat'],'InSitu_Idx_mat');

    clear Nucleus_Idx_cellID InSitu_Idx_cellID pixel_Nucleus pixel_InSitu;
    
    
%% Colocalize check in Nucleus
    
    % Final aim is to get pixel_overlap ratio (InSitu_pixel/Nucleus_pixel)

    % Make a list[nucleus_pixel_index, nucleus_cell_ID, InSitu_pixel_index,InSitu_cell_ID]
    % Variable is Nucleus_Idx_mat_merge

        [Lia,Locb] = ismember(InSitu_Idx_mat(:,1), Nucleus_Idx_mat(:,1)); % Check if InSitu_pixel is detected in Nucleus_pixel and index of them
        %Lia is the same dimension with InSitu_Idx_mat showing 0 (not detected), 1(detected)

        search_pixel = find(Lia);  % Get pixel position list of InSitu_pixel detected in Nucleus_pixel (removing non-detected pixel)

        cell_ind = zeros(length(Nucleus_Idx_mat),2);  % preallocation for next for loop

        for i=1:length(search_pixel)    % make list for Nucleus_mat which shows colocalizing InSitu pixel and cell number
            cell_ind(Locb(search_pixel(i)), 1:2) = InSitu_Idx_mat(search_pixel(i), 1:2);   
            % The InSitu pixel found detected in in Nucleus_pixel is registered to the cell_ind (Nucleus pixel_array)
        end

        Nucleus_Idx_mat_merge = cat(2, Nucleus_Idx_mat, cell_ind); % Make one array for colocalization (for data check)
        % Nucleus_Idx_mat_merge is made here
        varNames = {'Nucleus_pixel','Nucleus_cell_ID', 'Pixels_InSitu','InSitu_cell_ID'};        
        pixel_count_table = horzcat(Nucleus_Idx_mat, cell_ind);                 
        pixel_count_table = array2table(pixel_count_table, 'VariableNames',varNames);
        
    % Counting the number of pixels overlapping in each nuclear cells
        coloc_pixel_Nucleus=find(Nucleus_Idx_mat_merge(:,3));       % Get the position list of Nucleus pixel which has colocalized with InSitu
        coloc_pixel_Nucleuscell=Nucleus_Idx_mat(coloc_pixel_Nucleus,2);  % Get the list of cell number of Nucleus_pixel colocalized with InSitu
        %The code below is the same meaning
        %coloc_pixel_Nucleuscell2 = Nucleus_Idx_mat_merge(Nucleus_Idx_mat_merge(:,3) ~= 0, 2);

    % Count total pixels colocalized with InSitu in each cell;
        edges_one= cat(1, unique(coloc_pixel_Nucleuscell), max(coloc_pixel_Nucleuscell)+1);
        [N,edges] = histcounts(coloc_pixel_Nucleuscell,edges_one);
        % N shows the number of pixels colocalized in Nuclear_cells listed in edges
        % edges is the list of Nuclear_cell ID which have overlap with InSitu (N is the number of pixels overlapping).

        
    % Calculate pertange of colocalization
    
        pixel_number = EachCellPixel_Nucleus(edges(1:end-1));        
        % Get the list of the pixel_number of *only Nucleus cells colocalized*
        % Last edges does not have N value, so (1:end-1)
    
        coloc_Percent = N' ./ pixel_number' * 100;  % Actual calculation of ratio here
        % N is the number of pixels of in situ in cells, pixel_number is the number of pixels of all cell 
        
    % Make a table of the result  
        pixel_count = horzcat(edges(1:end-1), N', pixel_number', coloc_Percent);                 
        % Make array to easy to see the result of histocounts and colocalization
        
        % Make table for save
        varNames = {'Nucleus_ID','Pixels_InSitu','Pixels_Nucleus','Colocalization_ratio'};        
        pixel_count_table = table(edges(1:end-1), N', pixel_number', coloc_Percent, 'VariableNames',varNames);                 

        
    % Make a list of Nucleus cell and colocalization ratio
        Nucleus_cell_pixel_coloc = zeros(CC_Nucleus.NumObjects,2); % Preallocation of an array for the number of Nucleus
        Nucleus_cell_pixel_coloc(:,1) = EachCellPixel_Nucleus'; % First column is the number of pixels for each cells

        for i=1:length(coloc_Percent) % To remove last edge value
            Nucleus_cell_pixel_coloc(edges(i),2) = coloc_Percent(i);      
            %Generate matrix of 1)the number of pixel of each cell, 2)percent of colocalization with InSitu
        end

    % Saving and clear temporal variables
        disp(['Saving colocalize percentage result....']);
        save([out_dir '\' 'Nucleus_cell_pixel_coloc.mat'],'Nucleus_cell_pixel_coloc');
        save([out_dir '\' 'Nucleus_Idx_mat_merge.mat'],'Nucleus_Idx_mat_merge');
        save([out_dir '\' 'pixel_count.mat'],'pixel_count');
        writetable(pixel_count_table, [out_dir '\' 'pixel_count_table.csv'])
        
        clear cell_ind Lia Locb N edges coloc_pixel_Nucleus edges_one coloc_Percent search_pixel pixel_count_table ...
              Nucleus_Idx_mat varNames coloc_pixel_Nucleuscell InSitu_Idx_mat EachCellPixel_Nucleus EachCellPixel_InSitu ...
              pixel_count pixel_number;

%% Checking if InSitu staining overlap more than one Nucleus cell

    insitu_cellID = Nucleus_Idx_mat_merge(:,4) > 0;   %select colocarized pixel in Nucleus (logical =1)
    colocIdx_list = Nucleus_Idx_mat_merge(insitu_cellID, :);   %get Nucleus cell pixels overlap with insitu

    search_cell_Nucleus = unique(colocIdx_list(:,2)); % Get cell list of nucleus (unique)
    search_cell_insitu = unique(colocIdx_list(:,4));  % Get cell list of InSitu  (unique)

    for i=1:length(search_cell_insitu) % To search InSitu cells
    
        temp_cell = colocIdx_list(:,4) == search_cell_insitu((i),1);   % get logical position of each insitu cell pixel
        Nucleus_colcell = colocIdx_list(temp_cell, 2);             % get all Nucleus cell ID in one in situ cell
        check_unique = unique(Nucleus_colcell);                    % get unique cell ID number of Nucleus cell
    
        [Lia,Locb] = ismember(check_unique, search_cell_Nucleus(:,1));    % get the position of Nucleus cells in the search list
    
        if length(check_unique) == 1    %if the InSitu cell is colocalized with other Nuclear cells
        search_cell_Nucleus(Locb,2) =1; %if in situ cell colocalized with single Nucleus cells, then =1, and these cell will be kept
        else                            %if not
        search_cell_Nucleus(Locb,2) =0; %if in situ cell colocalized with many Nucleus cell, then =0, and these cell will be removed
        end 
        disp(['Object number=' num2str(i) 'is done!']);
    end

%% Checking if Nucleus staining overlap more than one insitu cell

    for i=1:length(search_cell_Nucleus)     % Search all Nuclear_cells
    
        temp_cell = colocIdx_list(:,2) == search_cell_Nucleus((i), 1);  % get logical position of Nucleus cell pixel
        insitu_colcell = colocIdx_list(temp_cell, 4);                   % get all InSitu cell ID in one in Nucleus cell
        check_unique = unique(insitu_colcell);                          % check unique number of in situ cell number
    
      [Lia,Locb] = ismember(check_unique, search_cell_insitu(:,1));   % get the position of in situ cell in the search list
    
        if length(check_unique) == 1   %if the cell is colocalized with other cells
        search_cell_insitu(Locb,2) =1; %if Nucleus cell colocalized with single in situ cells, then =1, and these cell will be kept
        else
        search_cell_insitu(Locb,2) =0; %if YORP celll colocalized with many in situ cell, then =0, and these cell will be removed
        end 
        disp(['Object number=' num2str(i) 'is done!']);
    end

    disp(['Saving colocalize cell calculation....']);
    save([out_dir '\' 'search_cell_Nucleus.mat'],'search_cell_Nucleus');
    save([out_dir '\' 'search_cell_insitu.mat'],'search_cell_insitu');
    save([out_dir '\' 'colocIdx_list.mat'],'colocIdx_list');
    
    % Here search_cell_insitu and search_cell_Nucleus are arrays of 
    % [cell_ID_list, single_colocarize_flag](1 is yes, 0 is no).
    
%% Summarize results

    % Extract the Nucleus cells which have single cell overlap with InSitu and vice versa
    colNucleus = search_cell_Nucleus(search_cell_Nucleus(:,2) == 1) ; % Extract cell IDs which is colocalized in sigle cell manner
    colinsitu = search_cell_insitu(search_cell_insitu(:,2) == 1) ;    % Extract cell IDs which is colocalized in sigle cell manner

    [Logical,Location] = ismember(colocIdx_list(:,2), colNucleus);    % find position of single cell positive Nucleus ID.
    colocIdx_check1 = colocIdx_list (Logical,:);                      % Get the Nucleus cells with single InSitu positive (but still other Nuclear_cell can overlap the same InSitu cells)
    [Logical,Location] = ismember(colocIdx_check1(:,4), colinsitu);   % find position of single cell positive both Nucleus and InSitu cells.
    colocIdx_check2 = colocIdx_check1 (Logical,:);                    % Get the Nucleus cells with single InSitu positive and no colocalize found in othe Nucleus cells

    coloc_cell_list = unique([colocIdx_check2(:,2), colocIdx_check2(:,4)], 'rows'); %colocalized cell ID list Nucleus(1) and InSitu cell(2) ID

    %
    Idx_percent = zeros(length(coloc_cell_list),2); % Preallocation
    Insitu_flag = zeros(length(Nucleus_cell_pixel_coloc),2); % Preallocation
    for i=1:length(coloc_cell_list) % For extracted cell list
        Idx_percent(i,:) = Nucleus_cell_pixel_coloc(coloc_cell_list(i), :); % Get the colocalized rate in pixel [num of pixel in each cell, coloc_rate]
        Insitu_flag(coloc_cell_list(i),1) = coloc_cell_list(i, 2); %In situ cell ID colocalized
        Insitu_flag(coloc_cell_list(i),2) = 1;  % Single cell staining confirmatin flag
    end    
    
    coloc_cell_list = cat(2, coloc_cell_list, Idx_percent);
    % the coloc_cell_list is an arrry of [Nucleus_cell_ID, InSitu_Cell_ID, num of pixel in each cell, coloc_rate]
    Nucleus_coloc_list = cat(2, Nucleus_cell_pixel_coloc, Insitu_flag);
    % the Nucleus_coloc_list is an arrry of [Nucleus_cell_pixel_num, colocalize_rate, InSitu_cell_ID, Flag_for_single_cell_overlap]
  
    disp(['Saving colocalize cell list....']);
    save([out_dir '\' 'coloc_cell_list.mat'],'coloc_cell_list');    
    %1:Nucleus cell ID, 2:InSitu cell ID, 3:Nucleus_Pixel number, 4:(InsituPixel/POPRO Pixel)*100 
    save([out_dir '\' 'Nucleus_coloc_list.mat'],'Nucleus_coloc_list');  
    %1:Nucleus pixel number, 2:colocalized percentage, 3:InSitu cell ID, 4:flag
    %for single colocalization

%% Merge with distance map

    Nucleus_coloc_distancemap_list = cat(2, Merge_distance, Nucleus_coloc_list);
    disp(['Saving distance-colocalize cell list....']);
    save([out_dir '\' 'Nucleus_coloc_distancemap_list.mat'],'Nucleus_coloc_distancemap_list');   
    %Make csv file with header
    %This code is from 'CSV with column headers'
    header_string ='ID,X,Y,Z,Distance,pixel_number,colc_percentage,InSitu_cell_ID, SingleFlag';
    filename=[out_dir '\' 'Nucleus_coloc_distancemap_list.csv'];   %'Merge_mean_table.csv'
    fid = fopen(filename,'w');
    fprintf(fid,'%s\r\n',header_string);
    fclose(fid);
    dlmwrite(filename, Nucleus_coloc_distancemap_list, '-append','delimiter',',', 'precision', 8);

    InSitu_distance = Nucleus_coloc_distancemap_list(Nucleus_coloc_distancemap_list(:,8) > 0, :);
    save([out_dir '\' 'InSitu_distance.mat'],'InSitu_distance');   
    header_string ='ID,X,Y,Z,Distance,pixel_number,colc_percentage,InSitu_cell_ID, SingleFlag';
    filename=[out_dir '\' 'InSitu_distance.csv'];   %'Merge_mean_table.csv'
    fid = fopen(filename,'w');
    fprintf(fid,'%s\r\n',header_string);
    fclose(fid);
    dlmwrite(filename, InSitu_distance, '-append','delimiter',',', 'precision', 8);
    
    
end  

