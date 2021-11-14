function mainProgramExecutor()
% Predictive distribution for Paper Helicopter with MCMC simulation
% 
% See also modelSelector, VV_executor, cdCalculator, fallTimeCalculator
 
%%% Importing the test data 
im_data = importdata('final_test_data.xlsx'); % Importing data
% im_data = importdata('Final paperheli_timeonlydata.xlsx'); % Importing data
sheet_names = fieldnames(im_data.data); % Get the sheet names of the xlsx file
run_sheet_number = 1; % Select the sheet number you want to run with
 
for kkk = 1:1:length(run_sheet_number)
% Get the time data and text data from the selected sheet
eval(['time_data = im_data.data.',sheet_names{run_sheet_number(kkk)},';']);
eval(['text_data = im_data.textdata.',sheet_names{run_sheet_number(kkk)},';']);
 
% Get the model type from the selected sheet
model_text = text_data{3};
model = model_text(strfind(model_text,':')+2);
cd_model_kind = {'linear','quadratic'};
simul_kind = 'pd' % 'pbox'
 
    for kkkk = 1:1:length(cd_model_kind)
        % Select the dynamics model (linear / quadratic)
        cd_model = cd_model_kind{kkkk}; 
        
        %%% Select the column index range of the data
        %%% 2.5m 1~6 / 4.5m 7~12 / 8.17m 13~end for S/M/L/L-add sheet
        %%% 2.4m 1~6 / 4.45m 7~12 for L-add2 sheet
        %%% 6.82m 1~6 / 10.67m 7~12 for S-add sheet
        column_range = 1:12; % 1st~12th index
        
        final_data = time_data(:,column_range); % Get the fall time data
        text_data = text_data(5,column_range); % Get the header data
        
        %%% Selecting the Calibration conditions
        % Select the column index for calibration conditions
        calib_index = [7 8 9]; 
        % Get the calibration conditions (height, no. of clips)
        [~,calib_clip_no,calib_h,~] = parse_text(text_data{1,calib_index(1)});    
         % Initialize the variable as zeros to put the data in      
        total_area_metric = zeros(1,length(column_range)); 
        for i=1:1:length(column_range)
            
            % Get the current fall time data
            time_data = final_data(:,i);
            % Get the current header data
            current_model = text_data{1,i}; 
            
            % Get the helicopter conditions from the header data
            [heli_no,clip_no,fall_height,~] = parse_text(current_model);
                       
            %%% Generate Calibration Cd data depeding on helicopter number
            switch heli_no
%                 case 0
%                     calib_data = final_data(:,calib_index(1));                               
%                     [rho, g, w, a] = modelSelector(model, calib_clip_no);
%                     calib_data_cd = cdCalculator(cd_model,w,rho,a,calib_h,g,calib_data);                    
                case {1,4}                   
                    calib_data = final_data(:,calib_index(1));            
                    [rho, g, w, a] = modelSelector(model, calib_clip_no);
                    calib_data_cd = cdCalculator(cd_model,w,rho,a,calib_h,g,calib_data);
                case {2,5}
                    calib_data = final_data(:,calib_index(2));               
                    [rho, g, w, a] = modelSelector(model, calib_clip_no);
                    calib_data_cd = cdCalculator(cd_model,w,rho,a,calib_h,g,calib_data);
                case {3,6}
                    calib_data = final_data(:,calib_index(3));                  
                    [rho, g, w, a] = modelSelector(model, calib_clip_no);
                    calib_data_cd = cdCalculator(cd_model,w,rho,a,calib_h,g,calib_data);
            end

     %%% Executing MCMC simulation
            [area_metric, tha] = VV_executor(model,fall_height,clip_no,...
                time_data,calib_data_cd,cd_model,simul_kind);
            total_area_metric(i) = area_metric;
 
            gcf; % Activate the current figure in MATLAB
            
            %%% Putting the title on the current graph depending on the conditions 
            switch i
                case num2cell(calib_index)
                % For the calibration conditions
                    title([current_model,' for calibration ',cd_model,' model']);
                otherwise
                % For the others
                    title([current_model,' ',cd_model,' model']);
            end
            %%% Saving the current graph to the png file        
            % Making the folder
            if isdir([model,' Model ',cd_model,'Cd table3 data ', simul_kind]) == 0
            % Check if there is the folder of which name is the same in the Current Folder
                mkdir([model,' Model ',cd_model,'Cd table3 data ', simul_kind]);
                % Generate the folder with this name 
            else
                % Do nothing
            end
            oldcd = cd([model,' Model ',cd_model,'Cd table3 data ', simul_kind]);
            % Changing the Current Folder to the folder 
            
            % Making PNG files
            temp_figure = findall(0,'type','figure');
            % Get the handler of all the figure
            filename = [model,num2str(heli_no),'_',num2str(fall_height),'_',num2str(clip_no)];
            % Select the filename for '.fig' and '.png' file
            saveas(gcf, filename, 'fig'); 
            % Saving the current figure as '.fig' file
            for ff = 1:1:length(temp_figure)
               print(temp_figure(ff), [filename,'.png'],'-dpng','-opengl');
               % Saving the current figure to the '.png' file
            end
            cd(oldcd);
            % Chainging the Current Folder to the previous folder
            close all 
            % Closing the all the figure showing on the MATLAB
            
            %%% Calculate mean of the mean/mean of the std/std of the mean/std of the std
            [mean(tha(:,1));std(tha(:,1));mean(tha(:,2));std(tha(:,2));]
        end
        % total_area_metric = reshape(total_area_metric,3,4);
        % total_area_metric = total_area_metric';
        save([model,' Model ',cd_model,'Cd predictive ', simul_kind],'total_area_metric','text_data');
        % Saving the 'total_area_metric' variable to the '.mat' file
    end
    
end

function [heli_no,clip_no,fall_height,height_char]=parse_text(current_model)
%%% Parse the column header into heli_no, clip_no, and fall height
%%% 
 
slash_index = strfind(current_model,'/');
if length(slash_index) == 2
%%% Current format of header S1/2.5 ~ S3/2.5    
    heli_no = str2double(current_model(2)); % no. of Helicopter
    clip_no = str2double(current_model(slash_index(2)+1:end));
    height_char = current_model(slash_index(1)+1:slash_index(2)-1);
    fall_height = str2double(height_char)*100; % Convert m to cm              
else
    %%% The previos format of header S1/2.5 ~ S6/2.5    
    heli_no = str2double(current_model(2)); % no. of Helicopter
    switch heli_no
        case {1,2,3}
            % 1~3 helicopter must have 1 clips
            clip_no = 1;
        case {4,5,6} 
            % 4~6 helicopter must have 2 clips
            clip_no = 2;
    end
    height_char = current_model(strfind(current_model,'/')+1:end);
    fall_height = str2double(height_char)*100; % Convert m to cm               
end