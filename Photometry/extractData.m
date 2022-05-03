function [] = extractData(GCAMP_Save_Dir)
% Select Location of Raw Excel File Directories
PathName_Folder = uigetdir();
cd(PathName_Folder);
indiv_files = dir;
indiv_files = indiv_files(3:end);
day_directory = indiv_files(1).folder;
day_names = {indiv_files.name};
% Loop through days and sessions
for day = 1:size(day_names,2)
    %% Raw Excel File Directories
    day_folder = [day_directory '\' day_names{day}];
    cd(day_folder)
    indiv_session_names_analog = dir('*analog.csv');
    indiv_session_names_photometry = dir('*photometry.csv');
    total_mice = length(indiv_session_names_analog);
    for mouse = 1:total_mice
        %% Initiate GCAMP Data Structure
        session_file_name_analog = indiv_session_names_analog(mouse).name;
        session_file_name_photometry = indiv_session_names_photometry(mouse).name;
        GCAMP.mouseID = session_file_name_analog(1:4);
        GCAMP.training_day = indiv_files(day).name;
        GCAMP.Criteria = indiv_files(day).name(1:4); % in ms
        GCAMP.onlyID = GCAMP.mouseID;
        %% Extract photometry and behavior data from .csv files
        [GCAMP] = Photometry_Signal(GCAMP, session_file_name_photometry);
        GCAMP.beh_data = Behavior_extract(session_file_name_analog);
        %% Save GCAMP Data Structure
        cd(GCAMP_Save_Dir)
        save(['GCAMP_' GCAMP.mouseID '_' GCAMP.training_day], 'GCAMP');
        cd(day_folder)
    end
end
end

