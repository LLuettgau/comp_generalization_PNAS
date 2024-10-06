%function to load questionnaire and task data from csv/xlsx files

function task_data = load_data(data_dir)

    cd(data_dir)

    task_data = readtable('all_data.csv');

end