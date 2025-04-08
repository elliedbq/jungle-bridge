%This block of code demonstrates how to data from
%an excel file into MATLAB
function data = load_excel_example()

    fpath = '.\';
    fname = 'RubberBandTemplate.xlsx';
   
    my_table = readtable([fpath,fname]);
   
    disp(my_table);

    row_range = 1:12;
    col_range = 3:6;
    disp(my_table(row_range,col_range));
    
    data_mat = table2array(my_table(row_range,col_range));

    disp(data_mat);

    data = data_mat;
end