% h5writeMatrix(f_h5,string_data,data)
function h5writeMatrix(f_h5,string_data,data)
try; 
    h5create(f_h5,string_data,size(data));
end
h5write(f_h5,string_data,data);