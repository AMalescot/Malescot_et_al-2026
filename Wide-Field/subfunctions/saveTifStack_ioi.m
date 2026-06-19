function saveTifStack_ioi(file,filename)


filename = strcat(filename, ".tif");

if isfile(filename) > 0
    delete(filename);
end

if contains(class(file),'double') == 1
    range = [min(file,[], "all")  max(file,[],"all")];    
    file = mat2gray(file,[-max(abs(range)) max(abs(range))]);
    file = im2uint16(file);
end

for i = 1:size(file,3)
    imwrite(file(:,:,i), filename ,'WriteMode', 'append');
end