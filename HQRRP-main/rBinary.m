function M = rBinary(file,row,column,type)
    % INPUT VARIABLES
    % file - file for which to read from.
    % row - number of rows to read.
    % column - number of columns to read.
    % type - 'd' or 'D' for double, and 'i' or 'I' for integer.
    % OUTPUT VARIABLES
    % M - place to store the data from the file.
    
    % open file for reading
    binFile = fopen(file,'r');
    
    % read data from binary file
    if upper(type) == 'D'
        % read data of type double from the file
        M = fread(binFile,[row column],'double');
    end
    if upper(type) == 'S'
        % read data of type double from the file
        M = fread(binFile,[row column],'single');
    end
    if upper(type) == 'I'
        % read data of type integer from the file
        M = fread(binFile,[row column],'int');
    end
    
    % close the file
    fclose(binFile);
end
