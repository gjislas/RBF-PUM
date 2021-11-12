function wBinary(M,file,prec) 
    % INPUT VARIABLES
    % M - contains data to be written
    % file - file for which to write to
    
    % open file for writing
    binFile = fopen(file,'w');
    % write data to file in double precision
    fwrite(binFile,M,prec);
    % close the file
    fclose(binFile);
end