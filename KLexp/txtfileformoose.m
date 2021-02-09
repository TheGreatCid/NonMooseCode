function txtfileformoose(name,x,y,data)
    
    xmin = min(x);
    xmax = max(x);
    

    fid = fopen(name,'wt');
    
    fprintf(fid,'AXIS X\n');
    
    for i = 1:length(x)
        fprintf(fid,'%f ',x(i));
    end
    
    fprintf(fid,'\nAXIS Y\n');
    
    for i = 1:length(y)
       fprintf(fid,'%f ',y(i)); 
    end
    
    fprintf(fid,'\nDATA\n');
    
    for i = 1:length(data)
        fprintf(fid,'%f ',data(i));
    end
    fclose(fid);


end