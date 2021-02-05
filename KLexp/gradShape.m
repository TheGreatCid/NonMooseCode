function G = gradShape(x)
                   
           G = [-x(2)-1 x(2)+1 1-x(2) x(2)-1;...
               1-x(1) x(1)+1 -1-x(1) x(1)-1];
end