function G = gradShape(x)
                   
           G = [-x(2)-1 x(2)+1 1-x(2)  x(2)-1;...
                      1-x(1)  x(1)+1 -1-x(1) x(1)-1]; %Should be multiplied by (1/4) I believe
end