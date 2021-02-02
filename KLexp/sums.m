%David Torres
%For use in translating RF isoparametric coords to FE mesh
function F = sums(x1,x2,j,FEel,el,nodes,FEnodes)
    Val1 = 0;
    Val2 = 0;
    for i = 1:4
        Val1 = Val1 + [(1/4)*(1-x1)*(1-x2'),(1/4)*(1+x1)*(1-x2'),(1/4)*(1+x1)*(1+x2'),(1/4)*(1-x1)*(1+x2')]*nodes(el(:,j),1) - nodes(el(:,1),1)'%need to find out how to find the location of element within RF element;
        Val2 = Val2 + [(1/4)*(1-x1)*(1-x2'),(1/4)*(1+x1)*(1-x2'),(1/4)*(1+x1)*(1+x2'),(1/4)*(1-x1)*(1+x2')]*nodes(el(:,j),2) - nodes(el(:,1),2)';
    end
    F(:,1) = Val1;
    F(:,2) = Val2;
end