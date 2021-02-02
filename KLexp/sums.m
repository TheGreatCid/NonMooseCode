%David Torres
%For use in translating RF isoparametric coords to FE mesh
function F = sums(x1,x2,m,FEel,el,nodes,FEnodes)
    for i = 1:4
        Val1 = [(1/4)*(1-x1)*(1-x2'),(1/4)*(1+x1)*(1-x2'),(1/4)*(1+x1)*(1+x2'),(1/4)*(1-x1)*(1+x2')]*el(:,m) + nodes(el(:,1),1)' + FEnodes(FEel(:,1),1)';
        Val2 = [(1/4)*(1-x1)*(1-x2'),(1/4)*(1+x1)*(1-x2'),(1/4)*(1+x1)*(1+x2'),(1/4)*(1-x1)*(1+x2')]*el(:,m) + nodes(el(:,1),2)' + FEnodes(FEel(:,1),2)';
    end
    F(:,1) = Val1;
    F(:,2) = Val2;
end