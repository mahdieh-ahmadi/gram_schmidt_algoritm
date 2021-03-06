function [independent , dependents , solves] = gram_schmidt_function(inputMatrix)
sizeMatrix = size(inputMatrix);
counter = sizeMatrix(2);
q=inputMatrix(:,1)/norm(inputMatrix(:,1));
dependents = [];
zero = zeros(sizeMatrix(1),1);
 

for i=2:1:counter
     ak = inputMatrix(:,i);
     qk = ak;
    for j=1:1:i-1
       qk = qk - ( q(:,j)'*ak).*q(:,j);
    end
    if (norm(qk)<0.01)
        dependents =[dependents , i];
        qk = zero;
        q =[q, qk];
    else
    q =[q, qk/norm(qk)];
    end
end


independent = [];
for i=1:1:counter
   if (sum(q(:,i)-zero) ~= 0)
       independent = [independent , q(:,i)];
   end
end

solves = [];
if(size(dependents) ~= 0 )
    coefficient = sym('coefficient' , [sizeMatrix(1),size(dependents)]);
    side = size(dependents);
    independentsize = size(independent);
    solves = zeros(sizeMatrix(1),side(2));
    for i=1:1:side(2)
        modle = 0;
        for j=1:1:size(independent)
            modle = modle + coefficient(j,i).*independent(:,j);
        end
        answers = solve(modle == inputMatrix(:,dependents(i)),coefficient(:,i).');
        c = struct2cell(answers);
        solves(:,i) = [c{:}];
    end
end



end

