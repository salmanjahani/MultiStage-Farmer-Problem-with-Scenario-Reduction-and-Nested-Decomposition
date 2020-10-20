function [ q ] = CalculateNewProbabilities( scenarios, probs, J )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    n=size(scenarios,1);
    q=zeros(n,1);
    dist=zeros(n,1);
        
    for i=1:n
        if (J(i)==1)
            for j=1:n
                if (J(j)==0)
                    dist(j)=cku(scenarios(i,:),scenarios(j,:));
                else
                    dist(j)=inf;
                end
            end
            [dumm, j_i]=min(dist);
            q(j_i)=q(j_i)+probs(i);
        else
            q(i)=q(i)+probs(i); 
        end
    end
    
end

