function [ o ] = cku( wk,wu )
    
    o=max([1,norm(wk-[2.5,3,20]),norm(wu-[2.5,3,20])])*norm(wk-wu);
    
end