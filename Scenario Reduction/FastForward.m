function [ J ] = FastForward ( scenarios, probs, iterations )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    n=size(scenarios,1);
    cku1=zeros(n,n);
    ckui_1=zeros(n,n);
    ckui=zeros(n,n);
    zu1=zeros(n,1);
    J1=ones(n,1);
    % Step 1
    for k=1:n
        for u=1:n
            cku1(k,u)=cku(scenarios(k,:),scenarios(u,:));
        end
    end
    for u=1:n
        zu1(u)=sum(probs.*cku1(:,u))-cku1(u,u);
    end
    [dumm,u1]=min(zu1);
    J1(u1)=0;
    
    ckui_1=cku1;
    Ji_1=J1;
    ui_1=u1;
    
for i=2:iterations
    %Step i
    for k=1:n
        if (Ji_1(k)==1)
            for u=1:n
                if (Ji_1(u)==1)
                    ckui(k,u)=min(ckui_1(k,u),ckui_1(k,ui_1));
                end
            end
        end
    end
    
    zui=zeros(n,1);
    for u=1:n
        if (Ji_1(u)==1)
            for k=1:n
                if (Ji_1(k)==1 && k~=u)
                    zui(u)=zui(u)+probs(k)*ckui(k,u);
                end
            end
        else
            zui(u)=inf;
        end
    end
    
    [dumm,ui]=min(zui);
    Ji=Ji_1;
    Ji(ui)=0;
    
    ckui_1=ckui;
    Ji_1=Ji;
    ui_1=ui;
end
J=Ji;
end

