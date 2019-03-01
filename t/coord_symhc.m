function [I]=coord_symhc(N,maxd)
    a1 = (0:N)';
    weight = max(1, double(a1));
    for d=2:maxd
        a2 = []; 
        weight2 = [];
        for k = 1:size(a1,1)
            tmp = (0:N/weight(k)); 
            a2 = [a2; [(ones(1,length(tmp))'*double(a1(k,:))),tmp']];
            tmp = max(1, double(tmp));
            weight2 = [weight2; weight(k)*tmp'];
        end %for
        a1 = a2;
        weight = weight2;
    end %for
    clear a2 weight weight2 ind ind0 tmp;
    if(prod(size(a1))>0)
        I = sym_set(a1);
    else
        I = zeros(0,maxd);
    end %if
end

function b=sym_set(a)
    tmp=[a];
    tmp=[tmp;-tmp(:,1),tmp(:,2:size(a,2))];
    tmp=unique(tmp, 'rows');
    for j=2:size(a,2)-1
        tmp=[tmp;tmp(:,1:j-1),-tmp(:,j),tmp(:,j+1:size(a,2))];
        tmp=unique(tmp, 'rows');
    end
    tmp=[tmp;tmp(:,1:size(a,2)-1),-tmp(:,size(a,2))];
    tmp=unique(tmp, 'rows');
    b=tmp;
end