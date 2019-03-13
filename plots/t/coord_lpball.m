function I = coord_lpball(N,maxd,p)
%
% I=coord_lpball(N,d,p)
%
%   computes lp-ball index sets containing all h\in\mathbb{Z}^d
%   with
%   \prod_{s=1}^d\max(1,|h_s|/gamma_s)\le N
%
%   N:      refinement
%   d:      dimension 
%   p:      parameter

    a1=(0:N)';
    weight=(double(a1)).^p;
    for d=2:maxd
        a2=[];
        weight2=[];
        for l=1:size(a1,1)
            if(weight(l)~=0)
                len = (N^p-weight(l)+eps)^(1/p);
            else
                len = N;
            end %if
            tmp=0:len;
            if(length(tmp)>0)
                a2 = [a2;[(ones(1,length(tmp))'*double(a1(l,:))),tmp']];
                tmp = (double(tmp)).^p;
                weight2 = [weight2;weight(l)+tmp'];
            end %if
        end %for
        a1=a2;
        weight=weight2;
    end %for
    clear a2 weight weight2 len tmp;
    I=sym_set(a1);
end %function

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