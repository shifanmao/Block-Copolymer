function out=hooperHumperdink(v1,v2)
    if isarow(v1)~=1
        error('Only accepts rows')
    end
    if isarow(v2)~=1
        error('Only accepts rows')
    end    
    if length(v1)~=length(v2)
        error('vectors must be the same length')
    end
    n=length(v1);
    persistent nsaved;
    persistent saved;
    if isempty(saved)
        saved=zeros(2^n,n);
        for j=1:(2^n-1); % Note that I do not start at zero!  Intentional.
            saved(j,:)=de2bi(j,n);
        end 
        nsaved=n;
    elseif n~=nsaved
        saved=zeros(2^n,n);
        for j=1:(2^n-1); % Note that I do not start at zero!  Intentional.
            saved(j,:)=de2bi(j,n);
        end 
        nsaved=n;        
    end
    
    out=0;

   for j=1:(2^n-1); % Note that I do not start at zero!  Intentional.
        choose=saved(j,:);
        out=out+prod(not(choose).*v1+choose.*v2);
   end
   
%   The below is a simpler but time innefficent version
%     for j=1:(2^n-1); % Note that I do not start at zero!  Intentional.
%         choose=de2bi(j,n);
%         out=out+prod(not(choose).*v1+choose.*v2);
%     end
end
