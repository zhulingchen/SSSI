function ind_xc=find_zero_crossings(s,flag)
% FIND_ZERO_CROSSINGS ... returns indicies of the zero crossings of a trace
%
% ind=find_zero_crossings(s,flag)
%
% s ... the trace
% flag ... 1 means find pos-to-neg zero crossings
%         -1 means find neg-to-pos zero crossings
%          0 means find both
% ind_xc ... vector of indices of the desired zero crossings
%
% Note: This function returns the index of the nearest sample that
% satisfies the condition. It does not interpolate. So the actual zero
% crossing will likely be just somewhere near the indicated sample. To
% increase accuracy, interpolate your input trace to ten or 100 times its
% original sample rate.
%
% G.F. Margrave CREWES 2010

ind=find(s>0);%find all positive samples
ind=ind(:);
if(max(ind)<length(s))
    ind=[ind;length(s)];%seems necessary to get the last zero crossing
end
ind2=find(diff(ind)>1);%finds discontinuities in the clustering of positive samples
indpn=ind(ind2);%positive to negative zc
indnp=ind(ind2+1);%negative to positive zc
switch flag
    case {1} 
        ind_xc=zeros(size(indpn));
        for k=1:length(indpn)
            if(abs(s(indpn(k)))<abs(s(indpn(k)+1)))
                ind_xc(k)=indpn(k);
            else
                ind_xc(k)=indpn(k)+1;
            end
        end
    case {-1}
        ind_xc=zeros(size(indnp));
        for k=1:length(indnp)
            if(abs(s(indnp(k)))<abs(s(indnp(k)-1)))
                ind_xc(k)=indnp(k);
            else
                ind_xc(k)=indnp(k)-1;
            end
        end
    case {0}
        ind_xc1=zeros(size(indpn));
        for k=1:length(indpn)
            if(abs(s(indpn(k)))<abs(s(indpn(k)+1)))
                ind_xc1(k)=indpn(k);
            else
                ind_xc1(k)=indpn(k)+1;
            end
        end
        ind_xc2=zeros(size(indnp));
        for k=1:length(indnp)
            if(abs(s(indnp(k)))<abs(s(indnp(k)-1)))
                ind_xc2(k)=indnp(k);
            else
                ind_xc2(k)=indnp(k)-1;
            end
        end
        ind_xc=sort([ind_xc1;ind_xc2]);
end
if(~isempty(ind_xc))        
    if(ind_xc(end)==length(s))
        ind_xc(end)=[];
    end
end
        