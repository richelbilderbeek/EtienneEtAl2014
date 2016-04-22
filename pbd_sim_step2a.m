function [ts,its] = pbd_sim_step2a(ti,tc,te,pa)
    % this computation needs:
    % - vector ti of speciation-initiation times
    % - vector tc of speciation-completion times
    % - vector te of extinction times
    % - vector pa of parent species

    ig=zeros(1,length(ti));
    ig(te==-1 & tc~=-1)=1;
    ig(te==-1 & tc==-1)=-1;
    if te(1)==-1,
        ig(1)=1;
    end

    igg=ig; % copy of table of good/incipient flags
    ppa=pa; % copy of table of parent indices
    its=0; % index that will run through table
    ts=zeros(1,1e3); % table of splitting times
    pp=zeros(1,1e3); % table of parent indices
    dd=zeros(1,1e3); % table of daughter indices
    idxs=find(igg~=0);
    while idxs(end)>1,
        [~,idx]=max(ti(idxs));
        di=idxs(idx); % daughter index
        pi=ppa(di); % parent index (can be negative!)
        if igg(abs(pi))==1 && pi>0 && igg(di)==-1,
            %disp('parent alive, good at event, good at present, daughter inc at present')
            igg(di)=0;
            igg(abs(pi))=1;
        elseif igg(abs(pi))==1 && pi>0 && igg(di)==1,
            %disp('parent alive, good at event, good at present, daughter good at present')
            igg(di)=0;
            igg(abs(pi))=1;
            its=its+1;
            ts(its)=ti(di);
            pp(its)=abs(pi);
            dd(its)=di;
        elseif igg(abs(pi))==1 && pi<0 && igg(di)==-1,
            %disp('parent alive, inc at event, good at present, daughter inc at present')
            igg(di)=0;
            igg(abs(pi))=1;
        elseif igg(abs(pi))==1 && pi<0 && igg(di)==1,
            %disp('parent alive, inc at event, good at present, daughter good at present')
            igg(di)=0;
            igg(abs(pi))=1;
            its=its+1;
            ts(its)=ti(di);
            pp(its)=abs(pi);
            dd(its)=di;
        elseif igg(abs(pi))==-1 && pi<0 && igg(di)==-1,
            %disp('parent alive, inc at event, inc at present, daughter inc at present')
            igg(di)=0;
            igg(abs(pi))=-1;
        elseif igg(abs(pi))==-1 && pi<0 && igg(di)==1,
            %disp('parent alive, inc at event, inc at present, daughter good at present')
            igg(di)=0;
            igg(abs(pi))=1;
        elseif igg(abs(pi))==0 && pi>0 && igg(di)==-1,
            %disp('parent dead, good at event, daughter inc at present')
            igg(di)=0;
            igg(abs(pi))=1;
        elseif igg(abs(pi))==0 && pi>0 && igg(di)==1,
            %disp('parent dead, good at event, daughter good at present')
            igg(di)=0;
            igg(abs(pi))=1;
        elseif igg(abs(pi))==0 && pi<0 && igg(di)==-1,
            %disp('parent dead, inc at event, daughter inc at present')
            igg(di)=0;
            igg(abs(pi))=-1;
        elseif igg(abs(pi))==0 && pi<0 && igg(di)==1,
            %disp('parent dead, inc at event, daughter good at present')
            igg(di)=0;
            igg(abs(pi))=1;
        end
        idxs=find(igg~=0);
end
