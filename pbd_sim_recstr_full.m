function [ti,tc,ig,pa,valid] = pbd_sim_recstr_full(pars)
la1=pars(1); % speciation initiation
la2=pars(3); % speciation completion
mu=pars(2); % extinction rate
T=0; Tend=pars(4); % simulation time

Ng=1; Ni=0; % initial condition
isp=1; % index that will run through table
ig=zeros(1,1e3); ig(1)=1; % table of good/incipient flags
pa=zeros(1,1e3); pa(1)=NaN; % table of parent indices
ti=zeros(1,1e3); ti(1)=-1; % table of speciation initiation times
tc=zeros(1,1e3); tc(1)=-1; % table of speciation completion times
denomvec = [la1*[Ng Ni] la2*Ni];
denom = sum(denomvec);
T = T - log(rand(1))/denom;
while T <= Tend,
    % ONE STEP OF PROCESS
    event = sample(denomvec,1);
    switch event
        case 1, % speciation-initiation in good species
            p0t=(la1-mu)/(la1-mu*exp(-(la1-mu)*(Tend-T)));
            if rand(1)<p0t,
                isp=isp+1;
                ig(isp)=-1;
                idx=find(ig==1);
                idx=idx(randi(Ng));
                pa(isp)=idx;
                ti(isp)=T;
                tc(isp)=-1;
                Ni=Ni+1;
            end                    
        case 2, % speciation-initiation in incipient species
            p0t=(la1-mu)/(la1-mu*exp(-(la1-mu)*(Tend-T)));
            if rand(1)<p0t,
                isp=isp+1;
                ig(isp)=-1;
                idx=find(ig==-1);
                idx=idx(randi(Ni));
                pa(isp)=-idx;
                ti(isp)=T;
                tc(isp)=-1;
                Ni=Ni+1;
            end                    
        case 3, % speciation-completion
            idx=find(ig==-1);
            idx=idx(randi(Ni));
            ig(idx)=1;
            tc(idx)=T;
            Ng=Ng+1;
            Ni=Ni-1;
    end
    if Ng+Ni==0,
        break,
    end
    denomvec = [la1*[Ng Ni] la2*Ni];
    denom = sum(denomvec);
    T = T - log(rand(1))/denom;
end
T=Tend;

valid=(Ng+Ni>0);