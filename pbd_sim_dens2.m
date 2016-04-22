function [results] = pbd_sim_dens2(runoncluster,mui,la2i,agei,seed,startmc,endmc)
  % Simulating number good and incipient species
  if ischar(runoncluster) runoncluster = str2num(runoncluster); end;
  if ischar(mui) mui = str2num(mui); end;
  if ischar(la2i) la2i = str2num(la2i); end;
  if ischar(agei) agei = str2num(agei); end;
  if ischar(startmc) startmc = str2num(startmc); end;
  if ischar(endmc) endmc = str2num(endmc); end;
  if ischar(seed) seed = str2num(seed); end;
  if runoncluster == 1
      outputdir = '';
  else
      outputdir = 'd:\data\ms\posterior\';
  end;
  rand('twister',seed);
  results = zeros(endmc,5);
  la1 = 0.5;
  muvec = [0,0.1,0.2];
  la2vec = [0.1,0.3,1,0.03];
  agevec = [5,10,15];
  stem_or_crown = 1;
  eps = 1E-6;
  outfilename2 = sprintf('%spbd_sim_dens2_%d-%d-%d.out',outputdir,mui,la2i,agei);
  for mc = 1:endmc  
      disp(mc);
      la2 = la2vec(la2i);
      mu = muvec(mui);
      age = agevec(agei);
      pars = [la1,mu,la2,age];
      tssim = durspec(log(pars));
      times = [];
      for ii = 0:stem_or_crown
          valid = 0;
          while not(valid)
             [ti,tc,ig,pa,valid] = pbd_sim_recstr_full(pars);
          end;        
          te = zeros(1,length(ti));
          te(1:sum(ig~=0)) = -1;
          [times1,~] = pbd_sim_step2a(ti,tc,te,pa);        
          times = [times,age-times1(times1>0)];        
      end;
      times = [age,times];
      times = sort(times(times > 0),'descend');
      outfile2 = fopen(outfilename2,'at');
      fprintf(outfile2,'%0.6f ',times);
      fprintf(outfile2,'\n');
      fclose(outfile2);
  end;
end
