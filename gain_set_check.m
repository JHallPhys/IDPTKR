function Hus=gain_set_check(k,N,gamma,nfq)

fcheck=fname_gset(k,N,gamma,'.mat') % The filename the matrix should have;

% return
% Now change to the new diretory and check if file exists
Hus=zeros(N,N);
parent_d = cd;
cd './Husimi_dat' % Directory where matrix is stored
if isfile(fcheck) % File exists
     'file exists'
    Hus = matfile(fcheck);
    Hus=Hus.Hus; % I think this step may be redundent 
    cd(parent_d)
else % File does not exist
    'file does not exist'
    cd(parent_d) % Go back to current directory to construct matrix
    Hus=zeros(N,N);
    tic
    for itt_efn = 1:round(nfq*N)
    itt_efn
    cd './Husimi_dat' % Directory where matrix is stored
    fname_efn=strcat('Husimi_Entropy_k',num2str(k),'_g0p001_N',num2str(N),'_single_efn',num2str(itt_efn),'_special');
    Hus_Entropy = matfile(fname_efn);
    Hus_Entropy=Hus_Entropy.Hus_Entropy; % I think this step may be redundent
    cd(parent_d)
    Hus=Hus+Hus_Entropy;
    
    end
    cd './Husimi_dat' % Directory where matrix is stored
    save(fcheck,'Hus'); % save it 
    'file created'
    cd(parent_d)
end




end