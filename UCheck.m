% Check if the matrix U is saved
% If it is then load it
% Otherwise calculate it and save it for later

function U=Ucheck(N,N_1,K_class,T,gamma,str_ext);


fcheck=fname_umatrix(K_class,N,imag(gamma),str_ext); % The filename the matrix should have;

% Now change to the new diretory and check if file exists

parent_d = cd;    
cd './MDat' % Directory where matrix is stored
if isfile(fcheck) % File exists
    U = matfile(fcheck);
    U=U.U; % I think this step may be redundent 
    'Matrix file exists';
else % File does not exist
    cd(parent_d) % Go back to current directory to construct matrix
    U=zeros(N,N); % Initialise Flouqet matrix
    [U,time1]=UMatrix(U,N,N_1,K_class,T,gamma); % Construct Flouqet matrix
    cd './MDat' % Directory where matrix is stored
    save(fcheck,'U'); % save it 
    'Matrix file created';
end
cd(parent_d)
end