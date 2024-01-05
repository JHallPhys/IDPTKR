% Check if the matrix U is saved
% If it is then load it
% Otherwise calculate it and save it for later

function [psi,En]=ECheck(U,N,N_1,K_class,T,gamma,str_ext);


fcheck=fname_eigenvalues(K_class,N,imag(gamma),str_ext); % The filename the matrix should have;

% Now change to the new diretory and check if file exists

parent_d = cd;    
cd './MDat' % Directory where matrix is stored
if isfile(fcheck) % File exists
    En = matfile(fcheck);
    En=En.En; % I think this step may be redundent
    psi = matfile(fcheck);
    psi=psi.psi; % I think this step may be redundent
    'Eigenvalue file exists';
else % File does not exist
    cd(parent_d) % Go back to current directory to construct matrix
    [psi,En]= schur(U); % psi are the Schur eigenfns and En matrix of eigs
%     Edat=[psi,En];
    cd './MDat' % Directory where matrix is stored
    save(fcheck,'psi','En'); % save it 
    'Eigenvalue file created';
end
cd(parent_d)
end