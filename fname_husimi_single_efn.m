function fname=fname_husimi_single_efn(K_class,N,gamma,nset,str_ext)

gamma=imag(gamma);
ds='_'; % seperator of strings
str_name='Husimi_Entropy';
str_kick='k';
str_kick_num= strrep(num2str(K_class),'.','p');
str_gamma='g';
str_gamma_num= strrep(num2str(gamma),'.','p');
str_N='N';
str_N_num=num2str(N);
str_subset='single_efn';
str_subset_num=num2str(nset);

fname=strcat(str_name,ds,str_kick,str_kick_num,ds,str_gamma,str_gamma_num,ds,str_N,str_N_num,ds,str_subset,str_subset_num); % The information about the N range
fname=strcat(fname,str_ext); % Add the extension

end