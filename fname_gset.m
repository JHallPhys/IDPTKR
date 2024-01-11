function fname=fname_gset(k,N,gamma,str_ext)

ds='_'; % seperator of strings
str_name='Husimi_gain_set';
str_kick='k';
str_kick_num= strrep(num2str(k),'.','p');
str_gamma='g';
str_gamma_num= strrep(num2str(gamma),'.','p');
str_N='N';
str_N_num=num2str(N);

fname=strcat(str_name,ds,str_kick,str_kick_num,ds,str_gamma,str_gamma_num,ds,str_N,str_N_num); % The information about the N range
fname=strcat(fname,str_ext); % Add the extension

end