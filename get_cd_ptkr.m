function CD=get_cd_ptkr(Norm_hm_av,sigma,nfq,N)

Norm_sort=Norm_hm_av(:);
Norm_sort=sort(Norm_sort,'descend');

CD=zeros(N,N);


Norm_single_state=Norm_hm_av;
Norm_single_state(Norm_single_state<=Norm_sort(round(nfq*N)*N))=0;
Norm_single_state(Norm_single_state>1)=1;
CD=imgaussfilt(Norm_single_state,sigma);


end