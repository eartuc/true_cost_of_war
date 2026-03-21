
clearvars -except Table1a Table1b Table2 Table3 TableA1 TableA2 TableA3 TableA4 TableB1 
tic
load data/paramters_iv.mat
load data/data_processed.mat

rng(123)

dur_all=[1:50 1000];

DN=25;

mig_r = IDP./MIG0;


%Assume that Luh and Don receive proportional IDP just to fill in zeros,
%this is not crucial.

mig_r(24,2)=mean(mig_r(1:23,2));
mig_r(25,1)=mean(mig_r(1:23,1));


s(1).index=1:23;
s(2).index=[1:4 7 9:13 15:20 22]; %robust1
s(3).index=[3 4 7]; %robust2
s(4).index=[ 9 11:13 15:18 20 22]; %robust3
s(5).index=[1:10]; %robust 4
s(6).index=[14:23]; %robust 5


MaxSim=5000;

hatU_an_simcase=zeros(3,2,51,2,MaxSim,6);
hatV_an_simcase=zeros(3,2,51,2,MaxSim,6);

new_welfare2_sig0=[];
new_welfare2_sig1=[];

for cas=1:6

%cas

clear hatU_full hatU_full_an w0_full u0_full dwU_full_an



for sim=1:MaxSim

    
NN=length(s(cas).index);
index_in = floor(rand(1,NN)*NN)+1;

if sim==1
    index_in=1:NN;
end

index_=s(cas).index(index_in);

for di=[1 10 51] %1:length(dur_all)
for bi=[1 2 3] %1:length(bta_all)
for si=[1 2]% 1:length(sigma_all)

        bta=bta_all(bi);
        sigma=sigma_all(si);
        dur=dur_all(di);

        theta=result_theta(bi,si,3);

        NN=length(index_);

        weight0=Ltm1(index_,5)/sum(Ltm1(index_,5));
        weight1=MIG0(index_,1)/sum(MIG0(index_,1));
        weight2=MIG0(index_,2)/sum(MIG0(index_,2));

        
        hatV_1 = log(repmat( [mig_r(24,1) mig_r(25,2)], [NN 1])) - log(mig_r(index_,1:2));

        hatU_2 = 0*log(repmat( [mig_r(24,1) mig_r(25,2)], [NN 1])) - log(mig_r(index_,1:2));
        
        hatV(bi,si,di,:)= sum(hatV_1(:,:).*[weight1 weight2]).*(1./theta);

        hatU(bi,si,di,:)= sum(hatU_2(:,:).*[weight1 weight2]).*(1./theta);

 
        hatU_full(bi,si,di,:,:)= (hatU_2(:,:)).*(1./theta);

        mul=(1-bta^dur)/(1-bta);


        hatV_an(bi,si,di,:)=hatV(bi,si,di,:)/mul;      
        hatU_an(bi,si,di,:)=hatU(bi,si,di,:)/mul;

        hatU_full_an(bi,si,di,:,:)=hatU_full(bi,si,di,:,:)/mul;
        
        w0(bi,si,di,:)=mean(rwaget(24:25,:),2); 
        u0(bi,si,di,:)=UU(mean(rwaget(24:25,:),2),sigma); 

        w0_full(bi,si,di,1,:)=mean(rwaget(24:25,:),2); 
        u0_full(bi,si,di,1,:)=UU(mean(rwaget(24:25,:),2),sigma); 

        dwV_an(bi,si,di,:)=100*(invUU(hatV_an(bi,si,di,:)+u0(bi,si,di,:),sigma)./w0(bi,si,di,:)-1);
        dwU_an(bi,si,di,:)=100*(invUU(hatU_an(bi,si,di,:)+u0(bi,si,di,:),sigma)./w0(bi,si,di,:)-1);
        
        dwU_full_an(bi,si,di,:,:)=100*(invUU(hatU_full_an(bi,si,di,:,:)+u0_full(bi,si,di,:,:),sigma)./w0_full(bi,si,di,:,:)-1);


end
end
end

hatU_an_simcase(:,:,:,:,sim,cas)=hatU_an(:,:,:,:);
hatV_an_simcase(:,:,:,:,sim,cas)=hatV_an(:,:,:,:);

dwV_an(3,:,2:end,:)=0;
dwU_an(3,:,2:end,:)=0;

bb=[1 2 3];
result_VU = [squeeze(hatV_an(bb,1,1,:)) squeeze(hatU_an(bb,1,1,:)); squeeze(hatV_an(bb,2,1,:)) squeeze(hatU_an(bb,2,1,:)) ];


result_UT=[];
result_VT=[];

for ii=[1 10 51]
    tmp1=[squeeze(dwU_an(bb,1,ii,:)); squeeze(dwU_an(bb,2,ii,:))]/100;
    result_UT=[result_UT tmp1];
end

for ii=[1 10 51]
    tmp1=[squeeze(dwV_an(bb,1,ii,:)); squeeze(dwV_an(bb,2,ii,:))]/100;
    result_VT=[result_VT tmp1];
end

result_VUsim(:,:,sim)=result_VU(:,:);
result_VTsim(:,:,sim)=result_VT(:,:);
result_UTsim(:,:,sim)=result_UT(:,:);


if sim==1 && cas==1

    welf90=squeeze(dwU_full_an(2,1,51,:,:));
    
    welf97=squeeze(dwU_full_an(1,1,51,:,:));

    migr2=mig_r(1:23,:);

    table_full_b90=[migr2(:,1) log(migr2(:,1)) welf90(:,1) migr2(:,2) log(migr2(:,2)) welf90(:,2)];
    
    table_full_b97=[migr2(:,1) log(migr2(:,1)) welf97(:,1) migr2(:,2) log(migr2(:,2)) welf97(:,2)];

    table_full_comb_sig0=[log(migr2(:,1)) welf97(:,1) welf90(:,1) log(migr2(:,2)) welf97(:,2) welf90(:,2) ];
    
    table_full_comb_sig0_avrs=mean(table_full_comb_sig0);
    table_full_comb_sig0_min=min(table_full_comb_sig0);
    table_full_comb_sig0_max=max(table_full_comb_sig0);
    table_full_comb_sig0=[table_full_comb_sig0; table_full_comb_sig0_avrs; table_full_comb_sig0_min; table_full_comb_sig0_max];

    welf90=squeeze(dwU_full_an(2,2,51,:,:));
    
    welf97=squeeze(dwU_full_an(1,2,51,:,:));

    table_full_comb_sig1=[log(migr2(:,1)) welf97(:,1) welf90(:,1) log(migr2(:,2)) welf97(:,2) welf90(:,2) ];
    table_full_comb_sig1_avrs=mean(table_full_comb_sig1);
    table_full_comb_sig1_min=min(table_full_comb_sig1);
    table_full_comb_sig1_max=max(table_full_comb_sig1);
    table_full_comb_sig1=[table_full_comb_sig1; table_full_comb_sig1_avrs; table_full_comb_sig1_min; table_full_comb_sig1_max];
end


end

for ii=1:6
for jj=1:6
    result_VTstd(ii,jj)=std(result_VTsim(ii,jj,:));
    result_UTstd(ii,jj)=std(result_UTsim(ii,jj,:));    
end
end

for ii=1:6
for jj=1:4
    result_VUstd(ii,jj)=std(result_VUsim(ii,jj,:));
    
end
end

s(cas).result_UTbig([1:2:12],:)=result_UTsim(:,:,1)*100;
s(cas).result_UTbig([2:2:12],:)=result_UTstd(:,:)*100;

s(cas).result_VTbig([1:2:12],:)=result_VTsim(:,:,1)*100;
s(cas).result_VTbig([2:2:12],:)=result_VTstd(:,:)*100;

s(cas).result_VUbig([1:2:12],:)=result_VUsim(:,:,1)*100;
s(cas).result_VUbig([2:2:12],:)=result_VUstd(:,:)*100;

new_welfare2_sig0=[new_welfare2_sig0; s(cas).result_UTbig([49 51 61 63])];
new_welfare2_sig1=[new_welfare2_sig1; s(cas).result_UTbig([49 51 61 63]+6)];



end


welfare_raw([1 3 5],:)=[squeeze(hatU_an_simcase(bb,1,1,:,1,1)) squeeze(hatU_an_simcase(bb,2,1,:,1,1))];
welfare_raw([2 4 6],:)=[squeeze(std(hatU_an_simcase(bb,1,1,:,:,1),[],5)) squeeze(std(hatU_an_simcase(bb,2,1,:,:,1),[],5))];



%disp('Table 3 - Rows III')
Table3(3,:)=new_welfare2_sig0(1,:);
Table3(13,:)=new_welfare2_sig1(1,:);
toc