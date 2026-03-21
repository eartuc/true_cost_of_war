
clearvars -except Table1a Table1b Table2 Table3 TableA1 TableA2 TableA3 TableA4  
tic
load data/paramters_iv.mat
load data/data_processed.mat


rng(123)
dur_all=[1:50 1000];
DN=25;

%migration ratio
mig_r = IDP./MIG0;

%Assume that Luh and Don receive proportional IDP just to fill in zeros,
%this is not needed.

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

for di=[1 10 51] 
for bi=[1 2 3] 
for si=[1 2]

        bta=bta_all(bi);
        sigma=sigma_all(si);
        dur=dur_all(di);

        theta=result_theta(bi,si,3);

        NN=length(index_);

        weight0=Ltm1(index_,5)/sum(Ltm1(index_,5));
        weight1a=MIG0(index_,1)/sum(MIG0(index_,1));
        weight2a=MIG0(index_,2)/sum(MIG0(index_,2));

        %EQ WEIGHTS
        weight1=(1./length(weight1a));
        weight2=(1./length(weight2a));

        hatV_1 = log(repmat( [mig_r(24,1) mig_r(25,2)], [NN 1])) - log(mig_r(index_,1:2));

        hatU_2 = 0*log(repmat( [mig_r(24,1) mig_r(25,2)], [NN 1])) - log(mig_r(index_,1:2));
        
        hatV(bi,si,di,:)= sum(hatV_1(:,:).*[weight1 weight2]).*(1./theta);

        hatU(bi,si,di,:)= sum(hatU_2(:,:).*[weight1 weight2]).*(1./theta);

        tmp1= sort(hatU_2(:,1)).*(1./theta);
        tmp2= sort(hatU_2(:,2)).*(1./theta);

        hatUsl(bi,si,di,:)= [tmp1(2) tmp2(2)];

        hatU_full(bi,si,di,:,:)= (hatU_2(:,:)).*(1./theta);

        mul=(1-bta^dur)/(1-bta);


        hatV_an(bi,si,di,:)=hatV(bi,si,di,:)/mul;      
        hatU_an(bi,si,di,:)=hatU(bi,si,di,:)/mul;
        hatUsl_an(bi,si,di,:)=hatUsl(bi,si,di,:)/mul;

        hatU_full_an(bi,si,di,:,:)=hatU_full(bi,si,di,:,:)/mul;
        
        w0(bi,si,di,:)=mean(rwaget(24:25,:),2); 
        u0(bi,si,di,:)=UU(mean(rwaget(24:25,:),2),sigma); 

        w0_full(bi,si,di,1,:)=mean(rwaget(24:25,:),2); 
        u0_full(bi,si,di,1,:)=UU(mean(rwaget(24:25,:),2),sigma); 

        dwV_an(bi,si,di,:)=100*(invUU(hatV_an(bi,si,di,:)+u0(bi,si,di,:),sigma)./w0(bi,si,di,:)-1);
        dwU_an(bi,si,di,:)=100*(invUU(hatU_an(bi,si,di,:)+u0(bi,si,di,:),sigma)./w0(bi,si,di,:)-1);
        
        dwUsl_an(bi,si,di,:)=100*(invUU(hatUsl_an(bi,si,di,:)+u0(bi,si,di,:),sigma)./w0(bi,si,di,:)-1);
        
        dwU_full_an(bi,si,di,:,:)=100*(invUU(hatU_full_an(bi,si,di,:,:)+u0_full(bi,si,di,:,:),sigma)./w0_full(bi,si,di,:,:)-1);


end
end
end

hatU_an_simcase(:,:,:,:,sim,cas)=hatU_an(:,:,:,:);
hatUsl_an_simcase(:,:,:,:,sim,cas)=hatUsl_an(:,:,:,:);
hatV_an_simcase(:,:,:,:,sim,cas)=hatV_an(:,:,:,:);

dwV_an(3,:,2:end,:)=0;
dwU_an(3,:,2:end,:)=0;
dwUsl_an(3,:,2:end,:)=0;

bb=[1 2 3];
result_VU = [squeeze(hatV_an(bb,1,1,:)) squeeze(hatU_an(bb,1,1,:)); squeeze(hatV_an(bb,2,1,:)) squeeze(hatU_an(bb,2,1,:)) ];


result_UT=[];
result_VT=[];
result_UTsl=[];

for ii=[1 10 51]
    tmp1=[squeeze(dwU_an(bb,1,ii,:)); squeeze(dwU_an(bb,2,ii,:))]/100;
    result_UT=[result_UT tmp1];
end

for ii=[1 10 51]
    tmp1=[squeeze(dwV_an(bb,1,ii,:)); squeeze(dwV_an(bb,2,ii,:))]/100;
    result_VT=[result_VT tmp1];
end

for ii=[1 10 51]
    tmp1=[squeeze(dwUsl_an(bb,1,ii,:)); squeeze(dwUsl_an(bb,2,ii,:))]/100;
    result_UTsl=[result_UTsl tmp1];
end

result_VUsim(:,:,sim)=result_VU(:,:);
result_VTsim(:,:,sim)=result_VT(:,:);
result_UTsim(:,:,sim)=result_UT(:,:);
result_UTslsim(:,:,sim)=result_UTsl(:,:);


if sim==1 && cas==1

    welf90=squeeze(dwU_full_an(2,1,51,:,:));
    
    welf97=squeeze(dwU_full_an(1,1,51,:,:));

    migr2=mig_r(1:23,:);

    table_full_b90=[migr2(:,1) log(migr2(:,1)) welf90(:,1) migr2(:,2) log(migr2(:,2)) welf90(:,2)];
    
    table_full_b97=[migr2(:,1) log(migr2(:,1)) welf97(:,1) migr2(:,2) log(migr2(:,2)) welf97(:,2)];

    
    table_full_comb_sig0=[log(migr2(:,1)) welf97(:,1) welf90(:,1) log(migr2(:,2)) welf97(:,2) welf90(:,2) ];    
    table_full_comb_sig0_avrs=sum(table_full_comb_sig0.*[weight1 weight1 weight1 weight2 weight2 weight2]);
    table_full_comb_sig0_median=median(table_full_comb_sig0);
    table_full_comb_sig0_min=min(table_full_comb_sig0);
    table_full_comb_sig0_max=max(table_full_comb_sig0);

    tmp1=sort(table_full_comb_sig0);
    tmp1(1,:)=[];
    tmp1(end,:)=[];
    table_full_comb_sig0_2min=min(tmp1);
    table_full_comb_sig0_2max=max(tmp1);
    
    table_full_comb_sig0=[table_full_comb_sig0; table_full_comb_sig0_avrs; table_full_comb_sig0_min; table_full_comb_sig0_max];
    
    block1=[table_full_comb_sig0_avrs; table_full_comb_sig0_median; table_full_comb_sig0_min; table_full_comb_sig0_max; table_full_comb_sig0_2min; table_full_comb_sig0_2max];
    block1(:,[1 4])=[];

    welf90=squeeze(dwU_full_an(2,2,51,:,:));    
    welf97=squeeze(dwU_full_an(1,2,51,:,:));

    table_full_comb_sig1=[log(migr2(:,1)) welf97(:,1) welf90(:,1) log(migr2(:,2)) welf97(:,2) welf90(:,2) ];
    table_full_comb_sig1_avrs=sum(table_full_comb_sig1.*[weight1 weight1 weight1 weight2 weight2 weight2]);
    table_full_comb_sig1_median=median(table_full_comb_sig1);
    table_full_comb_sig1_min=min(table_full_comb_sig1);
    table_full_comb_sig1_max=max(table_full_comb_sig1);

    tmp1=sort(table_full_comb_sig1);
    tmp1(1,:)=[];
    tmp1(end,:)=[];
    table_full_comb_sig1_2min=min(tmp1);
    table_full_comb_sig1_2max=max(tmp1);
    
    table_full_comb_sig1=[table_full_comb_sig1; table_full_comb_sig1_avrs; table_full_comb_sig1_min; table_full_comb_sig1_max];

    block2=[table_full_comb_sig1_avrs; table_full_comb_sig1_median; table_full_comb_sig1_min; table_full_comb_sig1_max; table_full_comb_sig1_2min; table_full_comb_sig1_2max];
    block2(:,[1 4])=[];


end



end %end sim

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

s(cas).result_UTslbig([1:2:12],:)=result_UTslsim(:,:,1)*100;
s(cas).result_UTslbig([2:2:12],:)=0;

s(cas).result_VTbig([1:2:12],:)=result_VTsim(:,:,1)*100;
s(cas).result_VTbig([2:2:12],:)=result_VTstd(:,:)*100;

s(cas).result_VUbig([1:2:12],:)=result_VUsim(:,:,1)*100;
s(cas).result_VUbig([2:2:12],:)=result_VUstd(:,:)*100;

s(cas).result_UTnewbig([1:3:18],:)=result_UTsim(:,:,1)*100;
s(cas).result_UTnewbig([3:3:18],:)=result_UTstd(:,:)*100;
s(cas).result_UTnewbig([2:3:18],:)=result_UTslsim(:,:,1)*100;


new_welfare2_sig0=[new_welfare2_sig0; s(cas).result_UTbig([49 51 61 63])];
new_welfare2_sig1=[new_welfare2_sig1; s(cas).result_UTbig([49 51 61 63]+6)];

end


welfare_raw([1 3 5],:)=[squeeze(hatU_an_simcase(bb,1,1,:,1,1)) squeeze(hatU_an_simcase(bb,2,1,:,1,1))];
welfare_raw([2 4 6],:)=[squeeze(std(hatU_an_simcase(bb,1,1,:,:,1),[],5)) squeeze(std(hatU_an_simcase(bb,2,1,:,:,1),[],5))];

welfare_raw_new([1 4 7],:)=[squeeze(hatU_an_simcase(bb,1,1,:,1,1)) squeeze(hatU_an_simcase(bb,2,1,:,1,1))];
welfare_raw_new([2 5 8],:)=[squeeze(hatUsl_an_simcase(bb,1,1,:,1,1)) squeeze(hatUsl_an_simcase(bb,2,1,:,1,1))];
welfare_raw_new([3 6 9],:)=[squeeze(std(hatU_an_simcase(bb,1,1,:,:,1),[],5)) squeeze(std(hatU_an_simcase(bb,2,1,:,:,1),[],5))];


%disp('Table 1b')
Table1b=welfare_raw_new;

%disp('Table 2')
Table2=s(1).result_UTbig;

new_welfare=[block1; block2];

new_welfare2_sig0old=[new_welfare2_sig0; block1([4 6 3],:) ];
new_welfare2_sig1old=[new_welfare2_sig1; block2([4 6 3],:) ];

new_welfare2_sig0=[new_welfare2_sig0; block1([5 4 3],:); zeros(2,4)  ];
new_welfare2_sig1=[new_welfare2_sig1; block2([5 4 3],:); zeros(2,4) ];


%disp('Table 3')
%disp('Note: Missing rows from welfare_weighted and welfare_distweighted scrips')
new_welfare2_sig0=new_welfare2_sig0([7 1 10 11 9 2 3 4 5 6 ],:);
new_welfare2_sig1=new_welfare2_sig1([7 1 10 11 9 2 3 4 5 6 ],:);

Table3=[new_welfare2_sig0; new_welfare2_sig1];

%disp('Tables A3 and A4 - Big tables in the appendix')
TableA3=table_full_comb_sig0;
TableA4=table_full_comb_sig1;
toc