
%mcmc only, no training
clear all
close all

%load model_motifA.mat
%load motifA_dataset.mat
%load random_rna.mat
load motif_A/motifA.mat
trset=readcell("motifA_trainingset.csv")

%trset=readcell("/Users/siminegroup/Desktop/Siba/Potts_Model_Exploration/oracle_sampling/oracle/motifB_trainingset.csv");
training_matrix=trset(2:end,1);
training_matrix=cell2mat(training_matrix);
tr=letters_to_numbers(training_matrix);

% Onehot the training set
for ii=1:50
    A = (tr(ii,:)'==1:4)';
    tr_onehot(ii,:)=A(:)';
end

samples=[]
ethresholds = [0,-3,-6,-9,-12,-15,-18,-21,-24,-100];
    for thread = 9
            thread
            ethreshold= ethresholds(thread);
            seed = randseq(30,'alphabet','dna');
            seed_e(thread) = mcmc_energy_full(JJ,h,letters_to_numbers(seed));
            %restricting sampling to energies -13 and higher inside mcmcmf:

            %load previously saved mcmc samples 
            mcmc_filename=['/Users/siminegroup/Desktop/Siba/Potts_Model_Exploration/motifA/mcmc_samples/replicate 2/samples/motifAcontacts_motifA_' num2str(thread) '.csv'];       
            corr_cont_filename=['/Users/siminegroup/Desktop/Siba/Potts_Model_Exploration/motifA/mcmc_samples/replicate 2/samples/rightcontacts/motifA_rightcontacts__motifAcontacts_motifA_' num2str(thread) '.csv'];
            mcmc_samples_f=readcell(mcmc_filename);
            corr_cont=readcell(corr_cont_filename);

            mcmc_samples_full=(mcmc_samples_f(2:end,2));
            corr_cont=(corr_cont(2:end,2));
            
            seq_letters=cell2mat(corr_cont);
            corr_cont_l=letters_to_numbers(seq_letters);
            
            seq_letters_all=cell2mat(mcmc_samples_full);
            mcmc_l=letters_to_numbers(seq_letters_all);
            [corr_cont_f] = mcmc_energy_full_all(JJ,h,corr_cont_l);
            [mcmc_all] = mcmc_energy_full_all(JJ,h,mcmc_l);

            %sort by increasing energies of accepted sequences
            [~,ascind]=sort(table2array(mcmc_all));
           
            %order sample list in order of energies 
            asc = mcmc_samples_full(ascind,:);
            sample_tot=asc(1:100,:)
            %order sequences with correct contacts
           
            %calculate sucess rate
            [tf,idx] = ismember(corr_cont,sample_tot);
            total_corr=sum(tf)
            success_rate(thread)=(total_corr/100);

             %diversity analysis
             for kk=1:100
%                 %one-hot encode the samples
                  onehot=sample_tot(kk,:);
                  seq_letters=cell2mat(onehot);
                  corr_cont_l=letters_to_numbers(seq_letters);
                  A = (corr_cont_l'==1:4)'
                  flatt=A(:)';
%                 %calculate cosine distance to training set
                 for rr=1:50
                     cos_d(rr) = 1-dot(double(flatt),double(tr_onehot(rr,:)))/norm(double(flatt))/norm(double(tr_onehot(rr,:)));
                 end
                 distance(thread,kk)=min(cos_d);
             end
            dist(thread)=mean(distance(thread,:));
    end
       %  [mcmc_samples_partial, corr_cont_p,dep,ep] = mcmcp(seed, 0.1,10000, JJ_average,h,contacts_of_interest);
         save('mcmc_full_threshold_scan','ascind','mcmc_samples_full',"corr_cont_f",'sample_tot') %,"def","ef","seed_e","samples","success_rate","distance","ethresholds");
%          save("replicates/dist_motifA_rep3",'dist')
%          save("replicates/success_rate_motifA_rep3",'success_rate')
    %let's assess the diversity of top 100 samples
% sample_seq(1-10,1-30,1-100) are the top 100 sequences 
% training_matrix is the training set: letters

    %Computes energy for every sequence:


function [mcmcef] = mcmc_energy_full_all(Jij,hi,seqq) 
  all_energies=[];

    for index=1:size(seqq,1)
        seq=seqq(index,:);
        energy=0;
          l=1;
            for i=1:size(seqq,2)-1
                for j=(i+1):size(seqq,2)
                    energy=energy+Jij(seq(i),seq(j),l);
                    l=l+1;
                    %energy=energy+Jij(sequence(j),sequence(i),i,j);
                end
            end

            for i=1:size(seqq,2)
                energy=energy+hi(seq(i),i);
            end
            mcmcef=energy*-1;
            all_energies=[all_energies;mcmcef];

    end
  mcmcef=table(all_energies);

end 


        
function [mcmcef] = mcmc_energy_full(Jij,hi,seq) 
        energy=0;
          l=1;
            for i=1:size(seq,2)-1
                for j=(i+1):size(seq,2)
                    energy=energy+Jij(seq(i),seq(j),l);
                    l=l+1;
                    %energy=energy+Jij(sequence(j),sequence(i),i,j);
                end
            end

            for i=1:size(seq,2)
                energy=energy+hi(seq(i),i);
            end
            mcmcef=energy*-1;
    end
        
function [map] = contactmap(bd, SequenceLength)
map=zeros(SequenceLength,SequenceLength);
for ii=1:SequenceLength
    if bd(ii)=='('
        dist=0;
        for jj=ii+1:SequenceLength
            if bd(jj)=='('
                dist=dist+1;
            elseif bd(jj)==')'
                dist=dist-1;
            end
            if dist==-1
                map(ii,jj)=1;
                map(jj,ii)=1;
                break
            end
        end
    end
end
end


    function new_seq = mutate_one_base(old_seq)
        new_seq = old_seq;
        count=0;
        while sum(new_seq - old_seq) == 0
            replacement=randi([1,4]);
            position=randi([1,30]);
            new_seq(position)=replacement;
            count=count+1;
        end
    end
    function new_seq = mutate_two_bases(old_seq)
        new_seq = old_seq;
        count=0;
        while sum(new_seq - old_seq) == 0
            replacement=randi([1,4],1,2);
            position=randi([2,30],1,2);
            new_seq(position(1))=replacement(1);
            new_seq(position(2))=replacement(2);
            count=count+1;
        end
    end
    function [mcmc_samples_full,de,e]= mcmcf(seed,tempr, steps,Jij,hi,e_threshold)
  %  correct_contacts_f = [];
        current_seq = letters_to_numbers(seed);
        accepted_count=0;
        de = 0;
        e=0;
       % print('running mcmc full')
        mcmc_samples_full=[];
        e1 = mcmc_energy_full(Jij,hi,current_seq);
       % sequences=[];

        for ii=1:steps
            %seq=letters_to_numbers(randseq(length(seed),'alphabet','rna'));
            seq = mutate_two_bases(current_seq);
            e2 = mcmc_energy_full(Jij,hi,seq);
            ar = exp(-(e2-e1)/tempr);
            if rand < ar &&e2>= e_threshold % -12.7706
                accepted_count= accepted_count +1;
                de(accepted_count) = e2-e1;
                e(accepted_count)=e2;
                e1=e2;
                mcmc_samples_full = [mcmc_samples_full, seq'];
                current_seq= seq;
                
               % sequences(ii)=current_seq
%                 phe_str=rnafold(numbers_to_letters(seq))
%                 map=contactmap(phe_str,length(seq));
%                 if map(8,23)==1 && map(9,22)==1 && map(10,21)==1
%                     correct_contacts_f(accepted_count,:)=[e2,3];
%                 elseif (map(8,23)==1 && map(9,22)==1 )|| (map(9,22)==1 && map(10,21)==1) || (map(8,23)==1 && map(10,21)==1)
%                     correct_contacts_f(accepted_count,:)=[e2,2];
%                 elseif (map(8,23)==1) || (map(9,22)==1) || (map(10,21)==1)
%                     correct_contacts_f(accepted_count,:)=[e2,1];
%                 else 
%                     correct_contacts_f(accepted_count,:)=[e2,0];
%                 end
             end

        end
    end
    
    
     function [SeqM_2]=numbers_to_letters(CompleteSeq)

        %Creating the matrix SeqM that has the same size as matrix CompleteSeq
        SeqM_2=(char(size(CompleteSeq)));
        %Going over all letters in the "CompleteSeq" matrix, and assigning them a number: 
        for i=1:(size(CompleteSeq,2)) 
            for j=1:(size(CompleteSeq,1))
                CompleteSeq(j,i);
                if CompleteSeq(j,i)==1;
                    SeqM_2(j,i)='A'; 
                elseif CompleteSeq(j,i)==2;
                    SeqM_2(j,i)='T';
                elseif CompleteSeq(j,i)==3;
                    SeqM_2(j,i)='C';
                elseif CompleteSeq(j,i)==4;
                    SeqM_2(j,i)='G';
                else
                    SeqM_2(j,i)=''; %if '-' is part of it change to 5     
                end  
            end
        end   

        %Taking the length L (number of columns) of the matrix:

     end
    
     function [SeqM,L]=letters_to_numbers(CompleteSeq)

        %Creating the matrix SeqM that has the same size as matrix CompleteSeq
        SeqM=zeros(size(CompleteSeq));

        %Going over all letters in the "CompleteSeq" matrix, and assigning them a number: 
        for i=1:(size(CompleteSeq,2)) 
            for j=1:(size(CompleteSeq,1))
                if CompleteSeq(j,i)=='A'
                    SeqM(j,i)=1; 
                elseif CompleteSeq(j,i)=='T'
                    SeqM(j,i)=2;
                elseif CompleteSeq(j,i)=='C'
                    SeqM(j,i)=3;
                elseif CompleteSeq(j,i)=='G'
                    SeqM(j,i)=4;
                else
                    SeqM(j,i)=1; %if '-' is part of it change to 5
                end       
            end
        end   

        %Taking the length L (number of columns) of the matrix:
        L=size(SeqM,2); 

    end