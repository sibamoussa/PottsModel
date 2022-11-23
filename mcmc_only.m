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


ethresholds = [0,-3,-6,-9,-12,-15,-18,-21,-24,-100];
    for thread = 1
            thread
            ethreshold= ethresholds(thread);
            seed = randseq(30,'alphabet','dna');
            seed_e(thread) = mcmc_energy_full(JJ,h,letters_to_numbers(seed));
            %restricting sampling to energies -13 and higher inside mcmcmf:
            [mcmc_samples_full,def,ef] = mcmcf(seed, 0.5, 100000,JJ,h,ethreshold);
            mcmc_samples_full=mcmc_samples_full';
            name=['motif_A/replicate3/mcmc_samples_full_' num2str(thread) '.mat']
            save(name,'mcmc_samples_full');
% 
%             name=['motif_A/mcmc_energies_full_' num2str(thread) '.mat']
%             save(name,'ef');
%             
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