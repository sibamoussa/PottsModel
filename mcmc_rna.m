%mcmc only, no training
load model_motifA.mat
load motifA_dataset.mat
training_matrix=trset(1:10:500,:);
tr=letters_to_numbers(training_matrix);

% Onehot the training set
for ii=1:50
    A = (tr(ii,:)'==1:4)';
    tr_onehot(ii,:)=A(:)';
end


ethresholds = [0,-3,-6,-9,-12,-15,-18,-21,-24,-100];
    for thread = 1:10
            thread
            ethreshold= ethresholds(thread);
            seed = randseq(30,'alphabet','rna');
            seed_e(thread) = mcmc_energy_full(JJ,h,letters_to_numbers(seed));
            %restricting sampling to energies -13 and higher inside mcmcmf:
            [mcmc_samples_full,corr_cont_f,def,ef] = mcmcf(seed, 0.5, 50000,JJ,h,ethreshold);
            [~,ascind]=sort(corr_cont_f,1);

            asc = corr_cont_f(ascind,:);
            samples(thread,:,:)=asc(1:100,:);
            rangede=ef(ascind);
            ranked_e(thread) = rangede(100,1);
            sample_seq_temp=mcmc_samples_full(:,ascind);
            sample_seq(thread,:,:) = sample_seq_temp(:,1:100);
            success_rate(thread)=sum(samples(thread,samples(thread,:,2)==4,2))/400;  
            %diversity analysis
            for kk=1:100
                %one-hot encode the samples
                A = (sample_seq(thread,:,kk)'==1:4)';
                flatt=A(:)';
                %calculate cosine distance to training set
                for rr=1:50
                    cos_d(rr) = 1-dot(double(flatt),double(tr_onehot(rr,:)))/norm(double(flatt))/norm(double(tr_onehot(rr,:)));
                end
                distance(thread,kk)=min(cos_d);
            end
            
    end
  %      [mcmc_samples_partial, corr_cont_p,dep,ep] = mcmcp(seed, 0.1,10000, JJ_average,h,contacts_of_interest);
        save('mcmc_full_threshold_scan','ascind','sample_seq','mcmc_samples_full',"corr_cont_f","def","ef","seed_e","samples","success_rate","distance","ethresholds");

%let's assess the diversity of top 100 samples
% sample_seq(1-10,1-30,1-100) are the top 100 sequences 
% training_matrix is the training set: letters








        
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
    function [mcmc_samples_full,correct_contacts_f,de,e]= mcmcf(seed,tempr, steps,Jij,hi,e_threshold)
  %  correct_contacts_f = [];
        current_seq = letters_to_numbers(seed);
        accepted_count=0;
        de = 0;
        e=0;
       % print('running mcmc full')
        mcmc_samples_full=[];
        e1 = mcmc_energy_full(Jij,hi,current_seq);

        for ii=1:steps
            %seq=letters_to_numbers(randseq(length(seed),'alphabet','rna'));
            seq = mutate_two_bases(current_seq);
            e2 = mcmc_energy_full(Jij,hi,seq);
            ar = expz(-(e2-e1)/tempr);
            if rand < ar &&e2>= e_threshold % -12.7706
                accepted_count= accepted_count +1;
                de(accepted_count) = e2-e1;
                e(accepted_count)=e2;
                e1=e2;
                mcmc_samples_full = [mcmc_samples_full, seq'];
                current_seq= seq;
                phe_str=rnafold(numbers_to_letters(seq));
                map=contactmap(phe_str,length(seq));
                if map(1,30)==1 && map(2,29)==1 && map(3,28)==1 && map(4,27)==1
                    correct_contacts_f(accepted_count,:)=[e2,4];
                elseif (map(1,30)==1 && map(2,29)==1 && map(3,28)==1) || (map(2,29)==1 && map(3,28)==1 && map(4,27)==1)
                    correct_contacts_f(accepted_count,:)=[e2,3];
                elseif (map(1,30)==1 && map(2,29)==1) || (map(2,29)==1 && map(3,28)==1)|| (map(3,28)==1 && map(4,27)==1)
                    correct_contacts_f(accepted_count,:)=[e2,2];
                elseif map(1,30)==1 || map(2,29)==1 || map(3,28)==1 || map(4,27)==1
                    correct_contacts_f(accepted_count,:)=[e2,1];
                else 
                    correct_contacts_f(accepted_count,:)=[e2,0];
                end
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
                    SeqM_2(j,i)='U';
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
                elseif CompleteSeq(j,i)=='U'
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