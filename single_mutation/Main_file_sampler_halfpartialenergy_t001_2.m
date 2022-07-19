clear all 
close all
training_matrix=readcell("originaldataset_filtered_nupack_50K.csv");
test_matrix=readcell("mcmc_sampler.csv");

%If you want to train on the *worst* performers, go to the function
%"top_treshold" and change to "ascend"

%Change the 0.25 for the top% treshold you want to train it on:

main(1,training_matrix,test_matrix);

 function main(cutoff, training_matrix,test_matrix)
        

        %cutoff for the graph:
        cutoff_graph=cutoff;
        %frequency theta:
        theta=0.1;
        %number of states:
        q=4;

        %Top % training array as sequences + their score:
        train_file_cut=creating_matrix(training_matrix,cutoff);

        %Using train_file_cut as the algorithm input (sequences): 
        train_fileplm=table2array(train_file_cut(:,1));
        [JJ,h,JJ_average]=plmDCA_tresholdtest(train_fileplm,theta,cutoff);

    

        %All train array as sequences:
        train_all_sequences=cell2mat(training_matrix(2:end,1));

        %All train array as numbers:
        train_Y=letters_to_numbers(train_all_sequences);


        %All test array as sequences:
        test_all_sequences=cell2mat(test_matrix(2:end,1));
        
        %adding second test matrix
       % test_matrix_2=readcell("mcmc_sampler.csv");
        %test_all_sequences_2=cell2mat(test_matrix_2(2:end,1));


        %All test array as numbers:
        test_Y=letters_to_numbers(test_all_sequences);
        %test_Y_2=letters_to_numbers(test_all_sequences_2);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %test_energy_table=compute_energy(JJ_average,h,test_Y,test_all_sequences);
        test_energy_table=compute_energy_motif(JJ_average,h,test_Y,test_all_sequences);
        %sibtest_energy_table_2=compute_energy_motif(JJ_average,h,test_Y_2,test_all_sequences_2);
 
        test_energies=table2array(test_energy_table(:,2));
       % test_energies_2=table2array(test_energy_table_2(:,2));

       % xlswrite('test_energies.xlsx',table2array(test_energy_table));
       % writetable(test_energy_table,'test_energies.xlsx')  

        test_score=cell2mat(test_matrix(2:end,2));
        test_energies=table2array(test_energy_table(:,2));
        all_test_data=[test_energies,test_score];

       %test_score_2=cell2mat(test_matrix_2(2:end,2));
       % test_energies_2=table2array(test_energy_table_2(:,2));
       % all_test_data_2=[test_energies_2,test_score_2];
        
        %train_energy_table=compute_energy(JJ_average,h,train_Y,train_all_sequences)
        train_energy_table=compute_energy_motif(JJ_average,h,train_Y,train_all_sequences)

        train_score=cell2mat(training_matrix(2:end,2));
        train_energies=table2array(train_energy_table(:,2));
        
        energy_of_trained_data=letters_to_numbers(table2array(train_file_cut(:,1)));
      % energy_training=compute_energy(JJ_average,h,energy_of_trained_data,table2array(train_file_cut(:,1)));
        energy_training=compute_energy_motif(JJ_average,h,energy_of_trained_data,table2array(train_file_cut(:,1)));
 
        %%% MCMC Sampling 
        %find the lowest emergy for mcmc seed 
        [min_vals, min_idx] = min(train_energy_table{:,2});
       % seq=train_energy_table(min_idx,1);

        size_table=(size(test_energy_table{:,2},1));
          
        for i=1:400:(size_table-97)
            seq=test_energy_table(i,1);
            seq_letters=table2array(seq); 
            seq=letters_to_numbers(seq_letters);
            [new_seq_list,energy_list]=mcmc(seq,JJ_average,h,seq_letters);
            new_seqs=numbers_to_letters(new_seq_list);
            new_seqs_1=horzcat(new_seqs{:});
            new_seqs=reshape((new_seqs_1), [],30);
            final_array=table(new_seqs,energy_list);
            %xlswrite(("partialE/post_sampling_mcmc_partial_energy" + int2str(i) + ".csv"),final_array); 
            writetable(final_array,("partialEE_original/post_sampling_mcmc_partial_energy_t001_" + int2str(i) + ".csv")); 
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Figures:
        figure(2)
        train_fig=scatter(train_score,train_energies,"ro");
        xlabel("Nupack Contacts");
        ylabel("Potts Energy");
        Pcoeff=corrcoef(train_score,train_energies);
        title({("Training data"),("Trained on top " + cutoff*100 + "%"),("Pearson Coefficient: " + Pcoeff(2,1))}) ;
        name_train=("Training top " + (cutoff *100) + "%");
     %  saveas(train_fig,"train.png")

        figure(3)
        %highe = (find(test_energies<-38 & test_energies>-39));
        test_fig=scatter(test_score,test_energies,"go",'filled');
        hold on
        %test_fig_2=scatter(test_score_2,test_energies_2,"ko");

        Pcoeff=corrcoef(test_score,test_energies);
        title({("Test data"),("Trained on top " + cutoff*100 + "%"),("Pearson Coefficient: " + Pcoeff(2,1))})  ;
        xlabel("Nupack Contacts");
        ylabel("Potts Energy");
        name_test=("Test top " + (cutoff *100) + "%");
      %  saveas(test_fig,"test.png")


        figure(4)
        %scatter(test_score_2,test_energies_2,"bo");
        hold on 
        scatter(test_score,test_energies,"go",'filled');
        hold on
        scatter(train_score,train_energies,"rs",'filled');

        [~,idx] = sort(train_score(:));
        [rows,cols]=size(train_score);
        pts=(rows+1)-(1-cutoff)*rows;
       % scatter(train_file_cut{2:end,2},train_file_cut{2:end,2},'ro','filled')
       % hold on
        %scatter(train_score(idx(1:pts)), train_energies(idx(1:pts)),'ko','filled')
      %  scatter(train_score((1:pts-1)), train_energies((1:pts-1)),'bo','filled')
             
%         %x(idx(1:pts))
        title("Training and test overlapped")
        xlabel("Nupack Contacts");
        ylabel("Potts Energy");
        legend('full test set from sampler','70 good test points from sampler','original training set with 4 contacts')
       % legend('unfiltered mcmc generated dataset','filtered mcmc dataset (4contacts)','original training set')

%        % legend('test all contacts','test filtered for 4 contacts','train')
        box off
        hold on

        
%         %find the mean of values for training 
%         size(train_energies((1:pts)))
%         avg_val = mean(train_energies((1:pts)))
%         train_energies_cut=train_energies((1:pts))
          
%         accuracy=size(train_energies_cut(rowstokeep))/size(train_energies)*100
%       
%         
        
        
%         avg_val = mean(train_energies(idx(1:pts)))
%         train_energies=train_energies(idx(1:pts))
%         rowstokeep=train_energies<avg_val
%         accuracy=size(train_energies(rowstokeep))/size(train_energies)*100
        
        %%%%%%
        
        
        %For the sequences from 1 to "size_to_eval", how many are in common
        %between the scores and energies:
        percent_eval=0.16
        size_to_eval=ceil(size(test_energies,1)*percent_eval);
        %energy_sorted=sortrows(all_test_data,1);
       % score_sorted=sortrows(all_test_data,2) %,"descend");
%         energy_sorted=all_test_data(:,1);
%         score_sorted=all_test_data(:,2); %,"descend");
%         int=intersect(energy_sorted(1:size_to_eval,2),score_sorted(1:size_to_eval,2));
%         percent_int=100*size(int,1)/size_to_eval    

        % figure(5)
        % idx=kmeans(test_energies,4);
        % scatter(test_energies,idx)
        
   
    end

    %Computes energy for every sequence:
       function [energy_table]=compute_energy(Jij,hi,Y,file)
        all_energies=[];

        for index=1:size(Y,1)
          sequence=Y(index,:);
          energy=0;
          l=1;
            for i=1:size(Y,2)-1
                for j=(i+1):size(Y,2)
                    energy=energy+Jij(sequence(i),sequence(j),l);
                    l=l+1;
                    %energy=energy+Jij(sequence(j),sequence(i),i,j);
                end
            end

            for i=1:size(Y,2)
                energy=energy+hi(sequence(i),i);
            end
            energy=energy*-1;
            all_energies=[all_energies;energy];

        end

        energy_table=table(file,all_energies);   
    end

    %Returns matrix that transformed n.a to their respective numbers:
    function [SeqM,L]=letters_to_numbers(CompleteSeq)

        %Creating the matrix SeqM that has the same size as matrix CompleteSeq
        SeqM=zeros(size(CompleteSeq));

        %Going over all letters in the "CompleteSeq" matrix, and assigning them a number: 
        for i=1:(size(CompleteSeq,2)) 
            for j=1:(size(CompleteSeq,1))
                if CompleteSeq(j,i)=="A"
                    SeqM(j,i)=1; 
                elseif CompleteSeq(j,i)=="T"
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

    %Returns a table with only the top% of sequences:
    function top_values_table = creating_matrix(full_file,percent_top)
     %assuming title to excel sheet so starts at 2:end

        sequences=mutated_sequences(cell2mat(full_file(2:end,1)));
        foldmean=cell2mat(full_file(2:end,2));

        all_values=table(sequences,foldmean);

        top_values_table= top_treshold(percent_top,all_values) ;
        %top_values_table=values_under(fold1,fold2,all_values);

    end

    function top_values_table=top_treshold(percent_top,all_values)

        %Change "descend" to "ascend" if taking the poorest performers:
        sorted_values =all_values % sortrows(all_values,2,"ascend");
        top_values_table= sorted_values(1:ceil(size(all_values,1)*percent_top),:)

    end

    function differences= mutated_sequences(CompleteSeq)

        differences=[];

        %Indexing column numbers that have rows with different bases:

        for i=1:size(CompleteSeq,2) 
           if size(unique(CompleteSeq(:,i)),1) ~=1 
               differences=[differences, CompleteSeq(:,i)];
           end     
        end
    end
    
    %returns energy for first and last 4 sections of sequence
    function [energy_table]=compute_energy_motif(Jij,hi,Y,file)
  all_energies=[];
    % Motif Energies
  contacts_of_interest=[1 2 3 4 27 28 29 30];
  for index=1:size(Y,1)
   sequence=Y(index,:);
   energy=0;
   l=1;
    for i=1:size(Y,2)-1
      for j=(i+1):size(Y,2)
        if any(contacts_of_interest==i) %&& any(contacts_of_interest==j)
          energy=energy+Jij(sequence(i),sequence(j),l);
          l=l+1;
        end
        %energy=energy+Jij(sequence(j),sequence(i),i,j);
      end
    end
    for i=1:size(Y,2)
      energy=energy+hi(sequence(i),i);
    end
    energy=energy*-1;
    all_energies=[all_energies;energy];
  end
  energy_table=table(file,all_energies);
end
     
    %%%%%%%%%%%%%%%%%%%%%%%
    
    %Converts numbers to letters 
    function [SeqM_2]=numbers_to_letters(CompleteSeq)

        %Creating the matrix SeqM that has the same size as matrix CompleteSeq
        SeqM_2=(strings(size(CompleteSeq)));
        %Going over all letters in the "CompleteSeq" matrix, and assigning them a number: 
        for i=1:(size(CompleteSeq,2)) 
            for j=1:(size(CompleteSeq,1))
                CompleteSeq(j,i);
                if CompleteSeq(j,i)==1;
                    SeqM_2(j,i)="A"; 
                elseif CompleteSeq(j,i)==2;
                    SeqM_2(j,i)="T";
                elseif CompleteSeq(j,i)==3;
                    SeqM_2(j,i)="C";
                elseif CompleteSeq(j,i)==4;
                    SeqM_2(j,i)="G";
                else
                    SeqM_2(j,i)=''; %if '-' is part of it change to 5     
                end  
            end
        end   

        %Taking the length L (number of columns) of the matrix:

    end

     %Computes energy for one sequence at time 
    function [energy]=compute_energy_seq(Jij,hi,Y,file)
        all_energies=[];

        for index=1:size(Y,1);
          sequence=Y(:);
          energy=0;
          l=1;
            for i=1:size(Y,2)-1;
                for j=(i+1):size(Y,2);
                    energy=energy+Jij(sequence(i),sequence(j),l);
                    l=l+1;
                end
            end

            for i=1:size(Y,2)
                energy=energy+hi(sequence(i),i);
            end
            energy=energy*-1;
            all_energies=[all_energies;energy];

        end

        energy_table=[file,all_energies];   
    end   

    %generate a single mutation 
    function mutated_seq=mutate_seqs(seq)  ;
        perm=randi(4);
        loc=randi(30);
        seq(loc)=perm;
        mutated_seq=seq;
    end 
   
    %returns energy for first and last 4 sections of sequence for each seq
    function [energy]=compute_energy_motif_seq(Jij,hi,Y,file)
  all_energies=[];
    % Motif Energies
  contacts_of_interest=[1 2 3 4 27 28 29 30];
  for index=1:size(Y,1)
   sequence=Y(index,:);
   energy=0;
   l=1;
    for i=1:size(Y,2)-1
      for j=(i+1):size(Y,2)
        if any(contacts_of_interest==i) %&& any(contacts_of_interest==j)
          energy=energy+Jij(sequence(i),sequence(j),l);
          l=l+1;
        end
        %energy=energy+Jij(sequence(j),sequence(i),i,j);
      end
    end
    for i=1:size(Y,2)
      energy=energy+hi(sequence(i),i);
    end
    energy=energy*-1;
    all_energies=[all_energies;energy];
  end
  energy_table=[file,all_energies];
    end   
    %probability of accepting new val
    function y=prob_density(seq,beta,JJ,hi,seq_letters)
        y=exp((-beta*compute_energy_motif_seq(JJ,hi,seq,seq_letters)));  % need to work on this here 
    end 
   
    %returns array with new sequences 
    function [new_seq_list,energy_list]=mcmc(seq,JJ,hi,seq_letters);
        n_iterations = 300000;
        k=8.617*10^-5;
        deg_C=0.01;
        constant=1;
        t=273+deg_C;
        barrier_height=k*t;
        beta=constant/barrier_height;
        counter = 1 ;
        new_seq_list=zeros(n_iterations,30);
        energy_list=zeros(n_iterations,1);
        % list mutated seq here 
           for i=1:n_iterations;
             mutated_seq=mutate_seqs(seq);
             if rand(1) < min(1,prob_density(mutated_seq,beta,JJ,hi,seq_letters)/prob_density(seq,beta,JJ,hi,seq_letters));
                 seq=mutated_seq;
                 energy=compute_energy_motif_seq(JJ,hi,seq,seq_letters);
                 new_seq_list(counter,:)=seq;
                 energy_list(counter,1)=compute_energy_motif_seq(JJ,hi,seq,seq_letters);
                 counter=counter+1;
             end 
           end 
          energy_list(any( energy_list==0, 2))=[];
    end 
       
    