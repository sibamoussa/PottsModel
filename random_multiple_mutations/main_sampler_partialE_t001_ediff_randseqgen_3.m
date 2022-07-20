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
          
        % MCMC Loop 
        for i=1:400:(size_table-97)
            seq=test_energy_table(i,1);
            seq_letters=table2array(seq); 
            seq=letters_to_numbers(seq_letters);
            [new_seq_list,energy_list,e_diff_accepted,e_diff_rejected]=mcmc(seq,JJ_average,h); %,seq_letters);
            
            %accepted sequences
            new_seqs=numbers_to_letters(new_seq_list);
            new_seqs_1=horzcat(new_seqs{:});
            new_seqs=reshape((new_seqs_1), [],30);
            final_array=table(new_seqs,energy_list);
            
            %accepted&rejectedenergydiff
             e_diff_accepted=num2cell(e_diff_accepted);       
             e_diff_rejected=num2cell(e_diff_rejected);       
           %  u_list=num2cell(u_list);       

           %  writecell(u_list,"partialEE_original/energies_"+ int2str(i) + ".xls",'Sheet',1,'Range','A1')
             writecell(e_diff_accepted,"partialEE/multiple/energies_"+ int2str(i) + ".xlsx",'Sheet',1,'Range','B1');
             writecell(e_diff_rejected,"partialEE/multiple/energies_"+ int2str(i) + ".xlsx",'Sheet',1,'Range','C1');

           % xlswrite("partialEE_original/post_sampling_mcmc_partial_energy_t001_" + int2str(i) + ".csv",final_array,'Sheet1'); 
            writetable(final_array,("partialEE/multiple/post_sampling_mcmc_partial_energy_t001_" + int2str(i) + ".csv")); 
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

        figure(3)
        test_fig=scatter(test_score,test_energies,"go",'filled');
        hold on
 
        Pcoeff=corrcoef(test_score,test_energies);
        title({("Test data"),("Trained on top " + cutoff*100 + "%"),("Pearson Coefficient: " + Pcoeff(2,1))})  ;
        xlabel("Nupack Contacts");
        ylabel("Potts Energy");
        name_test=("Test top " + (cutoff *100) + "%");
 

        figure(4)
         hold on 
        scatter(test_score,test_energies,"go",'filled');
        hold on
        scatter(train_score,train_energies,"rs",'filled');

        [~,idx] = sort(train_score(:));
        [rows,cols]=size(train_score);
        pts=(rows+1)-(1-cutoff)*rows;
 
        title("Training and test overlapped")
        xlabel("Nupack Contacts");
        ylabel("Potts Energy");
        legend('full test set from sampler','70 good test points from sampler','original training set with 4 contacts')
 
        box off
        hold on


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
    function [energy]=compute_energy_seq(Jij,hi,Y);
        all_energies=[];
        file=numbers_to_letters(Y);
        contacts_of_interest=[1 2 3 4 27 28 29 30];
        for index=1:size(Y,1);
          sequence=Y(index,:);
          energy=0;
          l=1;
            for i=1:size(Y,2)-1;
                for j=(i+1):size(Y,2);
                    if any(contacts_of_interest==i) %&& any(contacts_of_interest==j)
                        energy=energy+Jij(sequence(i),sequence(j),l);
                        l=l+1;
                    end 
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


    %random mutation generation
    function perm_mutated_seq=mutate_any_base(seq);
            new_seq=seq;
            num_mutations=randi([1,8]);
            positions=randperm(30,num_mutations);
            mutation_type=randi([1,4],1,num_mutations);
            new_seq(positions)=mutation_type;
            perm_mutated_seq=new_seq;
         
    end
    
    function [perm_mutated_seq,mutated_seq]=check_mut_seq(perm_mutated_seq,seq)
             if sum(perm_mutated_seq-seq)~=0
                mutated_seq=perm_mutated_seq;
             else;
                perm_mutated_seq=mutate_any_base(seq);
             end 
     end 
%      
    %probability of accepting new val
    function y=prob_density(seq,beta,JJ,hi)
        y=exp((-beta*compute_energy_seq(JJ,hi,seq)));  % need to work on this here 
    end 
   
    %returns array with new sequences 
    function [new_seq_list,energy_list,e_diff_accepted,e_diff_rejected]=mcmc(seq,JJ,hi) %,seq_letters);
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
        
          %e_difference
          reject_counter= 1;
          e_diff_accepted=zeros(n_iterations,1);
          e_diff_rejected=zeros(n_iterations,1);
          
        % list mutated seq here 
           for i=1:n_iterations;
             perm_mutated_seq=mutate_any_base(seq);
             mutated_seq=check_mut_seq(perm_mutated_seq,seq);
           % perm_mutated_seq=mutate_any_base(seq)
             e_diff=prob_density(mutated_seq,beta,JJ,hi)/prob_density(seq,beta,JJ,hi);
             e_diff_actual=log(e_diff)/(-beta);
             u=rand(1);

             if u <= min(1,e_diff);
                 seq=mutated_seq;
                 energy=compute_energy_seq(JJ,hi,seq);
                 new_seq_list(counter,:)=seq;
                 energy_list(counter,1)=energy;
                 e_diff_accepted(counter,1)=e_diff_actual;
                 counter=counter+1;
             else
                 e_diff_rejected(reject_counter,1)=e_diff_actual;
                 reject_counter=reject_counter+1; 

             end 
           end 
          energy_list(any( energy_list==0, 2))=[];

    end 