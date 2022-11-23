clear all 
close all

training_matrix=readcell("training_data.csv")
test_matrix=readcell("training_data.csv");

%If you want to train on the *worst* performers, go to the function
%"top_treshold" and change to "ascend"

%Change the 0.25 for the top% treshold you want to train it on:

main(1,training_matrix,test_matrix);

 function main(cutoff, training_matrix,test_matrix)
      %  test_matrix_contacts=readcell("training_data.xlsx");


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
      %  test_matrix_2=readcell("/Users/siminegroup/Desktop/submissions_20072022/fullE/fullEE/multiple/27201_rightpairs_1.csv");
       % test_all_sequences_2=cell2mat(test_matrix_contacts(2:end,1));


        %All test array as numbers:
        test_Y=letters_to_numbers(test_all_sequences);
      %  test_Y_2=letters_to_numbers(test_all_sequences_2);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      %  test_energy_table=compute_energy(JJ_average,h,test_Y,test_all_sequences);
        test_energy_table=compute_energy(JJ_average,h,test_Y,test_all_sequences);
    %    test_energy_table_2=compute_energy(JJ_average,h,test_Y_2,test_all_sequences_2);
 
        test_energies=table2array(test_energy_table(:,2));
    %    test_energies_2=table2array(test_energy_table_2(:,2));
        
        
        
%         
%         writetable(test_energy_table,'motif_A/test_energies_all.csv')  
%         writetable(test_energy_table_2,'motif_A/test_energies_contacts.csv')


     %   xlswrite('test_energies.xlsx',table2array(test_energy_table));
       % writetable(test_energy_table,'test_energies.xlsx')  

        test_score=cell2mat(test_matrix(2:end,2));
        test_energies=table2array(test_energy_table(:,2));
        all_test_data=[test_energies,test_score];

        
%        test_score_2=cell2mat(test_matrix_contacts(2:end,2));
%         test_energies_2=table2array(test_energy_table_2(:,2));
%         all_test_data_2=[test_energies_2,test_score_2];
%         
        %train_energy_table=compute_energy(JJ_average,h,train_Y,train_all_sequences)
        train_energy_table=compute_energy(JJ_average,h,train_Y,train_all_sequences)

        train_score=cell2mat(training_matrix(2:end,2));
        train_energies=table2array(train_energy_table(:,2));
        writematrix(train_energies,'motif_A/train_energies.csv')
%          for ii=1:435
%                JJ=JJ_average(:,:,ii);
%               name=['JJ_' num2str(ii) '.mat']
%               save(name,'JJ')
%         end
         save('motif_A/h_out.mat', 'h');
         JJ=JJ_average
         save('motif_A/JJ_out.mat', 'JJ');
  %          save('JJ_all.mat','JJ')
% 
        

        energy_of_trained_data=letters_to_numbers(table2array(train_file_cut(:,1)));
      %  energy_training=compute_energy(JJ_average,h,energy_of_trained_data,table2array(train_file_cut(:,1)));
        energy_training=compute_energy(JJ_average,h,energy_of_trained_data,table2array(train_file_cut(:,1)));

        
        
        %%% MCMC Sampling 
        %find the lowest emergy for mcmc seed 
        [min_vals, min_idx] = min(train_energy_table{:,2})
       % seq=train_energy_table(min_idx,1);

        size_table=(size(test_energy_table{:,2},1))


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Figures:
        figure(2)
        train_fig=scatter(train_score,train_energies,"ro");
        xlabel("Nupack Contacts");
        ylabel("Potts Energy");
        Pcoeff=corrcoef(train_score,train_energies);
        title({("Training data"),("Trained on top " + cutoff*100 + "%"),("Pearson Coefficient: " + Pcoeff(2,1))}) ;
        name_train=("Training top " + (cutoff *100) + "%");
     %   saveas(train_fig,"train.png")

        figure(3)
        %highe = (find(test_energies<-38 & test_energies>-39));
        test_fig=scatter(test_energies,test_energies,"go",'filled');
        hold on
        %test_fig_2=scatter(test_score_2,test_energies_2,"ko");

        Pcoeff=corrcoef(test_score,test_energies);
        title({("Test data"),("Trained on top " + cutoff*100 + "%"),("Pearson Coefficient: " + Pcoeff(2,1))})  ;
        xlabel("Nupack Contacts");
        ylabel("Partial Potts Energy (single contact)");
        name_test=("Test top " + (cutoff *100) + "%");
      %  saveas(test_fig,"test.png")


        figure(4)
   %     scatter(test_score_2,test_energies_2,"bo");
    %    hold on 
        scatter(test_score,test_energies,"go",'filled');
        sizes=size(test_energies)
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
        title("Training(47 Original) and test (MCMC Seed #27201, t=0.01) overlapped")
        xlabel("Nupack Contacts");
        ylabel("Partial Potts Energy (single contact)");
        legend('full test set from sampler','4 contacts test set from sampler','original training set with 4 contacts')
       % legend('unfiltered mcmc generated dataset','filtered mcmc dataset (4contacts)','original training set')

%        % legend('test all contacts','test filtered for 4 contacts','train')
        box off
        hold on


        %%%%%%
        
        
 
        
   
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
   compute=0;
    for i=1:size(Y,2)-1
      for j=(i+1):size(Y,2)
        if any(contacts_of_interest==i) %&& any(contacts_of_interest==j)
          energy=energy+Jij(sequence(i),sequence(j),l);
          l=l+1;
          compute=compute+1;
        end
        %energy=energy+Jij(sequence(j),sequence(i),i,j);
      end
    end
    for i=1:size(Y,2)
      energy=energy+hi(sequence(i),i);
      compute=compute+1;
    end
    energy=energy*-1;
    all_energies=[all_energies;energy];
  end
  energy_table=table(file,all_energies);
end
      
