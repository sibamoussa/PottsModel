
%--Getting outside functions:--
% fprintf('Compiling g_rC.c\n');
% mex -outdir functions functions/g_rC.c
% 
% cd minFunc_2012
% addpath(genpath(pwd))
% mexAll
% cd ..

%--Main function:--
function [JJ,h,JJ_test]=plmDCA(file,reweighting_threshold,cutoff)
    
%If should-be numericals are passed as strings, convert them.
    if (isstr(reweighting_threshold))
        reweighting_threshold = str2num(reweighting_threshold);
    end

%Minimization options
    options.method='lbfgs'; %Minimization scheme. Default: 'lbfgs', 'cg' for conjugate gradient (use 'cg' if out of RAM).
    options.Display='off';
    options.progTol=1e-7; %Threshold for when to terminate the descent. Default: 1e-9. 
    addpath(genpath(pwd));
    
%Read inputfile (removing inserts), remove duplicate sequences, and calculate weights and B_eff.

    [N,B_with_id_seq,q,Y]=return_alignment(file);
    %N : Number of i positions in one sequence
    %B : Number of sequences
    %q : Values possible at one positions i (spins)
    %Y : Matrix containing the data
    %Y=unique(Y,'rows'); %not sure if interferes w/ results?
    [B,N]=size(Y);

    %Calculating 1/m_a, the weights of each sequence:
    weights = ones(B,1);
    if reweighting_threshold>0.0
        fprintf('Starting to calculate weights \n...');
        tic
        %Reweighting in MATLAB:            
        weights = (1./(1+sum(squareform(pdist(Y,'hamm')<=reweighting_threshold))))';   
        
        fprintf('Finished calculating weights \n');
        toc
    end
    
    %Calculating the total weight:
    B_eff=sum(weights);
    fprintf('### N = %d B_with_id_seq = %d B = %d B_eff = %.2f q = %d\n',N,B_with_id_seq,B,B_eff,q);
   	
%Prepare inputs to optimizer.
    %Automatic specification of regularization strength based on B_eff. 
    %B_eff>500 means the standard regularization 0.01 is used, while B_eff<=500 means a higher regularization is chosen.
    
    if B_eff>500
        lambda_J=0.01;
    else
        lambda_J=0.1-(0.1-0.01)*B_eff/500;
    end 
    lambda_h=lambda_J;
    scaled_lambda_h=lambda_h*B_eff;   
    scaled_lambda_J=lambda_J*B_eff/2; %Divide by 2 to keep the size of the coupling regularizaion equivalent to symmetric variant of plmDCA.
    
    
    Y=int32(Y);q=int32(q);
    w=zeros(q+q^2*(N-1),N); %Matrix in which to store parameter estimates (column r will contain estimates from g_r).

%Run optimizer.
    tic
    for r=1:N %Going over all positions 1:N (length of a sequence)
        disp(strcat('Minimizing g_r for node r=',int2str(r)))       
        wr=min_g_r(Y,weights,N,q,scaled_lambda_h,scaled_lambda_J,r,options);
        w(:,r)=wr;
    end
    toc

%Extract the coupling estimates from w.
    %reshape(q + 1) to only isolate J_ij values
    %Explains why 20 values diff. from w to JJ (10 positions * 2 states =
    %20 h_i values)
    h=w(1:q,:); %col=position %row=q value h(q,i)
    JJ=reshape(w(q+1:end,:),q,q,N-1,N);
    Jtemp1=zeros(q,q,N*(N-1)/2);
    Jtemp2=zeros(q,q,N*(N-1)/2);
    l=1;
    for i=1:(N-1)
         for j=(i+1):N
            Jtemp1(:,:,l)=JJ(:,:,j-1,i); %J_ij as estimated from from g_i.
            Jtemp2(:,:,l)=JJ(:,:,i,j)'; %J_ij as estimated from from g_j.
            %Jtemp1 and Jtemp2 both have 1/2 of JJ; i.e Jtemp1 has info for
            %(1,5) while Jtemp2 has info for (5,1); both are at index l
            l=l+1;
        end
    end
     JJ_test=0.5*(Jtemp1+Jtemp2);
   
     
    
%A note on gauges: 
%The parameter estimates coming from g_r satisfy the gauge
%	lambda_J*sum_s Jtemp_ri(s,k) = 0
%	lambda_J*sum_k Jtemp_ri(s,k) = lambda_h*htemp_r(s)	
%	sum_s htemp_r(s) = 0.
%Only the couplings are used in what follows.
    
    
%Shift the coupling estimates into the Ising gauge.
    J1=zeros(q,q,N*(N-1)/2);
    J2=zeros(q,q,N*(N-1)/2);
    
    for l=1:(N*(N-1)/2) %For all possible combinations of i,j (divide by 2 to delete redundances like (i,j) (j,i)
        %Equation (26)
        %mean(Jtemp1(:,:,1),2) takes mean at a constant state k, while mean
        %(Jtemp1(:,:,1) takes mean at a constant state L 
        %Last term takes mean over all possible k & L
        J1(:,:,l)=Jtemp1(:,:,l)-repmat(mean(Jtemp1(:,:,l)),q,1)-repmat(mean(Jtemp1(:,:,l),2),1,q)+mean(mean(Jtemp1(:,:,l)));
        J2(:,:,l)=Jtemp2(:,:,l)-repmat(mean(Jtemp2(:,:,l)),q,1)-repmat(mean(Jtemp2(:,:,l),2),1,q)+mean(mean(Jtemp2(:,:,l)));
    end
    
%Take J_ij as the average of the estimates from g_i and g_j.
%See above eq. (16):
    J=0.5*(J1+J2);

%Calculate frob. norms FN_ij.
    NORMS=zeros(N,N); 
    l=1;
    for i=1:(N-1)
        for j=(i+1):N
            %Equation (27) using  Frobenius norm:
            NORMS(i,j)=norm(J(2:end,2:end,l),'fro');
            NORMS(j,i)=NORMS(i,j);
            l=l+1;
        end
    end                 
    
%Calculate scores CN_ij=FN_ij-(FN_i-)(FN_-j)/(FN_--), where '-'
%denotes average.
    %Equation(28):
    norm_means=mean(NORMS)*N/(N-1);
    norm_means_all=mean(mean(NORMS))*N/(N-1);
    CORRNORMS=NORMS-norm_means'*norm_means/norm_means_all;
    output=[];
    matrix_output=[];
    for i=1:(N-1)
        for j=(i+1):N
            output=[output;[i,j,CORRNORMS(i,j)]];
            matrix_output(i,j)=CORRNORMS(i,j);
            matrix_output(j,i)=matrix_output(i,j);
        end
    end
    
    figure(1);
    hm=heatmap(matrix_output);
    hm.Colormap=parula;
 %   hm.XLabel="position i"
 %   hm.YLabel="position j"
  %  title("Score heatmap of training sequences",'fontsize',20) %on random " + cutoff*100 + "% of contacts in sequence data");
    %save in excel
%     for k=1:size(matrix_output,3)
%         xlswrite('Matrix.xlsx',matrix_output(:,:,k),k)
%     end
% %     
%     excelfile=("Scores_Bottom_" + cutoff*100 + "%");
%     xlswrite(excelfile,output);
%     
%     saveas(hm,excelfile)
    
   
end

function [wr]=min_g_r(Y,weights,N,q,scaled_lambda_h,scaled_lambda_J,r,options)
%Creates function object for (regularized) g_r and minimizes it using minFunc.
    r=int32(r);
    funObj=@(wr)g_r(wr,Y,weights,N,q,scaled_lambda_h,scaled_lambda_J,r);
    wr0=zeros(q+q^2*(N-1),1);
    %minFunc is a more efficient "fminunc" function. Minimizes funObj
    wr=minFunc(funObj,wr0,options);
end

function [fval,grad] = g_r(wr,Y,weights,N,q,lambdah,lambdaJ,r)
%Evaluates (regularized) g_r using the mex-file.
    % wr contains the h and J values being evaluated
    % The first q values in that array are h(r) (since there are q possible
    % values for h_r
    % The rest of the array has all the J(r) values, which explains the
    % reshape lines:
	h_r=reshape(wr(1:q),1,q);
	J_r=reshape(wr(q+1:end),q,q,N-1);
	r=int32(r);
	[fval,grad1,grad2] = g_rC(Y-1,weights,h_r,J_r,[lambdah;lambdaJ],r);
    %g_rC returns two values: fval, which is the value of the function to
    %minimize, and grad1/2, which are the gradients of that function, to
    %increase efficiency of code
	grad = [grad1(:);grad2(:)];
end

function [N,B,q,Y] = return_alignment(inputfile)
%Reads alignment from inputfile, removes inserts and converts into numbers.

    align_full = (inputfile);
    N = size(align_full,2); %Lenght of one sequence (L)
    B = length(align_full); %Number of sequences/rows (M)
    Y = zeros(B,N);

    [Y,L]=letters_to_numbers(align_full);
    q=max(max(Y));
end

%Function "letters to numbers" returns the matrix "CompleteSeq" as the following numeric values;
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
                SeqM(j,i)=1;
            end       
        end
    end   
    
    %Taking the length L (number of columns) of the matrix:
    L=size(SeqM,2); 
    
end




% Copyright 2014 - by Magnus Ekeberg (magnus.ekeberg@gmail.com)
% All rights reserved
% 
% Permission is granted for anyone to copy, use, or modify this
% software for any uncommercial purposes, provided this copyright 
% notice is retained, and note is made of any changes that have 
% been made. This software is distributed without any warranty, 
% express or implied. In no event shall the author or contributors be 
% liable for any damage arising out of the use of this software.
% 
% The publication of research using this software, modified or not, must include 
% appropriate citations to:
%
% 	M. Ekeberg, C. Lövkvist, Y. Lan, M. Weigt, E. Aurell, Improved contact
% 	prediction in proteins: Using pseudolikelihoods to infer Potts models, Phys. Rev. E 87, 012707 (2013)
%
%	M. Ekeberg, T. Hartonen, E. Aurell, Fast pseudolikelihood
%	maximization for direct-coupling analysis of protein structure
%	from many homologous amino-acid sequences, J. Comput. Phys. 276, 341-356 (2014)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








