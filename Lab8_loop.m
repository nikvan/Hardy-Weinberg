% Clear all variable, Close all figures, Clear command window
clear all
close all
clc

% Define some constants that we will use to populate our gene pool
AA='AA';
aa='aa';
p_AA = 0.50;%default 0.75
p_aa = 0.50;%default 0.25


N=1000; %originally 1000 % Population Size

% Preallocate Space to hold our Gene Pool.
% Repmat creates N copies of '  ' in an N by 1 vector.
Pool=repmat('  ',N,1); 

% Loop N times. 
for i = 1:N
    % If (i/N)>p_aa is true then we add AA to our pool if not add aa
    % This means that the first 25% of then entries will be aa, then
    % the rest will be AA.
    if((i/N)>p_aa)
        Pool(i,:) = AA; % Since AA is a 2 length array of characters we need to use (i,:) notation.
    else
        Pool(i,:) = aa;
    end
end

% Verify that we have created the initial gene pool properly
% Pool(:,1) == 'A' returns 1 if 'A' is in the first position.
% Pool(:,2) == 'A' returns 1 if 'A' is in the first position.
% & returns 1 if both sides are 1 and 0 otherwise
% 1 & 1 -> 1
% 0 & 0 -> 0
% 0 & 1 -> 0
% 1 & 0 -> 0
% sum(v) adds all the numbers in the vector v. 
% Dividing by N gives us the percentage of each in the initial gene pool
start_AA = sum(Pool(:,1) == 'A' & Pool(:,2) == 'A')/N
start_aa = sum(Pool(:,1) == 'a' & Pool(:,2) == 'a')/N

% Use the square law to predict frequency of genotypes in future
% generations
next_AA_predicted=p_AA^2
next_Aa_predicted=2*p_AA*p_aa
next_aa_predicted=p_aa^2

%% Begin Simulation
G=100;  % originally 100 % Number of generations to simulate
% Use repmat to create an N by 1 by G sized matrix to hold the genes for
% each future generation. 
Generations=repmat('  ',N,1,G);
 
% We will use a for-loop where each loop iteration represents 1 generation
for gen = 1:G
    
    % Use repmat to create an N by 1 sized matrix to hold the genes for
    % a single generation.
    OffSpring=repmat('  ',N,1);
    
    % In this for-loop each iteration represents an individual in the next
    % generation.
    for i = 1:N
        mate1 = randi(N); % Randomly select first parent
        mate2 = randi(N); % Randomly select second parent      
        flip1 = round(rand); % Randomly select allele provided by first parent (coin flip)
        flip2 = round(rand); % Randomly select allele provided by first parent (coin flip)
        if(flip1)
            Allele1 = Pool(mate1,1); % If "heads" use first allele from first parent
        else
            Allele1 = Pool(mate1,2); % If "tails" use second allele from first parent
        end
        if(flip2)
            Allele2 = Pool(mate2,1); % If "heads" use first allele from second parent
        else
            Allele2 = Pool(mate2,2); % If "tails" use second allele from second parent
        end
        OffSpring(i,:) = strcat(Allele1,Allele2); % combine two alleles to create offspring.
    end
    Generations(:,:,gen) = OffSpring; % Add vector of OffSpring to Generations.
    
    % Calculate the number of each genotype in OffSpring
    % sum(v) adds all the numbers in the vector v. 
    % Dividing by N gives us the percentage of each in the offspring
    next_AA(gen) = sum(OffSpring(:,1) == 'A' & OffSpring(:,2) == 'A')/N
    next_Aa_exact = sum(OffSpring(:,1) == 'A' & OffSpring(:,2) == 'a')/N;
    next_aA_exact = sum(OffSpring(:,1) == 'a' & OffSpring(:,2) == 'A')/N;
    next_Aa(gen) = next_Aa_exact + next_aA_exact;
    next_aa(gen) = sum(OffSpring(:,1) == 'a' & OffSpring(:,2) == 'a')/N
end



%% Generate Figures
figure(1)
plot(next_AA,'linewidth',2)
hold on;
plot(next_AA_predicted*ones(size(next_AA)),'linewidth',2)
hold off
title('Frequency of AA Genotype')
xlabel('Generation')
ylabel('Frequency')
legend('Simulated Frequency', 'Predicted Frequency')

figure(2)
plot(next_Aa,'linewidth',2)
hold on;
plot(next_Aa_predicted*ones(size(next_Aa)),'linewidth',2)
hold off
title('Frequency of Aa Genotype')
xlabel('Generation')
ylabel('Frequency')
legend('Simulated Frequency', 'Predicted Frequency')

figure(3)
plot(next_aa,'linewidth',2)
hold on;
plot(next_aa_predicted*ones(size(next_aa)),'linewidth',2)
hold off
title('Frequency of aa Genotype')
xlabel('Generation')
ylabel('Frequency')
legend('Simulated Frequency', 'Predicted Frequency')

% figure(4)
% plot(next_aa./N) %divide by population N to get percentage?
% hold on
% plot(next_Aa./N)
% plot(next_AA./N)
% plot(next_aa_predicted*ones(size(next_aa)),'linewidth',2)
% plot(next_Aa_predicted*ones(size(next_Aa)),'linewidth',2)
% plot(next_AA_predicted*ones(size(next_AA)),'linewidth',2)
% legend('aa Calculated','Aa Calculated','AA Calculated','aa Predicted','Aa Predicted','AA Predicted');
