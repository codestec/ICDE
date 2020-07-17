% the combined population is in the infeasible situation, the function will
% select mu individuals from the population

function number=replacement(fit,total_elem_num,mu)

%% input parameters
    % fit -- the values of the fitness functions
    % total_elem_num -- the total element number
    % mu -- the number of the selected individuals
%% output parameters
    % number -- the index of the selected individuals

%the index of all the individuals in the population
temp_number=1:total_elem_num;

%popsize--the number of the individuals that has been selected
popsize=0;

%number -- the set of the indexs
number=[];

%select mu individuals from the population
while popsize<mu
    
    % R -- the set of the rank values of all the individuals
    R=[];
    
    %compute the rank values of individuals
    for i=1:length(temp_number)
        R(i)=length(find((fit(i,1)>=fit(:,1)) & (fit(i,2)>=fit(:,2)) & ((fit(i,1)>fit(:,1)) | (fit(i,2)>fit(:,2)))));
    end
    
    % find individuals whoes rank value is equal to 0, i.e., the Pareto front
    temp=(find(R==0));
    
    %sort individuals in the Pareto front by the degree of constraint violations
    [mouse,nouse]=sort(fit(temp,2));
    temp=temp(nouse);
    
    %select the half of individuals in the Pareto front after sorting
    if round(length(temp)*1/2)>=1
        temp=temp(1:round(length(temp)*1/2));
    end
    
    %put the indexs of selected individuals into the set 'number'
    number=[number temp_number(temp)];
    %delete the indexs of selected individuals from 'temp_number'
    temp_number(temp)=[];
    %delete the fitness values of selected individuals from 'fit'
    fit(temp,:)=[];
    %recalculate the number of selected individuals
    popsize=popsize+length(temp);
    
end
%select the first mu individuals
number=number(1:mu);
