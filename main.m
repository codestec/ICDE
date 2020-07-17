clc;clear all;
tic;

format long;format compact;
rand('seed',sum(100*clock));

fprintf('improved_ICDE, table2 - table5 ，5e+3 \n');
% the set of all the problems that will be tested
problem_index = 1:24;

%choose the problem
for problem_sel = 16
    
    %the problem that will be tested
    problem = problem_index(problem_sel)
    
    switch problem
        
        case 1
            %variable bounds
            lu=[0 0 0 0 0 0 0 0 0 0 0 0 0; 1 1 1 1 1 1 1 1 1 100 100 100 1];
            %the number of the decision variables
            n=13;
            %the array of parameters p,q,r
            A=[];
            %\ the obejective function value of the best known solution is : -15 \
            
        case 2
            lu=[zeros(1,20); 10*ones(1,20)]; 
            n=20;A=[];%\\\-0.803619\
            
        case 3
            lu=[zeros(1,10); ones(1,10)];
            n=10;A=[];%\\\-1\
            
        case 4
            lu=[78 33 27 27 27; 102 45 45 45 45];
            n=5;A=[];%-30665.539\
            
        case 5
            lu=[0 0 -0.55 -0.55; 1200 1200 0.55 0.55];
            n=4;A=[];%\\5126.4981\
            
        case 6
            lu=[13 0; 100 100];
            n=2;A=[];%\\-6961.81388\
            
        case 7
            lu=[-10*ones(1,10); 10*ones(1,10)];
            n=10;A=[];%\\24.306\
            
        case 8
            lu=[0 0; 10 10];
            n=2;A=[];%\-0.095825\
            
        case 9
            lu=[-10 -10 -10 -10 -10 -10 -10; 10 10 10 10 10 10 10];
            n=7;A=[];%\680.6300573\
            
        case 10
            lu=[100 1000 1000 10 10 10 10 10; 10000 10000 10000 1000 1000 1000 1000 1000];
            n=8;A=[];%\7049.2480\
            
        case 11
            lu=[-1 -1; 1 1];
            n=2;A=[];%\\0.75\
            
        case 12
            lu=[0 0 0; 10 10 10];
            n=3;
            l=1;
            for i=1:9
                for j=1:9
                    for k=1:9
                        A(l,:)=[i j k];
                        l=l+1;
                    end
                end
            end  %\\-1\
            
        case 13
            lu=[-2.3 -2.3 -3.2 -3.2 -3.2; 2.3 2.3 3.2 3.2 3.2];
            n=5;A=[];%\\0.0539498\
            
        case 14
            lu=[zeros(1,10); 10*ones(1,10)];
            n=10;A=[]; %\\-47.7648884595\
            
        case 15
            lu=[zeros(1,3);10*ones(1,3)];
            n=3;A=[]; %\\961.7150222899\
            
        case 16
            lu=[704.4148 68.6 0 193 25;906.3855 288.88 134.75 287.0966 84.1988];
            n=5;A=[];%\\-1.9051552586\
            
        case 17
            lu=[0 0 340 340 -1000 0; 400 1000 420 420 1000 0.5236];
            n=6;A=[];%\\8853.5338748065\
            
        case 18
            lu=[-10 -10 -10 -10 -10 -10 -10 -10 0; 10 10 10 10 10 10 10 10 20];
            n=9;A=[];%\\-0.8660254038\
            
        case 19
            lu=[zeros(1,15); 10*ones(1,15)];
            n=15;A=[];%\\32.6555929502\
            
        case 20
            lu=[zeros(1,24);10*ones(1,24)];
            n=24;A=[];%\\0.2049794002\
            
        case 21
            lu=[0 0 0 100 6.3 5.9 4.5;1000 40 40 300 6.7 6.4 6.25];
            n=7;A=[];%\\193.7245100700\
            
        case 22
            lu=[0 0 0 0 0 0 0 100 100 100.01 100 100 0 0 0  0.01 0.01 -4.7 -4.7 -4.7 -4.7 -4.7;20000 10^6 10^6 10^6 4*10^7 ...
                4*10^7 4*10^7 299.99 399.99 300 400 600 500 500 500 300 400 6.25 6.25 6.25 6.25 6.25];
            n=22;A=[];%\\236.4309755040\
            
        case 23
            lu=[0 0 0 0 0 0 0 0 0.01; 300 300 100 200 100 300 100 200 0.03];
            n=9;A=[];%\\-400.0551\
            
        case 24
            lu=[0 0; 3 4];
            n=2;A=[];%\\\-5.5080132716\
            
    end

    %the size of the parent population
    mu=70;
    %the size of the offspring population
    lambda=210;
    
    %the maximum number of the fitness evaluations(FES)
    MAX_FES=5000;
    
    %the array used to store the objective function values of the best individuals in the final parent population in each run
    outcome_1=[];
    fea_outcome_1 = [];
    infea_outcome_1 = [];
    %the array used to store the feasibility proportion of the final parent population in each run
    outcome_2=[];
    fea_outcome_2 = [];
    infea_outcome_2 = [];
    
    x_best_array = [];
    fea_x_best_arr = [];
    infea_x_best_arr = [];
    num_array = [];
    fea_num_arr = [];
    infea_num_arr = [];
    fea_average_v = [];
    infea_average_v = [];
    average_v = [];
    fun_err_arr = [];
    
    %store the function values of the best solution of the 24 test funcitons
    fbest=[ -15.0000000000, -0.8036191042, -1.0005001000, -30665.5386717834,...
        5126.4967140071, -6961.8138755802, 24.3062090681, -0.0958250415,...
        680.6300573745, 7049.2480205286, 0.7499000000, -1.0000000000,...
        0.0539415140, -47.7648884595, 961.7150222899, -1.9051552586,...
        8853.5338748065, -0.8660254038, 32.6555929502, 0.2049794002,...
        193.7245100700, 236.4309755040, -400.0551000000, -5.5080132716 ];
    
    %the tolarence value for the equality constraints 
    delta=0.0001;
   
    %the total number of generations in each run
    total_gen=ceil(5000/lambda)+1;
    %the value of the threshold parameter
    r=0.6;
    %the threshold generation
    threshold_gen=ceil(r*total_gen);
    
    %set the value of two flags ( 'equality_constraint'& 'inequality_constraint')that indicate whether the problems have equality constraints and inequality constraints
    %equality_constraint=1： the problem has equality constraints
    %equality_constraint=0： the problem has no equality constraints
    %inequality_constraint=1： the problem has inequality constraints
    %inequality_constraint=0： the problem has no inequality constraints
    if  problem==3   ||  problem==5    ||  problem==11 ||  problem==13 ||  problem==14 ||  problem==15|| problem==17 ||  problem==20 ||  problem==21 ||  problem==22 ||  problem==23
        equality_constraint = 1;
        if problem==3   ||  problem==11 ||  problem==13 ||  problem==14 ||  problem==15 || problem==17
            inequality_constraint=0;
        else
            inequality_constraint=1;
        end
    else
        equality_constraint = 0;
        inequality_constraint=1;
    end
    
    %the total number of the independent runs
    total_time=25;
    %initialize the value of the variable 'time'
    time=1;
    
    %do total_gen independent runs
    while time<=total_time
     
        fprintf('time=%d\n',time);

        %generate the initial population
        x=ones(mu,1)*lu(1,:)+rand(mu,n).*(ones(mu,1)*(lu(2,:)-lu(1,:)));
        
        %initialize the value of the variable 'gen'
        gen=1;
        
%         fprintf('gen=%d\n',gen);

        %determine the criterion to compute the degree of the constraint violation
        normalization=1;
        
        %evaluate the initial population
        %and use the first criterion to compute the degree of the constraint violation
        [fit,g,h]=fitness(x,problem,delta,A,normalization);
        
        %if the difference between the maximum value of the constraint violation and the minimum value of the constraint violation is not larger than 200, then normalization=2
        if max(max([g,h]))-min(max([g,h])) <= 200
            normalization=2;
            % reevaluate the initial population using the second criterion to compute the degree of the constraint violation
            fit=recalculate_fitness(fit,problem,delta,normalization,g,h);
        end

        %compute the feasibility proportion of the initial population
        percent=length(find(fit(:,2)==0))/mu;
        
        %modify the value of FES
        FES=mu;

        while FES<=MAX_FES

            %generate an offspring population from the parent population by DE
            temp_x=select_reproduce(x,mu,lu,n,gen,threshold_gen,total_gen);
            %evaluate the offspring population
            [temp_fit,temp_g,temp_h]=fitness(temp_x,problem,delta,A,normalization);
            
            %modify the value of FES
            FES=FES+lambda;
            
            %combine the parent population and the offspring population
            x(((mu+1):(mu+lambda)),:)=temp_x;
            fit(((mu+1):(mu+lambda)),:)=temp_fit;
            if normalization==1
                if inequality_constraint==1
                    g(((mu+1):(mu+lambda)),:)=temp_g;
                end
                if equality_constraint==1
                    h(((mu+1):(mu+lambda)),:)=temp_h;
                end
            end
            
            %if the number of the generation is larger then 1 and the
            %combined population is an infeasible population, then do as follows
            if gen>1 && length(find(fit(:,2)>0))==size(fit,1)
                
                %total_num is the total number of the individuals in the population 'non_elite_x'
                total_num=size(non_elite_x,1);
                %generate a number randomly between 1 and total_num
                rand_num=randint(1,1,[1,total_num]);
                %initialize 'ind_index'
                ind_index=zeros(1,rand_num);
                %initialize 'temp_index'
                temp_index=(1:total_num);

                %generate 'rand_num' different numbers between 1 and
                %'total_num' randomly and put them into the set 'ind_index'
                for i=1:rand_num
                    position=floor(total_num*rand)+1;
                    ind_index(i)=temp_index(position);
                    temp_index(position)=[];
                    total_num=total_num-1;
                end

                %put these individuals selected randomly into the combined population
                x(((lambda+mu+1):(lambda+mu+rand_num)),:)=non_elite_x(ind_index,:);
                fit(((lambda+mu+1):(lambda+mu+rand_num)),:)=non_elite_fit(ind_index,:);
                if normalization==1
                    if inequality_constraint==1
                        g(((lambda+mu+1):(lambda+mu+rand_num)),:)=non_elite_g(ind_index,:);
                    end
                    if equality_constraint==1
                        h(((lambda+mu+1):(lambda+mu+rand_num)),:)=non_elite_h(ind_index,:);
                    end
                end

            end
            
            %if normalization=1, then reevaluate the combined population
            if normalization == 1
                fit=recalculate_fitness(fit,problem,delta,normalization,g,h);
            end

            %compute the feasibility proportion of the final combined population
            total_ind_num=size(fit,1);
            percent=length(find(fit(:,2)==0))/total_ind_num;

            %the combined population has only infeasible individuals
            if percent == 0

                %select mu individuals from the combined population by the ATM strategy
                nouse=replacement(fit,total_ind_num,mu);

                %store all the remaining individuals into the population 'non_elite_x'
                non_elite_x=x;
                non_elite_x(nouse,:)=[];
                non_elite_fit=fit;
                non_elite_fit(nouse,:)=[];
                if normalization==1
                    if inequality_constraint==1
                        non_elite_g=g;
                        non_elite_g(nouse,:)=[];
                    end
                    if equality_constraint==1
                        non_elite_h=h;
                        non_elite_h(nouse,:)=[];
                    end
                end

                %obtain the new parent population
                x=x(nouse,:);
                fit=fit(nouse,:);
                if normalization==1
                    if inequality_constraint==1
                        g=g(nouse,:);
                    end
                    if equality_constraint==1
                        h=h(nouse,:);
                    end
                end

            %the combined population has both infeasible individuals and feasible individuals
            elseif percent>0 && percent<1

                %if normalization==2, normalize the degree of the constraint violation
                if normalization==2
                    
                    %find the feasible individuals and infeasible individuals
                    temp_1=find(fit(:,2)==0);
                    temp_2=find(fit(:,2)>0);
                    
                    %seperate the value of the objective function and the degree of the constraint violation
                    fit_1=fit(:,1);
                    fit_2=fit(:,2);
                    
                    %find the minmum value and the maxmum value of the objective function of the feasible individuals
                    min_fit=min(fit_1(temp_1));
                    max_fit=max(fit_1(temp_1));
                    
                    %find the minmum value and the maxmum value of the degree 
                    %of the constraint violation of the infeasible individuals
                    min_vio=min(fit_2(temp_2));
                    max_vio=max(fit_2(temp_2));
                    
                    %normalize the value of the objective function of all the individuals
                    fit_1(temp_2)=max(fit_1(temp_2),ones(size(temp_2,1),1)*(percent*min_fit+(1-percent)*max_fit));
                    fit_1=(fit_1-min(fit_1))./((max(fit_1)-min(fit_1))+1E-30);
                    
                    %normalize the value of the degree of the constraint violation of all the individuals
                    fit_2(temp_2)=(fit_2(temp_2)-min_vio)./(max_vio-min_vio+1E-30);
                    
                    %compute the final value of the fitness function of all the individuals
                    fit_v=fit_1+fit_2;
                 
                %if normalization==1, do not normalize the degree of the constraint violation
                else

                    %find the feasible individuals and infeasible individuals
                    temp_1=find(fit(:,2)==0);
                    temp_2=find(fit(:,2)>0);
                    
                    %seperate the value of the objective function and the degree of the constraint violation
                    fit_1=fit(:,1);
                    fit_2=fit(:,2);
                    
                    %find the minmum value and the maxmum value of the objective function of the feasible individuals
                    min_fit=min(fit_1(temp_1));
                    max_fit=max(fit_1(temp_1));
                    
                    %normalize the value of the objective function of all the individuals
                    fit_1(temp_2)=max(fit_1(temp_2),ones(size(temp_2,1),1)*(percent*min_fit+(1-percent)*max_fit));
                    fit_1=(fit_1-min(fit_1))./((max(fit_1)-min(fit_1))+1E-30);
                    
                    %compute the final value of the fitness function of all the individuals
                    fit_v=fit_1+fit_2;
                    
                end

                %sort all the individuals in the non-descending order according to fit_v
                [fitVal,fitIndex]=sort(fit_v);

                %select the first mu individuals to construct the new parent population
                x=x(fitIndex(1:mu),:);
                fit=fit(fitIndex(1:mu),:);
                if normalization==1
                    if inequality_constraint==1
                        g=g(fitIndex(1:mu),:);
                    end
                    if equality_constraint==1
                        h=h(fitIndex(1:mu),:);
                    end
                end
                
            %the combined population has only feasible individuals
            else

                %sort all the individuals in the non-descending order
                %according to the value of their objective function
                [fitVal,fitIndex]=sort(fit(:,1));

                %select the first mu individuals to construct the new parent population
                x=x(fitIndex(1:mu),:);
                fit=fit(fitIndex(1:mu),:);
                if normalization==1
                    if inequality_constraint==1
                        g=g(fitIndex(1:mu),:);
                    end
                    if equality_constraint==1
                        h=h(fitIndex(1:mu),:);
                    end
                end
            end

            %modify the value of 'gen'
            gen=gen+1;
            
%             fprintf('gen=%d\n',gen);
            
        end
        
        %find the feasible individuals in the final population
        num_3=find(fit(:,2)==0);
        %the final population has feasible individuals

        if ~isempty(num_3)
            
            fea_x = x( num_3 , : );
            fea_fit = fit( num_3 , 1 );
            
            % 找出最小的目标函数值
            [ min_value, min_nouse ] = min( fea_fit );
            
            % 最优个体
            x_best = x( min_nouse , : );            
            
            fea_x_best_arr = [ fea_x_best_arr ; x_best ];
            
            fea_num_arr = [ fea_num_arr , 0 ];
            
            fea_average_v = [ fea_average_v , 0 ];
            
            outcome_1 = [ outcome_1, min_value ];
            fea_outcome_1 = [ fea_outcome_1 , min_value ];
            
            percent = length(num_3) / size(fit,1);
            outcome_2 = [outcome_2 percent];
            fea_outcome_2 = [ fea_outcome_2 , percent ];
            
%             %save the value of the objective function of the best individuals in the final population
%             outcome_1 = [outcome_1, min(fit(num_3,1))];
%             %save the feasibility proportion of the final population
%             percent = length(num_3) / size(fit,1);
%             outcome_2=[outcome_2 percent];
         
        %the final population has no feasible individuals
        else

            [ m_ , n_ ] = min( fit(:,2) );
            
            x_best = x( n_ , : );
            
            infea_x_best_arr = [ infea_x_best_arr ; x_best ];
            
            [ f_1, g_1, h_1 ] = fitness( x_best, problem, delta, A, normalization );
            
            num_1 = 0;
            if ~isempty(g_1)
                num_1 = num_1 + length( find( g_1>0 ) );
            end
            if ~isempty(h_1)
                num_1 = num_1 + length( find( abs(h_1)>0.0001 ) );
            end
            infea_num_arr = [ infea_num_arr , num_1 ];
            
            new_g = max( 0 , g_1 );
            new_h = zeros(1,size(h_1,2));
            index_h = find( abs(h_1) > 0.0001 );
            new_h(index_h) = abs( h_1(index_h) );
            violations = [ new_g , new_h ];
            average_vio = sum( violations ) / ( size(g_1,2) + size(h_1,2) );
            infea_average_v = [ infea_average_v , average_vio ];            
            
            outcome_1 = [ outcome_1 , fit(n_,1) ];
            infea_outcome_1 = [ infea_outcome_1 , outcome_1 ];
            
            outcome_2 = [ outcome_2 , 0 ];
            infea_outcome_2 = [ infea_outcome_2 , 0 ];
            
        end
        
        %modify the value of 'time'
        time = time + 1;
        
    end
    
    %all the feasibility proportion of the final population
    all_outcome_2 = outcome_2
    %all the value of the objective function of the best individuals in the final population
    all_outcome_1 = outcome_1
    %the maximum value
    max_outcome_1 = max(outcome_1)
    %the minimum value
    min_outcome_1 = min(outcome_1)
    %the median value
    median_outcome_1 = median(outcome_1)
    %the mean value
    mean_outcome_1 = mean(outcome_1)
    %the standard value
    std_outcome_1 = std(outcome_1)
   
    %***************************************************************
    fprintf('\nthe function value of the theoretical best solution: ');
    f_xbest = fbest(problem)
    fprintf('\nthe function error values: ');
    fun_errs = outcome_1 - f_xbest;
    fea_fun_errs = fea_outcome_1 - f_xbest;
    infea_fun_errs = infea_outcome_1 - f_xbest;
    
    total_n = 0;
    
    if ~isempty(fea_x_best_arr)
        [ val , ind ] = sort( fea_fun_errs );
        total_n = total_n + size( fea_x_best_arr , 1 );
        
        x_best_array = [ x_best_array ;  fea_x_best_arr(ind,:) ];
        
        average_v = [ average_v , fea_average_v(ind) ];
        
        fun_err_arr = [ fun_err_arr , fea_fun_errs(ind) ];
        
        num_array = [ num_array , fea_num_arr(ind) ];
        
    end 
    if ~isempty(infea_x_best_arr)
        [ val_in , ind_in ] = sort( infea_average_v );
        total_n = total_n + size( infea_x_best_arr , 1 );
        x_best_array = [ x_best_array ; infea_x_best_arr( ind_in , : ) ];
        average_v = [ average_v , infea_average_v(ind_in) ];
        fun_err_arr = [ fun_err_arr , infea_fun_errs(ind_in) ];
        num_array = [ num_array , infea_num_arr(ind_in) ];
    end
        
    %find the best value of the fucntion error values
    %sort all the values in fun_errs
    
    fprintf('\nthe best solution: ');
    Best_x = x_best_array( 1 , : )
    fprintf('\nthe best value of the function error values: ');    
    Best_fun_err = fun_err_arr( 1 )
    fprintf('\nthe number of constraints violated by the best solution: ');
    Best_num_vio_cons = num_array( 1 )
    
    fprintf('\nthe median solution: ');
    Median_x = x_best_array( 13 , : )
    fprintf('\nthe median value of the function error: ');
    Median_fun_err = fun_err_arr( 13 )
    fprintf('\nthe number of constraints violated by the median solution: ');
    Median_num_vio_cons = num_array( 13 )
    
    fprintf('\nthe worst solution: ');
    Worst_x = x_best_array( 25 , : )
    fprintf('\nthe worst value of the function error values: ');    
    Worst_fun_err = fun_err_arr( 25 )
    fprintf('\nthe number of constraints violated by the worst solution: ');
    Worst_num_vio_cons = num_array( 25 )

    fprintf('\nthe mean value of the function error: ');    
    Mean_fun_err = mean( fun_err_arr )
    
    fprintf('\nthe standard deviation of the function error: ');
    Std_deviation = std( fun_err_arr )
    
    fprintf('\nthe standard error of the function error: ');
    % Std_error = Std_deviation / sqrt(25) , the sample size = 25
    Std_error = Std_deviation / 5
    
    fprintf('\nthe number of violated constraints at the median solution c: ');
    if num_array(13)==0
        fprintf('\nthe number of violated constraints at the median solution c:(0,0,0) ');
    elseif num_array(13)<0
        fprintf('\n when computing c,num_array(13)<0, error');
        return
    elseif num_array(13)>0
        [ f_2, g_2, h_2 ] = fitness( Median_x, problem, delta, A, normalization );
        num_large = 0;
        num_between = 0;
        num_small = 0;
        if ~isempty(g_2)
            num_large = num_large + length( find( g_2>1.0 ) );
            num_between = num_between + length( find( 0.01<g_2 & g_2<=1.0 ) );
            num_small = num_small + length( find( 0.0001<g_2 & g_2<=0.01 ) );
        end
        if ~isempty(h_2)
            num_large = num_large + length( find( abs(h_2)>1.0 ) );
            num_between = num_between + length( find( 0.01<abs(h_2) & abs(h_2)<=1.0 ) );
            num_small = num_small + length( find( 0.0001<abs(h_2) & abs(h_2)<=0.01 ) );
        end
        fprintf('\nthe number of constraints c: ( %d , %d , %d ).',num_large,num_between,num_small);
    end
    
    %compuate the average violations 'v'
    fprintf('\nthe mean violations of all constrains: v =  %d .',average_v(13));
    
    fprintf('\nall the 25 solutions: ');
    x_best_array
    
    fprintf('\n all the numbers of constraints violated by solutions: ');
    num_array
    
    fprintf('\nall the function error values of 25 solutions:\n');
    fun_err_arr
    
end

%the time of the run
run_time = toc
