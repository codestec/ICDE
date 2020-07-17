
% recalculate the degree of constraint violations of individuals

function old_fit=recalculate_fitness(old_fit,problem,delta,normalization,g,h)

%% input parameters
    % old_fit -- the old values
    % problem -- the problem tested
    % delta -- the tolerance value
    % normalization -- the value of the flag
    % g -- degree of the inequality constraint violations
    % h -- degree of the equality constraint violations
%% output parameters
    % old_fit -- the new calculated values
    
%compute the size of the population p
popsize=size(old_fit,1);

%choose the problem
switch problem

    case 1

        %compute the degree of constraint violations according to
        %'normalization'
        if normalization==1

            old_fit(:,2)=sum(max(0,g)./(ones(popsize,1)*(max(max(0,g))+1E-30)),2)./size(g,2);

        else

            old_fit(:,2)=sum(max(0,g),2);

        end

    case 2
        
        %compute the degree of constraint violations according to
        %'normalization'
        if normalization==1

            old_fit(:,2)=sum(max(0,g)./(ones(popsize,1)*(max(max(0,g))+1E-30)),2)./size(g,2);

        else

            old_fit(:,2)=sum(max(0,g),2);

        end

    case 3
        
        %compute the degree of constraint violations according to
        %'normalization'
        if normalization==1

            old_fit(:,2)=sum(max(0,abs(h)-delta)./(ones(popsize,1)*(max(max(0,abs(h)-delta))+1E-30)),2)./(size(h,2));

        else

            old_fit(:,2)=sum(max(0,abs(h)-delta),2);

        end

    case 4

        %compute the degree of constraint violations according to
        %'normalization'
        if normalization==1

            old_fit(:,2)=sum(max(0,g)./(ones(popsize,1)*(max(max(0,g))+1E-30)),2)./size(g,2);

        else

            old_fit(:,2)=sum(max(0,g),2);

        end

    case 5
        
        %compute the degree of constraint violations according to
        %'normalization'
        if normalization==1

            old_fit(:,2)=(sum(max(0,g)./(ones(popsize,1)*(max(max(0,g))+1E-30)),2)+sum(max(0,abs(h)-delta)./(ones(popsize,1)*(max(max(0,abs(h)-delta))+1E-30)),2))./(size(g,2)+size(h,2));

        else

            old_fit(:,2)=sum(max(0,g),2)+sum(max(0,abs(h)-delta),2);

        end

    case 6
        
        %compute the degree of constraint violations according to
        %'normalization'
        if normalization==1

            old_fit(:,2)=sum(max(0,g)./(ones(popsize,1)*(max(max(0,g))+1E-30)),2)./size(g,2);

        else

            old_fit(:,2)=sum(max(0,g),2);

        end

    case 7
        
        %compute the degree of constraint violations according to
        %'normalization'
        if normalization==1

            old_fit(:,2)=sum(max(0,g)./(ones(popsize,1)*(max(max(0,g))+1E-30)),2)./size(g,2);

        else

            old_fit(:,2)=sum(max(0,g),2);

        end

    case 8
        
        %compute the degree of constraint violations according to
        %'normalization'
        if normalization==1

            old_fit(:,2)=sum(max(0,g)./(ones(popsize,1)*(max(max(0,g))+1E-30)),2)./size(g,2);

        else

            old_fit(:,2)=sum(max(0,g),2);

        end

    case 9
        
        %compute the degree of constraint violations according to
        %'normalization'
        if normalization==1

            old_fit(:,2)=sum(max(0,g)./(ones(popsize,1)*(max(max(0,g))+1E-30)),2)./size(g,2);

        else

            old_fit(:,2)=sum(max(0,g),2);

        end

    case 10
        
        %compute the degree of constraint violations according to
        %'normalization'
        if normalization==1

            old_fit(:,2)=sum(max(0,g)./(ones(popsize,1)*(max(max(0,g))+1E-30)),2)./size(g,2);

        else

            old_fit(:,2)=sum(max(0,g),2);

        end

    case 11
        
        %compute the degree of constraint violations according to
        %'normalization'
        if normalization==1

            old_fit(:,2)=sum(max(0,abs(h)-delta)./(ones(popsize,1)*(max(max(0,abs(h)-delta))+1E-30)),2)./(size(h,2));

        else

            old_fit(:,2)=sum(max(0,abs(h)-delta),2);

        end

    case 12
        
        %compute the degree of constraint violations according to
        %'normalization'
        if normalization==1

            old_fit(:,2)=sum(max(0,g)./(ones(popsize,1)*(max(max(0,g))+1E-30)),2)./size(g,2);

        else

            old_fit(:,2)=sum(max(0,g),2);

        end

    case 13
        
        %compute the degree of constraint violations according to
        %'normalization'
        if normalization==1

            old_fit(:,2)=sum(max(0,abs(h)-delta)./(ones(popsize,1)*(max(max(0,abs(h)-delta))+1E-30)),2)./(size(h,2));

        else

            old_fit(:,2)=sum(max(0,abs(h)-delta),2);

        end

    case 14
        
        %compute the degree of constraint violations according to
        %'normalization'
        if normalization==1

            old_fit(:,2)=sum(max(0,abs(h)-delta)./(ones(popsize,1)*(max(max(0,abs(h)-delta))+1E-30)),2)./(size(h,2));

        else

            old_fit(:,2)=sum(max(0,abs(h)-delta),2);

        end

    case 15
        
        %compute the degree of constraint violations according to
        %'normalization'
        if normalization==1

            old_fit(:,2)=sum(max(0,abs(h)-delta)./(ones(popsize,1)*(max(max(0,abs(h)-delta))+1E-30)),2)./(size(h,2));

        else

            old_fit(:,2)=sum(max(0,abs(h)-delta),2);

        end

    case 16

        %compute the degree of constraint violations according to
        %'normalization'
        if normalization==1

            old_fit(:,2)=sum(max(0,g)./(ones(popsize,1)*(max(max(0,g))+1E-30)),2)./size(g,2);

        else

            old_fit(:,2)=sum(max(0,g),2);

        end

    case 17
        
        %compute the degree of constraint violations according to
        %'normalization'
        if normalization==1

            old_fit(:,2)=sum(max(0,abs(h)-delta)./(ones(popsize,1)*(max(max(0,abs(h)-delta))+1E-30)),2)./(size(h,2));

        else

            old_fit(:,2)=sum(max(0,abs(h)-delta),2);

        end

    case 18

        %compute the degree of constraint violations according to
        %'normalization'
        if normalization==1

            old_fit(:,2)=sum(max(0,g)./(ones(popsize,1)*(max(max(0,g))+1E-30)),2)./size(g,2);

        else

            old_fit(:,2)=sum(max(0,g),2);

        end

    case 19
        
        %compute the degree of constraint violations according to
        %'normalization'
        if normalization==1

            old_fit(:,2)=sum(max(0,g)./(ones(popsize,1)*(max(max(0,g))+1E-30)),2)./size(g,2);

        else

            old_fit(:,2)=sum(max(0,g),2);

        end

    case 20
        
        %compute the degree of constraint violations according to
        %'normalization'
        if normalization==1

            old_fit(:,2)=(sum(max(0,g)./(ones(popsize,1)*(max(max(0,g))+1E-30)),2)+sum(max(0,abs(h)-delta)./(ones(popsize,1)*(max(max(0,abs(h)-delta))+1E-30)),2))./(size(g,2)+size(h,2));

        else

            old_fit(:,2)=sum(max(0,g),2)+sum(max(0,abs(h)-delta),2);

        end

    case 21
        
        %compute the degree of constraint violations according to
        %'normalization'
        if normalization==1

            old_fit(:,2)=(sum(max(0,g)./(ones(popsize,1)*(max(max(0,g))+1E-30)),2)+sum(max(0,abs(h)-delta)./(ones(popsize,1)*(max(max(0,abs(h)-delta))+1E-30)),2))./(size(g,2)+size(h,2));

        else

            old_fit(:,2)=sum(max(0,g),2)+sum(max(0,abs(h)-delta),2);

        end

    case 22
        
        %compute the degree of constraint violations according to
        %'normalization'
        if normalization==1

            old_fit(:,2)=(sum(max(0,g)./(ones(popsize,1)*(max(max(0,g))+1E-30)),2)+sum(max(0,abs(h)-delta)./(ones(popsize,1)*(max(max(0,abs(h)-delta))+1E-30)),2))./(size(g,2)+size(h,2));

        else

            old_fit(:,2)=sum(max(0,g),2)+sum(max(0,abs(h)-delta),2);

        end

    case 23
        
        %compute the degree of constraint violations according to
        %'normalization'
        if normalization==1
            old_fit(:,2)=(sum(max(0,g)./(ones(popsize,1)*(max(max(0,g))+1E-30)),2)+sum(max(0,abs(h)-delta)./(ones(popsize,1)*(max(max(0,abs(h)-delta))+1E-30)),2))./(size(g,2)+size(h,2));
        else
            old_fit(:,2)=sum(max(0,g),2)+sum(max(0,abs(h)-delta),2);
        end
        
    case 24
        
        %compute the degree of constraint violations according to
        %'normalization'
        if normalization==1
            old_fit(:,2)=sum(max(0,g)./(ones(popsize,1)*(max(max(0,g))+1E-30)),2)./size(g,2);
        else
            old_fit(:,2)=sum(max(0,g),2);
        end
        
end
