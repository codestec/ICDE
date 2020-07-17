% evolve the parent population by (mu+lambda)-DE and generate the offspring population 

function temp_x=select_reproduce(x,mu,lu,n,gen,threshold_gen,total_gen)
%% input parameters
    % x -- the parent population
    % mu -- the size of the parent population
    % lu -- the range of all the parameters
    % n -- the dimention of the variables
    % gen -- the current generation
    % threshold_gen -- the threshold generation
    % total_gen -- the total generation
%% output parameters
    % temp_x -- the generated offspring population
    
% set the values of these parameters F,CR,pm
F=0.8;
CR=0.9;
pm=0.05;

%initialize the set 'temp_x'
temp_x=zeros(mu*3,n);

%evolve all the parent individuals by (mu+lambda)-DE
for i=1:mu

   %% "rand/1" strategy ....%
    
    %select two different individuals which are different from the current
    %individual from the parent population randomly
    nouse_=1:mu;
    nouse_(i)=[];

    temp=floor(rand*(mu-1))+1;
    nouse(1)=nouse_(temp);
    nouse_(temp)=[];

    temp=floor(rand*(mu-2))+1;
    nouse(2)=nouse_(temp);
    nouse_(temp)=[];

    temp=floor(rand*(mu-3))+1;
    nouse(3)=nouse_(temp);

    %do the mutation operation 'rand/1'
    v1=x(nouse(1),:)+F.*(x(nouse(2),:)-x(nouse(3),:));

    %modify the values of some elements violating the boudary constraint,
    %by reflecting them back from the violated boundary
    w=find(v1<lu(1,:));
    if ~isempty(w)
        v1(1,w)=2.*lu(1,w)-v1(1,w);
        w1=find(v1(1,w)>lu(2,w));
        if ~isempty(w1)
            v1(1,w(w1))=lu(2,w(w1));
        end
    end
    y=find(v1>lu(2,:));
    if ~isempty(y)
        v1(1,y)=2.*lu(2,y)-v1(1,y);
        y1=find(v1(1,y)<lu(1,y));
        if ~isempty(y1)
            v1(1,y(y1))=lu(1,y(y1));
        end
    end

    %do the binomial crossover operation
    j_rand=floor(rand*n)+1;
    t=rand(1,n)<CR;
    t(1,j_rand)=1;
    t_=1-t;
    %obtain the first offspring individual
    u1=t.*v1+t_.*x(i,:);

    %% "current to rand/1" & 'current to best/1 'strategy
    
    %select three different individuals which are different from the current
    %individual from the parent population randomly
    nouse_=1:mu;
    nouse_(i)=[];

    temp=floor(rand*(mu-1))+1;
    nouse(1)=nouse_(temp);
    nouse_(temp)=[];

    temp=floor(rand*(mu-2))+1;
    nouse(2)=nouse_(temp);
    nouse_(temp)=[];

    temp=floor(rand*(mu-3))+1;
    nouse(3)=nouse_(temp);
    
    %do the mutation operation according to the relationship between gen
    %and threshold_gen
    if gen <= threshold_gen
        v3=x(i,:)+rand*(x(nouse(3),:)-x(i,:))+F.*(x(nouse(1),:)-x(nouse(2),:));
    else
%         v3=x(i,:)+rand*(x(1,:)-x(i,:))+F.*(x(nouse(1),:)-x(nouse(2),:));
        v3=x(i,:)+F*(x(1,:)-x(i,:))+F.*(x(nouse(1),:)-x(nouse(2),:));
    end

    %modify the values of some elements violating the boudary constraint,
    %by reflecting them back from the violated boundary
    w=find(v3<lu(1,:));
    if ~isempty(w)
        v3(1,w)=2.*lu(1,w)-v3(1,w);
        w1=find(v3(1,w)>lu(2,w));
        if ~isempty(w1)
            v3(1,w(w1))=lu(2,w(w1));
        end
    end
    y=find(v3>lu(2,:));
    if ~isempty(y)
        v3(1,y)=2.*lu(2,y)-v3(1,y);
        y1=find(v3(1,y)<lu(1,y));
        if ~isempty(y1)
            v3(1,y(y1))=lu(1,y(y1));
        end
    end

    %do the IBGA operation,and obtain the second offspring individual
    u3=v3;
    if rand<pm
        temp=floor(rand*n)+1;
        mutrange=(lu(2,temp)-lu(1,temp))*(1-gen/total_gen)^6;
        if rand<0.5
            sign_=1;
        else
            sign_=-1;
        end
        u3(1,temp)=u3(1,temp)+mutrange*sign_*2^(-round(rand*15));
        if u3(1,temp)>lu(2,temp)
            u3(1,temp)=lu(2,temp);
        elseif u3(1,temp)<lu(1,temp)
            u3(1,temp)=lu(1,temp);
        end
    end
        
    %% "rand/2" strategy%
    
    %select five different individuals which are different from the current
    %individual from the parent population randomly
    nouse_=1:mu;
    nouse_(i)=[];

    temp=floor(rand*(mu-1))+1;
    nouse(1)=nouse_(temp);
    nouse_(temp)=[];

    temp=floor(rand*(mu-2))+1;
    nouse(2)=nouse_(temp);
    nouse_(temp)=[];

    temp=floor(rand*(mu-3))+1;
    nouse(3)=nouse_(temp);
    nouse_(temp)=[];

    temp=floor(rand*(mu-4))+1;
    nouse(4)=nouse_(temp);
    nouse_(temp)=[];

    temp=floor(rand*(mu-5))+1;
    nouse(5)=nouse_(temp);

    %do the mutation operation 'rand/1'
    v5=x(nouse(1),:)+F.*(x(nouse(2),:)-x(nouse(3),:))+F.*(x(nouse(4),:)-x(nouse(5),:));

    %modify the values of some elements violating the boudary constraint,
    %by reflecting them back from the violated boundary
    w=find(v5<lu(1,:));
    if ~isempty(w)
        v5(1,w)=2.*lu(1,w)-v5(1,w);
        w1=find(v5(1,w)>lu(2,w));
        if ~isempty(w1)
            v5(1,w(w1))=lu(2,w(w1));
        end
    end
    y=find(v5>lu(2,:));
    if ~isempty(y)
        v5(1,y)=2.*lu(2,y)-v5(1,y);
        y1=find(v5(1,y)<lu(1,y));
        if ~isempty(y1)
            v5(1,y(y1))=lu(1,y(y1));
        end
    end

    j_rand=floor(rand*n)+1;
    t=rand(1,n)<CR;
    t(1,j_rand)=1;
    t_=1-t;
    %obtain the third offspring individual
    u5=t.*v5+t_.*x(i,:);

    %save these three offspring individuals into the offspring population
    %temp_x
    temp_x(3*i-2:3*i,:)=[u1;u3;u5];
    
end
