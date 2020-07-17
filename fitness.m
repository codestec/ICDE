%evaluate all the individuals in the population p, and return its values of
%the fitness function, and the degree of the inequality and equality
%constraint violation

function [fit,g,h] = fitness(p,problem,delta,A,normalization)

%% input parameters
    % p -- the population evaluated
    % problem -- the problem tested
    % delta -- the tolarence value
    % A -- the array of parameters p,q,r
    % normalization -- the flag to indicate the way that will be seleced
%% output parameters
    % fit -- values of fitness functions and degree of constraint violation
    % g -- degree of the inequality constraint violations
    % h -- degree of the equality constraint violations

%compute the size of the population p
popsize=size(p,1);

%initialize the set fit
fit=zeros(popsize,2);

%initialize the set of degree of equality constraint violations
h=[];

%initialize the set of degree of inequality constraint violations
g=[];

%choose the problem
switch problem

    case 1

        %calculate the degree of inequality constraint violations
        g(:,1)=2*p(:,1)+2*p(:,2)+p(:,10)+p(:,11)-10;
        g(:,2)=2*p(:,1)+2*p(:,3)+p(:,10)+p(:,12)-10;
        g(:,3)=2*p(:,2)+2*p(:,3)+p(:,11)+p(:,12)-10;
        g(:,4)=-8*p(:,1)+p(:,10);
        g(:,5)=-8*p(:,2)+p(:,11);
        g(:,6)=-8*p(:,3)+p(:,12);
        g(:,7)=-2*p(:,4)-p(:,5)+p(:,10);
        g(:,8)=-2*p(:,6)-p(:,7)+p(:,11);
        g(:,9)=-2*p(:,8)-p(:,9)+p(:,12);
        
        %compute the value of the objective function
        f1=5*sum(p(:,1:4),2)-5*sum(p(:,1:4).^2,2)-sum(p(:,5:13),2);
        
        %compute the degree of constraint violations according to 'normalization'
        if normalization==1

            f2=sum(max(0,g)./(ones(popsize,1)*(max(max(0,g))+1E-30)),2)./size(g,2);

        else

            f2=sum(max(0,g),2);

        end

    case 2

        %calculate the degree of inequality constraint violations
        g(:,1)=0.75-prod(p,2);
        g(:,2)=sum(p,2)-7.5*size(p,2);
        
        %compute the value of the objective function
        f1=-abs(sum((cos(p).^4),2)-2*prod((cos(p).^2),2))./(sqrt(sum(((ones(size(p,1),1)*[1:size(p,2)]).*(p.^2)),2))+1e-30);
        
        %compute the degree of constraint violations according to 'normalization'
        if normalization==1

            f2=sum(max(0,g)./(ones(popsize,1)*(max(max(0,g))+1E-30)),2)./size(g,2);

        else

            f2=sum(max(0,g),2);

        end

    case 3

        %calculate the degree of equality constraint violations
        h(:,1)=sum(p.^2,2)-1;
        
        %compute the value of the objective function
        f1=-(10.^0.5)^10*prod(p,2);
        
        %compute the degree of constraint violations according to 'normalization'
        if normalization==1

            f2=sum(max(0,abs(h)-delta)./(ones(popsize,1)*(max(max(0,abs(h)-delta))+1E-30)),2)./(size(h,2));

        else

            f2=sum(max(0,abs(h)-delta),2);

        end

    case 4

        %calculate the degree of inequality constraint violations
        g(:,1)=+85.334407+0.0056858*p(:,2).*p(:,5)+0.0006262*p(:,1).*p(:,4)-0.0022053*p(:,3).*p(:,5)-92;
        g(:,2)=-85.334407-0.0056858*p(:,2).*p(:,5)-0.0006262*p(:,1).*p(:,4)+0.0022053*p(:,3).*p(:,5);
        g(:,3)=+80.51249+0.0071317*p(:,2).*p(:,5)+0.0029955*p(:,1).*p(:,2)+0.0021813*p(:,3).^2-110;
        g(:,4)=-80.51249-0.0071317*p(:,2).*p(:,5)-0.0029955*p(:,1).*p(:,2)-0.0021813*p(:,3).^2+90;
        g(:,5)=+9.300961+0.0047026*p(:,3).*p(:,5)+0.0012547*p(:,1).*p(:,3)+0.0019085*p(:,3).*p(:,4)-25;
        g(:,6)=-9.300961-0.0047026*p(:,3).*p(:,5)-0.0012547*p(:,1).*p(:,3)-0.0019085*p(:,3).*p(:,4)+20;
        
        %compute the value of the objective function
        f1=5.3578547*p(:,3).^2+0.8356891*p(:,1).*p(:,5)+37.293239*p(:,1)-40792.141;
        
        %compute the degree of constraint violations according to 'normalization'
        if normalization==1

            f2=sum(max(0,g)./(ones(popsize,1)*(max(max(0,g))+1E-30)),2)./size(g,2);

        else

            f2=sum(max(0,g),2);

        end

    case 5

        %calculate the degree of inequality constraint violations
        g(:,1)=-p(:,4)+p(:,3)-0.55;
        g(:,2)=-p(:,3)+p(:,4)-0.55;
        
        %calculate the degree of equality constraint violations
        h(:,1)=1000*sin(-p(:,3)-0.25)+1000*sin(-p(:,4)-0.25)+894.8-p(:,1);
        h(:,2)=1000*sin(p(:,3)-0.25)+1000*sin(p(:,3)-p(:,4)-0.25)+894.8-p(:,2);
        h(:,3)=1000*sin(p(:,4)-0.25)+1000*sin(p(:,4)-p(:,3)-0.25)+1294.8;
        
        %compute the value of the objective function
        f1=3*p(:,1)+0.000001*p(:,1).^3+2*p(:,2)+0.000002/3*p(:,2).^3;
%求最大值
% f1 = -( 3*p(:,1)+0.000001*p(:,1).^3+2*p(:,2)+0.000002/3*p(:,2).^3 ) ;
        
        %compute the degree of constraint violations according to 'normalization'
        if normalization==1

            f2=(sum(max(0,g)./(ones(popsize,1)*(max(max(0,g))+1E-30)),2)+...
                sum(max(0,abs(h)-delta)./(ones(popsize,1)*(max(max(0,abs(h)-delta))+1E-30)),2))./(size(g,2)+size(h,2));

        else

            f2=sum(max(0,g),2)+sum(max(0,abs(h)-delta),2);

        end

    case 6

        %calculate the degree of inequality constraint violations
        g(:,1)=-(p(:,1)-5).^2-(p(:,2)-5).^2+100;
        g(:,2)=(p(:,1)-6).^2+(p(:,2)-5).^2-82.81;
        
        %compute the value of the objective function
        f1=(p(:,1)-10).^3+(p(:,2)-20).^3;
        
        %compute the degree of constraint violations according to
        %'normalization'
        if normalization==1

            f2=sum(max(0,g)./(ones(popsize,1)*(max(max(0,g))+1E-30)),2)./size(g,2);

        else

            f2=sum(max(0,g),2);

        end

    case 7

        %calculate the degree of inequality constraint violations
        g(:,1)=-105+4*p(:,1)+5*p(:,2)-3*p(:,7)+9*p(:,8);
        g(:,2)=10*p(:,1)-8*p(:,2)-17*p(:,7)+2*p(:,8);
        g(:,3)=-8*p(:,1)+2*p(:,2)+5*p(:,9)-2*p(:,10)-12;
        g(:,4)=3*(p(:,1)-2).^2+4*(p(:,2)-3).^2+2*p(:,3).^2-7*p(:,4)-120;
        g(:,5)=5*p(:,1).^2+8*p(:,2)+(p(:,3)-6).^2-2*p(:,4)-40;
        g(:,6)=p(:,1).^2+2*(p(:,2)-2).^2-2*p(:,1).*p(:,2)+14*p(:,5)-6*p(:,6);
        g(:,7)=0.5*(p(:,1)-8).^2+2*(p(:,2)-4).^2+3*p(:,5).^2-p(:,6)-30;
        g(:,8)=-3*p(:,1)+6*p(:,2)+12*(p(:,9)-8).^2-7*p(:,10);
        
        %compute the value of the objective function
        f1=p(:,1).^2+p(:,2).^2+p(:,1).*p(:,2)-14*p(:,1)-16*p(:,2)+(p(:,3)-10).^2+4*(p(:,4)-5).^2+(p(:,5)-3).^2+2*(p(:,6)-1).^2+5*p(:,7).^2+7*(p(:,8)-11).^2+2*(p(:,9)-10).^2+(p(:,10)-7).^2+45;
        
        %compute the degree of constraint violations according to
        %'normalization'
        if normalization==1

            f2=sum(max(0,g)./(ones(popsize,1)*(max(max(0,g))+1E-30)),2)./size(g,2);

        else

            f2=sum(max(0,g),2);

        end

    case 8

        %calculate the degree of inequality constraint violations
        g(:,1)=p(:,1).^2-p(:,2)+1;
        g(:,2)=1-p(:,1)+(p(:,2)-4).^2;
        
        %compute the value of the objective function
        f1=-(sin(2*pi*p(:,1)).^3).*sin(2*pi*p(:,2))./((p(:,1).^3.*(p(:,1)+p(:,2)))+1e-17);
        
        %compute the degree of constraint violations according to
        %'normalization'
        if normalization==1

            f2=sum(max(0,g)./(ones(popsize,1)*(max(max(0,g))+1E-30)),2)./size(g,2);

        else

            f2=sum(max(0,g),2);

        end

    case 9

        %calculate the degree of inequality constraint violations
        g(:,1)=-127+2*p(:,1).^2+3*p(:,2).^4+p(:,3)+4*p(:,4).^2+5*p(:,5);
        g(:,2)=-282+7*p(:,1)+3*p(:,2)+10*p(:,3).^2+p(:,4)-p(:,5);
        g(:,3)=-196+23*p(:,1)+p(:,2).^2+6*p(:,6).^2-8*p(:,7);
        g(:,4)=4*p(:,1).^2+p(:,2).^2-3*p(:,1).*p(:,2)+2*p(:,3).^2+5*p(:,6)-11*p(:,7);
        
        %compute the value of the objective function
        f1=(p(:,1)-10).^2+5*(p(:,2)-12).^2+p(:,3).^4+3*(p(:,4)-11).^2+10*p(:,5).^6+7*p(:,6).^2+p(:,7).^4-4*p(:,6).*p(:,7)-10*p(:,6)-8*p(:,7);
        
        %compute the degree of constraint violations according to
        %'normalization'
        if normalization==1

            f2=sum(max(0,g)./(ones(popsize,1)*(max(max(0,g))+1E-30)),2)./size(g,2);

        else

            f2=sum(max(0,g),2);

        end

    case 10

        %calculate the degree of inequality constraint violations
        g(:,1)=-1+0.0025*(p(:,4)+p(:,6));
        g(:,2)=-1+0.0025*(p(:,5)+p(:,7)-p(:,4));
        g(:,3)=-1+0.01*(p(:,8)-p(:,5));
        g(:,4)=-p(:,1).*p(:,6)+833.33252*p(:,4)+100*p(:,1)-83333.333;
        g(:,5)=-p(:,2).*p(:,7)+1250*p(:,5)+p(:,2).*p(:,4)-1250*p(:,4);
        g(:,6)=-p(:,3).*p(:,8)+1250000+p(:,3).*p(:,5)-2500*p(:,5);
        
        %compute the value of the objective function
        f1=p(:,1)+p(:,2)+p(:,3);
        
        %compute the degree of constraint violations according to
        %'normalization'
        if normalization==1

            f2=sum(max(0,g)./(ones(popsize,1)*(max(max(0,g))+1E-30)),2)./size(g,2);

        else

            f2=sum(max(0,g),2);

        end

    case 11

        %calculate the degree of equality constraint violations
        h(:,1)=p(:,2)-p(:,1).^2;
        
        %compute the value of the objective function
        f1=p(:,1).^2+(p(:,2)-1).^2;
        
        %compute the degree of constraint violations according to
        %'normalization'
        if normalization==1

            f2=sum(max(0,abs(h)-delta)./(ones(popsize,1)*(max(max(0,abs(h)-delta))+1E-30)),2)./(size(h,2));

        else

            f2=sum(max(0,abs(h)-delta),2);

        end

    case 12

        %compute the value of the objective function
        f1=-(100-(p(:,1)-5).^2-(p(:,2)-5).^2-(p(:,3)-5).^2)/100;
        
        %calculate the degree of inequality constraint violations
        for j=1:popsize
            g(j,1)=max(min(sum((ones(9*9*9,1)*p(j,:)-A).^2,2))-0.0625,0);
        end
        
        %compute the degree of constraint violations according to
        %'normalization'
        if normalization==1

            f2=sum(max(0,g)./(ones(popsize,1)*(max(max(0,g))+1E-30)),2)./size(g,2);

        else

            f2=sum(max(0,g),2);

        end

    case 13

        %calculate the degree of equality constraint violations
        h(:,1)=p(:,1).^2+p(:,2).^2+p(:,3).^2+p(:,4).^2+p(:,5).^2-10;
        h(:,2)=p(:,2).*p(:,3)-5*p(:,4).*p(:,5);
        h(:,3)=p(:,1).^3+p(:,2).^3+1;
        
        %compute the value of the objective function
        f1=exp(p(:,1).*p(:,2).*p(:,3).*p(:,4).*p(:,5));
        
        %compute the degree of constraint violations according to
        %'normalization'
        if normalization==1

            f2=sum(max(0,abs(h)-delta)./(ones(popsize,1)*(max(max(0,abs(h)-delta))+1E-30)),2)./(size(h,2));

        else

            f2=sum(max(0,abs(h)-delta),2);

        end

    case 14

        c=[-6.089 -17.164 -34.054 -5.914 -24.721 -14.986 -24.1 -10.708 -26.662 -22.179];
        
        %calculate the degree of equality constraint violations
        h(:,1)=p(:,1)+2*p(:,2)+2*p(:,3)+p(:,6)+p(:,10)-2;
        h(:,2)=p(:,4)+2*p(:,5)+p(:,6)+p(:,7)-1;
        h(:,3)=p(:,3)+p(:,7)+p(:,8)+2*p(:,9)+p(:,10)-1;
        
        %compute the value of the objective function
        f1=sum(p.*(repmat(c,popsize,1)+log(1E-30+p./(1E-30+repmat(sum(p,2),1,10)))),2);
        
        %compute the degree of constraint violations according to
        %'normalization'
        if normalization==1

            f2=sum(max(0,abs(h)-delta)./(ones(popsize,1)*(max(max(0,abs(h)-delta))+1E-30)),2)./(size(h,2));

        else

            f2=sum(max(0,abs(h)-delta),2);

        end

    case 15

        %calculate the degree of equality constraint violations
        h(:,1)=p(:,1).^2+p(:,2).^2+p(:,3).^2-25;
        h(:,2)=8*p(:,1)+14*p(:,2)+7*p(:,3)-56;
        
        %compute the value of the objective function
        f1=1000-p(:,1).^2-2*p(:,2).^2-p(:,3).^2-p(:,1).*p(:,2)-p(:,1).*p(:,3);
        
        %compute the degree of constraint violations according to
        %'normalization'
        if normalization==1

            f2=sum(max(0,abs(h)-delta)./(ones(popsize,1)*(max(max(0,abs(h)-delta))+1E-30)),2)./(size(h,2));

        else

            f2=sum(max(0,abs(h)-delta),2);

        end

    case 16
        
        %calculate the degree of inequality constraint violations
        y1=p(:,2)+p(:,3)+41.6;
        c1=0.024*p(:,4)-4.62;
        y2=12.5./c1+12;
        c2=0.0003535*p(:,1).^2+0.5311*p(:,1)+0.08705*y2.*p(:,1);
        c3=0.052*p(:,1)+78+0.002377*y2.*p(:,1);
        y3=c2./c3;
        y4=19*y3;
        c4=0.04782*(p(:,1)-y3)+0.1956*(p(:,1)-y3).^2./p(:,2)+0.6376*y4+1.594*y3;
        c5=100*p(:,2);
        c6=p(:,1)-y3-y4;
        c7=0.950-c4./c5;
        y5=c6.*c7;
        y6=p(:,1)-y5-y4-y3;
        c8=(y5+y4)*0.995;
        y7=c8./y1;
        y8=c8/3798;
        c9=y7-0.0663*y7./y8-0.3153;
        y9=96.82./c9+0.321*y1;
        y10=1.29*y5+1.258*y4+2.29*y3+1.71*y6;
        y11 = 1.71 * p(:,1) - 0.452 * y4 + 0.580 * y3;
        c10 = 12.3 / 752.3;
        c11 = 1.75 * y2.*0.995.* p(:,1);
        c12 = 0.995 * y10 + 1998.0;
        y12 = c10*p(:,1) + (c11./ c12);
        y13 = c12 - 1.75 *y2;
        y14= 3623.0 + 64.4 *p(:,2) + 58.4 * p(:,3) + (146312.0./ (y9 + p(:,5)));
        c13 = 0.995 * y10 + 60.8 * p(:,2) + 48 * p(:,4) - 0.1121 * y14 - 5095.0;
        y15 = y13./ c13;
        y16= 148000.0 - 331000.0 * y15 + 40.0 * y13 - 61.0 .*y15.*y13;
        c14 = 2324 * y10 - 28740000 * y2;
        y17 = 14130000 - 1328.0 * y10 - 531.0 * y11 + (c14./c12);
        c15 = (y13./y15) - (y13/ 0.52);
        c16 = 1.104 - 0.72 * y15;
        c17= y9 + p(:,5);

        g(:,1)=0.28/0.72.*y5-y4;
        g(:,2)=p(:,3)-1.5*p(:,2);
        g(:,3)=3496.*y2./c12-21;
        g(:,4)=110.6+y1-62212./c17;
        g(:,5)=213.1-y1;
        g(:,6)=y1-405.23;
        g(:,7)=17.505-y2;
        g(:,8)=y2-1053.6667;
        g(:,9)=11.275-y3;
        g(:,10)=y3-35.03;
        g(:,11)=214.228-y4;
        g(:,12)=y4-665.585;
        g(:,13)=7.458-y5;
        g(:,14)=y5-584.463;
        g(:,15)=0.961-y6;
        g(:,16)=y6-265.916;
        g(:,17)=1.612-y7;
        g(:,18)=y7-7.046;
        g(:,19)=0.146-y8;
        g(:,20)=y8-0.222;
        g(:,21)=107.99-y9;
        g(:,22)=y9-273.366;
        g(:,23)=922.693-y10;
        g(:,24)=y10-1286.105;
        g(:,25)=926.832-y11;
        g(:,26)=y11-1444.046;
        g(:,27)=18.766-y12;
        g(:,28)=y12-537.141;
        g(:,29)=1072.163-y13;
        g(:,30)=y13-3247.039;
        g(:,31)=8961.448-y14;
        g(:,32)=y14-26844.086;
        g(:,33)=0.063-y15;
        g(:,34)=y15-0.386;
        g(:,35)=71084.33-y16;
        g(:,36)=-140000+y16;
        g(:,37)=2802713-y17;
        g(:,38)=y17-12146108;
        
        %compute the value of the objective function
        f1=0.000117 *y14+0.1365 + 0.00002358 * y13 + 0.000001502 * y16 + 0.0321 * y12 ...
            + 0.004324 * y5 + 0.0001 * (c15 ./ c16) + 37.48 * (y2./c12)-0.0000005843 * y17;
        
        %compute the degree of constraint violations according to
        %'normalization'
        if normalization==1

            f2=sum(max(0,g)./(ones(popsize,1)*(max(max(0,g))+1E-30)),2)./size(g,2);

        else

            f2=sum(max(0,g),2);

        end

    case 17

        %calculate the degree of equality constraint violations
        h(:,1)=-p(:,1)+300-p(:,3).*p(:,4)./131.078.*cos(1.48477-p(:,6))+0.90798.*p(:,3).^2./131.078.*cos(1.47588);
        h(:,2)=-p(:,2)-p(:,3).*p(:,4)./131.078.*cos(1.48477+p(:,6))+0.90798.*p(:,4).^2./131.078.*cos(1.47588);
        h(:,3)=-p(:,5)-p(:,3).*p(:,4)./131.078.*sin(1.48477+p(:,6))+0.90798.*p(:,4).^2./131.078.*sin(1.47588);
        h(:,4)=200-p(:,3).*p(:,4)./131.078.*sin(1.48477-p(:,6))+0.90798.*p(:,3).^2./131.078.*sin(1.47588);
        
        %compute the value of the objective function
        f1=30.*p(:,1).*(p(:,1)<300)+31.*p(:,1).*(p(:,1)>=300)+28.*p(:,2).*(p(:,2)<100)+29.*p(:,2).*(p(:,2)>=100&p(:,2)<200)+30.*p(:,2).*(p(:,2)>=200&p(:,2)<1000);
        
        %compute the degree of constraint violations according to
        %'normalization'
        if normalization==1

            f2=sum(max(0,abs(h)-delta)./(ones(popsize,1)*(max(max(0,abs(h)-delta))+1E-30)),2)./(size(h,2));

        else

            f2=sum(max(0,abs(h)-delta),2);

        end

    case 18

        
        %calculate the degree of inequality constraint violations
        g(:,1)=p(:,3).^2+p(:,4).^2-1;
        g(:,2)=p(:,9).^2-1;
        g(:,3)=p(:,5).^2+p(:,6).^2-1;
        g(:,4)=p(:,1).^2+(p(:,2)-p(:,9)).^2-1;
        g(:,5)=(p(:,1)-p(:,5)).^2+(p(:,2)-p(:,6)).^2-1;
        g(:,6)=(p(:,1)-p(:,7)).^2+(p(:,2)-p(:,8)).^2-1;
        g(:,7)=(p(:,3)-p(:,5)).^2+(p(:,4)-p(:,6)).^2-1;
        g(:,8)=(p(:,3)-p(:,7)).^2+(p(:,4)-p(:,8)).^2-1;
        g(:,9)=p(:,7).^2+(p(:,8)-p(:,9)).^2-1;
        g(:,10)=p(:,2).*p(:,3)-p(:,1).*p(:,4);
        g(:,11)=-p(:,3).*p(:,9);
        g(:,12)=p(:,5).*p(:,9);
        g(:,13)=p(:,6).*p(:,7)-p(:,5).*p(:,8);
        
        %compute the value of the objective function
        f1=-0.5*(p(:,1).*p(:,4)-p(:,2).*p(:,3)+p(:,3).*p(:,9)-p(:,5).*p(:,9)+p(:,5).*p(:,8)-p(:,6).*p(:,7));
        
        %compute the degree of constraint violations according to
        %'normalization'
        if normalization==1

            f2=sum(max(0,g)./(ones(popsize,1)*(max(max(0,g))+1E-30)),2)./size(g,2);

        else

            f2=sum(max(0,g),2);

        end

    case 19

        %calculate the degree of equality constraint violations
        a=[-16 2 0 1 0;
            0 -2 0 0.4 2;
            -3.5 0 2 0 0;
            0 -2 0 -4 -1;
            0 -9 -2 1 -2.8;
            2 0 -4 0 0;
            -1 -1 -1 -1 -1;
            -1 -2 -3 -2 -1;
            1 2 3 4 5;
            1 1 1 1 1];
        b=[-40 -2 -0.25 -4 -4 -1 -40 -60 5 1];
        c=[30 -20 -10 32 -10;
            -20 39 -6 -31 32;
            -10 -6 10 -6 -10;
            32 -31 -6 39 -20;
            -10 32 -10 -20 30];
        d=[4 8 10 6 2];
        e=[-15 -27 -36 -18 -12];
        g(:,1)=-2*sum(repmat(c(1:5,1)',popsize,1).*p(:,11:15),2)-3*d(1).*p(:,11).^2-e(1)+sum(repmat(a(1:10,1)',popsize,1).*p(:,1:10),2);
        g(:,2)=-2*sum(repmat(c(1:5,2)',popsize,1).*p(:,11:15),2)-3*d(2).*p(:,12).^2-e(2)+sum(repmat(a(1:10,2)',popsize,1).*p(:,1:10),2);
        g(:,3)=-2*sum(repmat(c(1:5,3)',popsize,1).*p(:,11:15),2)-3*d(3).*p(:,13).^2-e(3)+sum(repmat(a(1:10,3)',popsize,1).*p(:,1:10),2);
        g(:,4)=-2*sum(repmat(c(1:5,4)',popsize,1).*p(:,11:15),2)-3*d(4).*p(:,14).^2-e(4)+sum(repmat(a(1:10,4)',popsize,1).*p(:,1:10),2);
        g(:,5)=-2*sum(repmat(c(1:5,5)',popsize,1).*p(:,11:15),2)-3*d(5).*p(:,15).^2-e(5)+sum(repmat(a(1:10,5)',popsize,1).*p(:,1:10),2);
        
        %compute the value of the objective function
        f1=sum(repmat(c(1:5,1)',popsize,1).*p(:,11:15),2).*p(:,11)+sum(repmat(c(1:5,2)',popsize,1).*p(:,11:15),2).*p(:,12)...
            +sum(repmat(c(1:5,3)',popsize,1).*p(:,11:15),2).*p(:,13)+sum(repmat(c(1:5,4)',popsize,1).*p(:,11:15),2).*p(:,14)...
            +sum(repmat(c(1:5,5)',popsize,1).*p(:,11:15),2).*p(:,15)+2*sum(repmat(d,popsize,1).*p(:,11:15).^3,2)...
            -sum(repmat(b,popsize,1).*p(:,1:10),2);
        
        %compute the degree of constraint violations according to
        %'normalization'
        if normalization==1

            f2=sum(max(0,g)./(ones(popsize,1)*(max(max(0,g))+1E-30)),2)./size(g,2);

        else

            f2=sum(max(0,g),2);

        end

    case 20

        a=[0.0693 0.0577 0.05 0.2 0.26 0.55 0.06 0.1 0.12 0.18 0.1 0.09...
            0.0693 0.0577 0.05 0.2 0.26 0.55 0.06 0.1 0.12 0.18 0.1 0.09];
        b=[44.094 58.12 58.12 137.4 120.9 170.9 62.501 84.94 133.425 82.507 46.07 60.097...
            44.094 58.12 58.12 137.4 120.9 170.9 62.501 84.94 133.425 82.507 46.07 60.079];
        c=[123.7 31.7 45.7 14.7 84.7 27.7 49.7 7.1 2.1 17.7 0.85 0.64];
        d=[31.244 36.12 34.784 92.7 82.7 91.6 56.708 82.7 80.8 64.517 49.4 49.1];
        e=[0.1 0.3 0.4 0.3 0.6 0.3];
        
        %calculate the degree of inequality constraint violations
        g(:,1)=(p(:,1)+p(:,13))./(sum(p,2)+e(1));
        g(:,2)=(p(:,2)+p(:,14))./(sum(p,2)+e(2));
        g(:,3)=(p(:,3)+p(:,15))./(sum(p,2)+e(3));
        g(:,4)=(p(:,7)+p(:,19))./(sum(p,2)+e(4));
        g(:,5)=(p(:,8)+p(:,20))./(sum(p,2)+e(5));
        g(:,6)=(p(:,9)+p(:,21))./(sum(p,2)+e(6));
        
        %calculate the degree of equality constraint violations
        h(:,1)=p(:,13)./(b(13)*(sum(p(:,13:24)./repmat(b(13:24),popsize,1),2)))-c(1)*p(:,1)./(40*b(1)*(sum(p(:,1:12)./repmat(b(1:12),popsize,1),2)));
        h(:,2)=p(:,14)./(b(14)*(sum(p(:,13:24)./repmat(b(13:24),popsize,1),2)))-c(2)*p(:,2)./(40*b(2)*(sum(p(:,1:12)./repmat(b(1:12),popsize,1),2)));
        h(:,3)=p(:,15)./(b(15)*(sum(p(:,13:24)./repmat(b(13:24),popsize,1),2)))-c(3)*p(:,3)./(40*b(3)*(sum(p(:,1:12)./repmat(b(1:12),popsize,1),2)));
        h(:,4)=p(:,16)./(b(16)*(sum(p(:,13:24)./repmat(b(13:24),popsize,1),2)))-c(4)*p(:,4)./(40*b(4)*(sum(p(:,1:12)./repmat(b(1:12),popsize,1),2)));
        h(:,5)=p(:,17)./(b(17)*(sum(p(:,13:24)./repmat(b(13:24),popsize,1),2)))-c(5)*p(:,5)./(40*b(5)*(sum(p(:,1:12)./repmat(b(1:12),popsize,1),2)));
        h(:,6)=p(:,18)./(b(18)*(sum(p(:,13:24)./repmat(b(13:24),popsize,1),2)))-c(6)*p(:,6)./(40*b(6)*(sum(p(:,1:12)./repmat(b(1:12),popsize,1),2)));
        h(:,7)=p(:,19)./(b(19)*(sum(p(:,13:24)./repmat(b(13:24),popsize,1),2)))-c(7)*p(:,7)./(40*b(7)*(sum(p(:,1:12)./repmat(b(1:12),popsize,1),2)));
        h(:,8)=p(:,20)./(b(20)*(sum(p(:,13:24)./repmat(b(13:24),popsize,1),2)))-c(8)*p(:,8)./(40*b(8)*(sum(p(:,1:12)./repmat(b(1:12),popsize,1),2)));
        h(:,9)=p(:,21)./(b(21)*(sum(p(:,13:24)./repmat(b(13:24),popsize,1),2)))-c(9)*p(:,9)./(40*b(9)*(sum(p(:,1:12)./repmat(b(1:12),popsize,1),2)));
        h(:,10)=p(:,22)./(b(22)*(sum(p(:,13:24)./repmat(b(13:24),popsize,1),2)))-c(10)*p(:,10)./(40*b(10)*(sum(p(:,1:12)./repmat(b(1:12),popsize,1),2)));
        h(:,11)=p(:,23)./(b(23)*(sum(p(:,13:24)./repmat(b(13:24),popsize,1),2)))-c(11)*p(:,11)./(40*b(11)*(sum(p(:,1:12)./repmat(b(1:12),popsize,1),2)));
        h(:,12)=p(:,24)./(b(24)*(sum(p(:,13:24)./repmat(b(13:24),popsize,1),2)))-c(12)*p(:,12)./(40*b(12)*(sum(p(:,1:12)./repmat(b(1:12),popsize,1),2)));
        h(:,13)=sum(p,2)-1;
        h(:,14)=sum(p(:,1:12)./repmat(d(1:12),popsize,1),2)+0.7302*530*14.7/40*sum(p(:,13:24)./repmat(b(13:24),popsize,1),2)-1.671;
        
        %compute the value of the objective function
        f1=sum(repmat(a,popsize,1).*p,2);
        
        %compute the degree of constraint violations according to
        %'normalization'
        if normalization==1

            f2=(sum(max(0,g)./(ones(popsize,1)*(max(max(0,g))+1E-30)),2)+...
                sum(max(0,abs(h)-delta)./(ones(popsize,1)*(max(max(0,abs(h)-delta))+1E-30)),2))./(size(g,2)+size(h,2));

        else

            f2=sum(max(0,g),2)+sum(max(0,abs(h)-delta),2);

        end

    case 21

        %calculate the degree of inequality constraint violations
        g(:,1)=-p(:,1)+35*p(:,2).^0.6+35*p(:,3).^0.6;
        
        %calculate the degree of equality constraint violations
        h(:,1)=-300*p(:,3)+7500*p(:,5)-7500*p(:,6)-25*p(:,4).*p(:,5)+25*p(:,4).*p(:,6)+p(:,3).*p(:,4);
        h(:,2)=100*p(:,2)+155.365*p(:,4)+2500*p(:,7)-p(:,2).*p(:,4)-25*p(:,4).*p(:,7)-15536.5;
        h(:,3)=-p(:,5)+log(-p(:,4)+900);
        h(:,4)=-p(:,6)+log(p(:,4)+300);
        h(:,5)=-p(:,7)+log(-2*p(:,4)+700);
        
        %calculate the degree of equality constraint violations
        f1=p(:,1);
        
        %compute the degree of constraint violations according to
        %'normalization'
        if normalization==1

            f2=(sum(max(0,g)./(ones(popsize,1)*(max(max(0,g))+1E-30)),2)+...
                sum(max(0,abs(h)-delta)./(ones(popsize,1)*(max(max(0,abs(h)-delta))+1E-30)),2))./(size(g,2)+size(h,2));

        else

            f2=sum(max(0,g),2)+sum(max(0,abs(h)-delta),2);

        end

    case 22

        %calculate the degree of inequality constraint violations
        g(:,1)=-p(:,1)+p(:,2).^0.6+p(:,3).^0.6+p(:,4).^0.6;
        
        %calculate the degree of equality constraint violations
        h(:,1)=p(:,5)-100000*p(:,8)+10^7;
        h(:,2)=p(:,6)+100000*p(:,8)-100000*p(:,9);
        h(:,3)=p(:,7)+100000*p(:,9)-5*10^7;
        h(:,4)=p(:,5)+100000*p(:,10)-3.3*10^7;
        h(:,5)=p(:,6)+100000*p(:,11)-4.4*10^7;
        h(:,6)=p(:,7)+100000*p(:,12)-6.6*10^7;
        h(:,7)=p(:,5)-120*p(:,2).*p(:,13);
        h(:,8)=p(:,6)-80*p(:,3).*p(:,14);
        h(:,9)=p(:,7)-40*p(:,4).*p(:,15);
        h(:,10)=p(:,8)-p(:,11)+p(:,16);
        h(:,11)=p(:,9)-p(:,12)+p(:,17);
        h(:,12)=-p(:,18)+log(p(:,10)-100);
        h(:,13)=-p(:,19)+log(-p(:,8)+300);
        h(:,14)=-p(:,20)+log(p(:,16));
        h(:,15)=-p(:,21)+log(-p(:,9)+400);
        h(:,16)=-p(:,22)+log(p(:,17));
        h(:,17)=-p(:,8)-p(:,10)+p(:,13).*p(:,18)-p(:,13).*p(:,19)+400;
        h(:,18)=p(:,8)-p(:,9)-p(:,11)+p(:,14).*p(:,20)-p(:,14).*p(:,21)+400;
        h(:,19)=p(:,9)-p(:,12)-4.60517*p(:,15)+p(:,15).*p(:,22)+100;
        
        %calculate the degree of equality constraint violations
        f1=p(:,1);
        
        %compute the degree of constraint violations according to
        %'normalization'
        if normalization==1

            f2=(sum(max(0,g)./(ones(popsize,1)*(max(max(0,g))+1E-30)),2)+...
                sum(max(0,abs(h)-delta)./(ones(popsize,1)*(max(max(0,abs(h)-delta))+1E-30)),2))./(size(g,2)+size(h,2));

        else

            f2=sum(max(0,g),2)+sum(max(0,abs(h)-delta),2);

        end

    case 23

        %calculate the degree of inequality constraint violations
        g(:,1)=p(:,9).*p(:,3)+0.02.*p(:,6)-0.025.*p(:,5);
        g(:,2)=p(:,9).*p(:,4)+0.02.*p(:,7)-0.015.*p(:,8);
        
        %calculate the degree of equality constraint violations
        h(:,1)=p(:,1)+p(:,2)-p(:,3)-p(:,4);
        h(:,2)=0.03.*p(:,1)+0.01.*p(:,2)-p(:,9).*(p(:,3)+p(:,4));
        h(:,3)=p(:,3)+p(:,6)-p(:,5);
        h(:,4)=p(:,4)+p(:,7)-p(:,8);
        
        %calculate the degree of equality constraint violations
        f1=-9.*p(:,5)-15.*p(:,8)+6.*p(:,1)+16.*p(:,2)+10.*(p(:,6)+p(:,7));
        
        %compute the degree of constraint violations according to
        %'normalization'
        if normalization==1

            f2=(sum(max(0,g)./(ones(popsize,1)*(max(max(0,g))+1E-30)),2)+...
                sum(max(0,abs(h)-delta)./(ones(popsize,1)*(max(max(0,abs(h)-delta))+1E-30)),2))./(size(g,2)+size(h,2));

        else

            f2=sum(max(0,g),2)+sum(max(0,abs(h)-delta),2);

        end

    case 24

        %calculate the degree of inequality constraint violations
        g(:,1)=-2*p(:,1).^4+8*p(:,1).^3-8*p(:,1).^2+p(:,2)-2;
        g(:,2)=-4*p(:,1).^4+32*p(:,1).^3-88*p(:,1).^2+96*p(:,1)+p(:,2)-36;
        
        %calculate the degree of equality constraint violations
        f1=-p(:,1)-p(:,2);
        
        %compute the degree of constraint violations according to
        %'normalization'
        if normalization==1

            f2=sum(max(0,g)./(ones(popsize,1)*(max(max(0,g))+1E-30)),2)./size(g,2);

        else

            f2=sum(max(0,g),2);

        end

end

%save the values obtained
fit(1:popsize,1)=f1;
fit(1:popsize,2)=f2;
