%% CAGWOA
function [Leader_pos,Convergence_curve]=CAGWOA(SearchAgents_no,MaxFEs,lb,ub,dim,fobj)

% initialize position vector and score for the leader
Leader_pos=zeros(1,dim);
Leader_score=inf; %change this to -inf for maximization problems

%% parameters in ADN
L0=0.35;
L=L0;
Lmin=1e-8;
Padn=0.2;
Rou=0.4;

%% COSI
Ps = lhsdesign(SearchAgents_no,dim);
Positions = zeros(SearchAgents_no,dim);
for i=1:SearchAgents_no
    for j=1:dim
        b=rand;
        Positions(i,j)=cos(pi*(1/2-Ps(i,j))*b)*(ub-lb)+lb;
        if  Positions(i,j)>ub||Positions(i,j)<lb
            Positions(i,j)=lb+(ub-lb)*rand;
        end
    end
end

Convergence_curve=[];
FEs=0;
t=1;

% Main loop
while  FEs < MaxFEs
    for i=1:size(Positions,1)
        % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        % Calculate objective function for each search agent
        fitness=fobj(Positions(i,:));
        FEs=FEs+1;
        % Update the leader
        if fitness<Leader_score % Change this to > for maximization problem
            Leader_score=fitness; % Update alpha
            Leader_pos=Positions(i,:);
        end
    end

    a=2-FEs*((2)/MaxFEs);
    a2=-1+FEs*((-1)/MaxFEs);

    % Update the Position of search agents
    for i=1:size(Positions,1)
        r1=rand();
        r2=rand();
        A=2*a*r1-a;
        C=2*r2;
        b=1;
        l=(a2-1)*rand+1;
        p = rand();
        abc=2*exp(FEs/MaxFEs);
        for j=1:size(Positions,2)
            if p<0.5
                %% GS
                if(tan(pi*(rand-0.5))>(1-FEs/MaxFEs))
                    rnum=ceil(SearchAgents_no*rand);while rnum==i,rnum=ceil(SearchAgents_no*rand);end;
                    rnum1=ceil(SearchAgents_no*rand);while rnum1==i || rnum1==rnum,rnum1=ceil(SearchAgents_no*rand);end;
                    %  Positions(i,j)=bcd*Leader_pos(j)+2*abc*rand*(Positions(rnum,j)-Positions(rnum1,j));
                    Positions(i,j)=Leader_pos(j)+2*abc*rand*(Positions(rnum,j)-Positions(rnum1,j));
                end
                if abs(A)>=1
                    rand_leader_index = floor(SearchAgents_no*rand()+1);
                    X_rand = Positions(rand_leader_index, :);
                    D_X_rand=abs(C*X_rand(j)-Positions(i,j));
                    Positions(i,j)=X_rand(j)-A*D_X_rand;
                else
                    D_Leader=abs(C*Leader_pos(j)-Positions(i,j));
                    Positions(i,j)=Leader_pos(j)-A*D_Leader;
                end
            else
                distance2Leader=abs(Leader_pos(j)-Positions(i,j));
                Positions(i,j)=distance2Leader*exp(b.*l).*cos(l.*2*pi)+Leader_pos(j);
            end
        end
    end

    %% ADN
    [Leader_pos,Leader_score,FEs,L] = ADN(Leader_pos,Leader_score,Padn,L0,L,Lmin,Rou,dim,fobj,FEs);
    Convergence_curve(t)=Leader_score;
    t=t+1;
end
end
