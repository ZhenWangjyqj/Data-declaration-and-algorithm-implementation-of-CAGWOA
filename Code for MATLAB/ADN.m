function [bestPositions,Destination_fitness,FEs,L] = ADN(bestPositions,Destination_fitness,Padn,L0,L,Lmin,Rou,dim,fobj,FEs)

Nb=zeros(2*dim,dim);
Nb_fit=zeros(2*dim,1);
BestPosition=bestPositions;
BestFitness=Destination_fitness;

if rand()<Padn
    for k=1:dim
        for d=1:dim
            if k==d
                Nb(2*k-1,d)=bestPositions(d)+L;
                Nb(2*k,d)=bestPositions(d)-L;
            else
                Nb(2*k-1,d)=bestPositions(d);
                Nb(2*k,d)=bestPositions(d);
            end
        end
        Nb_fit(2*k-1)=fobj(Nb(2*k-1,:));
        Nb_fit(2*k)=fobj(Nb(2*k,:));
        FEs=FEs+2;
    end
    [bst_val,bst_id]=min(Nb_fit);
    if bst_val<BestFitness
        BestFitness=bst_val;
        BestPosition=Nb(bst_id,:);
    else
        L=L*Rou;
        if L<Lmin
            L=L0*rand();
        end
    end
end
if (BestFitness<Destination_fitness)
    bestPositions = BestPosition;
    Destination_fitness = BestFitness;
end
end
