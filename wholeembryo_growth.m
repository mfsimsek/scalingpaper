clear;
SomiteSize=5; % Initial Somite Size in Cell #
Vg=[0.084    0.084    0.084 0.045	0.044	0.07	0.038	0.032	0.032	0.033	0.055	0.069	0.063	0.06	0.064	0.054	0.051	0.043	0.019	0.016	0.007	0.023	0.044	0.032];
dataY=[49
    48
    47.9
    47.2
    48.8
    49.8
    49.2
    50.9
    50.2
    49.8
    45.5
    44.6
    41
    36.3
    33.4
    29.9
    29.5
    29.4
    28.1
    ];

K=85; % PSM Size in Cell #
tailbud=15; % # of Tailbud Cells
T=71*length(Vg); % Simulation Duration
PSM=1:K; % Position Vector
time=1:T; % Time Vector

PPMatrices=ones(2,2*K,T);

for c=1:2
    if c==1
        cond=0;
    else
        cond=42;
    end
    lock=0;
    pbd=-0.2;
    rd=-0.01; % FGF RNA Decay Rate
    pfp=0.5; % Free FGF Protein Translation Rate
    pfd=-0.1; % Free FGF Protein Decay Rate
    prp=0.1; % FGF Repressor Protein Translation Rate
    prd=-1; % FGF Repressor Protein Decay Rate
    pd=-0.1; % Un-phosphorylated Signaling Protein Decay Rate
    pp=10; % Un-phosphorylated Signaling Protein Translation Rate
    rp=5; % FGF RNA Transcription Rate
    trlatedel=14; % FGF Protein Translation Time Delay
    interndel=14; % Bound Receptor Internalization Time Delay
    reprdel=28; % Repressor Protein Translation Time Delay
    pb=0.01;
    ppd=0.03;
    ppp=0.015;
    Diff=1.8;
    
    
    
    r0=0; % Initial FGF RNA
    pf0=0; % Initial Free FGF Protein
    pb0=0; % Initial Bound FGF Protein
    pr0=0; % Initial FGF Repressor Protein
    pp0=0; % Initial Phosphorylated Signaling Protein
    p0=0; % Initial Un-phosphorylated Signaling Protein
    
    RMatrix=r0*ones(K,1)*ones(1,T); % 4D FGF RNA Matrix
    PFMatrix=pf0*ones(K,1)*ones(1,T); % 4D FGF Protein Matrix
    PBMatrix=pb0*ones(K,1)*ones(1,T); % 4D Bound FGF Receptor Protein Matrix
    PRMatrix=pr0*ones(K,1)*ones(1,T); % 4D FGF Repressor Protein Matrix
    PMatrix=p0*ones(K,1)*ones(1,T); % 4D Un-phosphorylated Signaling Protein Matrix
    PPMatrix=pp0*ones(2*K,1)*ones(1,T); % 4D Phosphorylated Signaling Protein Matrix
    StatMatrix=PPMatrix;
    
    ry0matrix=r0*ones(K,1);
    pfy0matrix=pf0*ones(K,1);
    pby0matrix=pb0*ones(K,1);
    pry0matrix=pr0*ones(K,1);
    ppy0matrix=pp0*ones(K,1);
    y0matrix=p0*ones(K,1);
    
    for ti=0:T-1
        %Euler's Method for dy/dt=f(y,t)
        clear y t
        N=100;
        for k=1:tailbud
            tf=ti+1; % initial and final time
            h=(tf-ti)/N; % Time step
            ry=zeros(N+1,1);
            pfy=zeros(N+1,1);
            pby=zeros(N+1,1);
            pry=zeros(N+1,1);
            py=zeros(N+1,1);
            ppy=zeros(N+1,1);
            
            if ti==0
                ry(1)=ry0matrix(k);
                pfy(1)=pfy0matrix(k);
                pby(1)=pby0matrix(k);
                pry(1)=pry0matrix(k);
                py(1)=y0matrix(k);
                ppy(1)=ppy0matrix(k);
            else
                ry(1)=RMatrix(k,ti);
                pfy(1)=PFMatrix(k,ti);
                pby(1)=PBMatrix(k,ti);
                pry(1)=PRMatrix(k,ti);
                py(1)=PMatrix(k,ti);
                ppy(1)=PPMatrix(k,ti);
            end
            for n=1:N
                ry(n+1)=ry(n)+h*rp+h*rd*ry(n);
                
                if ti>interndel&&ti>trlatedel
                    pfy(n+1)=pfy(n)+h*(pfp*RMatrix(k,ti-trlatedel)+pfd*pfy(n)-pb*PFMatrix(k,ti-interndel));
                    pby(n+1)=pby(n)+h*(pb*PFMatrix(k,ti-interndel)+pbd*pby(n));
                else
                    pfy(n+1)=pfy(n);
                    pby(n+1)=pby(n);
                end
                if ti>reprdel
                    pry(n+1)=pry(n)+h*(prp*PPMatrix(k,ti-reprdel)+prd*pry(n));
                else
                    pry(n+1)=pry(n);
                end
                ppy(n+1)=ppy(n)+h*(ppp*pby(n)*py(n)+pd*ppy(n)-ppd*pry(n)*ppy(n));
                py(n+1)=py(n)+h*(pp+ppd*pry(n)*ppy(n)+pd*py(n)-ppp*pby(n)*py(n));
                
            end
            
            ry0=ry(N+1);
            pfy0=pfy(N+1);
            pby0=pby(N+1);
            pry0=pry(N+1);
            py0=py(N+1);
            ppy0=ppy(N+1);
            
            ry0matrix(k)=ry0;
            pfy0matrix(k)=pfy0;
            pby0matrix(k)=pby0;
            pry0matrix(k)=pry0;
            y0matrix(k)=py0;
            ppy0matrix(k)=ppy0;
            
            RMatrix(k,tf)=ry0;
            PFMatrix(k,tf)=pfy0;
            PBMatrix(k,tf)=pby0;
            PRMatrix(k,tf)=pry0;
            PMatrix(k,tf)=py0;
            PPMatrix(k,tf)=ppy0;
        end
        
    end
    
    r0=ry0matrix(tailbud);
    pf0=pfy0matrix(tailbud);
    pb0=pby0matrix(tailbud);
    pr0=pry0matrix(tailbud);
    pp0=ppy0matrix(tailbud);
    p0=y0matrix(tailbud);
    
    for k=tailbud+1:K
        N=100; % Number of time steps
        ry0=r0;
        pfy0=pf0;
        pby0=pb0;
        pry0=pr0;
        py0=p0;
        ppy0=pp0;
        
        past=floor((k-tailbud)/Vg(2));
        for ti=1:past
            %Euler's Method for dy/dt=f(y,t)
            clear y t
            tf=ti+1; % initial and final time
            h=(tf-ti)/N; % Time step
            ry=zeros(N+1,1);
            pfy=zeros(N+1,1);
            pby=zeros(N+1,1);
            pry=zeros(N+1,1);
            py=zeros(N+1,1);
            ppy=zeros(N+1,1);
            
            ry(1)=ry0;
            pfy(1)=pfy0;
            pby(1)=pby0;
            pry(1)=pry0;
            py(1)=py0;
            ppy(1)=ppy0;
            
            for n=1:N
                ry(n+1)=ry(n)+h*rd*ry(n);
                
                if ti>interndel&&ti>trlatedel
                    pfy(n+1)=pfy(n)+h*(pfp*RMatrix(k,ti-trlatedel)+pfd*pfy(n)-pb*PFMatrix(k,ti-interndel));
                    pby(n+1)=pby(n)+h*(pb*PFMatrix(k,ti-interndel)+pbd*pby(n));
                else
                    pfy(n+1)=pfy(n);
                    pby(n+1)=pby(n);
                end
                
                if ti>reprdel
                    pry(n+1)=pry(n)+h*(prp*PPMatrix(k,ti-reprdel)+prd*pry(n));
                else
                    pry(n+1)=pry(n);
                end
                ppy(n+1)=ppy(n)+h*(ppp*pby(n)*py(n)+pd*ppy(n)-ppd*pry(n)*ppy(n));
                py(n+1)=py(n)+h*(pp+ppd*pry(n)*ppy(n)+pd*py(n)-ppp*pby(n)*py(n));
                
            end
            ry0=ry(N+1);
            pfy0=pfy(N+1);
            pby0=pby(N+1);
            pry0=pry(N+1);
            py0=py(N+1);
            ppy0=ppy(N+1);
            
            ry0matrix(k)=ry0;
            pfy0matrix(k)=pfy0;
            pby0matrix(k)=pby0;
            pry0matrix(k)=pry0;
            y0matrix(k)=py0;
            ppy0matrix(k)=ppy0;
            
            RMatrix(k,tf)=ry0;
            PFMatrix(k,tf)=pfy0;
            PBMatrix(k,tf)=pby0;
            PRMatrix(k,tf)=pry0;
            PMatrix(k,tf)=py0;
            PPMatrix(k,tf)=ppy0;
        end
    end
    
    for ti=0:T-1
        %Euler's Method for dy/dt=f(y,t)
        clear y t
        
        for k=1
            tf=ti+1; % initial and final time
            h=(tf-ti)/N; % Time step
            ry=zeros(N+1,1);
            pfy=zeros(N+1,1);
            pby=zeros(N+1,1);
            pry=zeros(N+1,1);
            py=zeros(N+1,1);
            ppy=zeros(N+1,1);
            
            if ti==0
                ry(1)=ry0matrix(k);
                pfy(1)=pfy0matrix(k);
                pby(1)=pby0matrix(k);
                pry(1)=pry0matrix(k);
                py(1)=y0matrix(k);
                ppy(1)=ppy0matrix(k);
            else
                ry(1)=RMatrix(k,ti);
                pfy(1)=PFMatrix(k,ti);
                pby(1)=PBMatrix(k,ti);
                pry(1)=PRMatrix(k,ti);
                py(1)=PMatrix(k,ti);
                ppy(1)=PPMatrix(k,ti);
            end
            for n=1:N
                ry(n+1)=ry(n)+h*(rp+rd*ry(n));
                if ti>trlatedel&&ti>interndel
                    pfy(n+1)=pfy(n)+h*(pfp*RMatrix(k,ti-trlatedel)+pfd*pfy(n)-pb*PFMatrix(k,ti-interndel)+Diff*(PFMatrix(k+1,ti)-pfy(n)));
                else
                    pfy(n+1)=pfy(n);
                end
                if ti>interndel
                    pby(n+1)=pby(n)+h*(pb*PFMatrix(k,ti-interndel)+pbd*pby(n));
                else
                    pby(n+1)=pby(n);
                end
                if ti>reprdel
                    pry(n+1)=pry(n)+h*(prp*PPMatrix(k,ti-reprdel)+prd*pry(n));
                else
                    pry(n+1)=pry(n);
                end
                ppy(n+1)=ppy(n)+h*(ppp*pby(n)*py(n)+pd*ppy(n)-ppd*pry(n)*ppy(n));
                py(n+1)=py(n)+h*(pp+ppd*pry(n)*ppy(n)+pd*py(n)-ppp*pby(n)*py(n));
                
            end
            
            ry0=ry(N+1);
            pfy0=pfy(N+1);
            pby0=pby(N+1);
            pry0=pry(N+1);
            py0=py(N+1);
            ppy0=ppy(N+1);
            
            ry0matrix(k)=ry0;
            pfy0matrix(k)=pfy0;
            pby0matrix(k)=pby0;
            pry0matrix(k)=pry0;
            y0matrix(k)=py0;
            ppy0matrix(k)=ppy0;
            
            RMatrix(k,tf)=ry0;
            PFMatrix(k,tf)=pfy0;
            PBMatrix(k,tf)=pby0;
            PRMatrix(k,tf)=pry0;
            PMatrix(k,tf)=py0;
            PPMatrix(k,tf)=ppy0;
        end
        
        for k=2:K-1
            tf=ti+1; % initial and final time
            h=(tf-ti)/N; % Time step
            ry=zeros(N+1,1);
            pfy=zeros(N+1,1);
            pby=zeros(N+1,1);
            pry=zeros(N+1,1);
            py=zeros(N+1,1);
            ppy=zeros(N+1,1);
            
            if ti==0
                ry(1)=ry0matrix(k);
                pfy(1)=pfy0matrix(k);
                pby(1)=pby0matrix(k);
                pry(1)=pry0matrix(k);
                py(1)=y0matrix(k);
                ppy(1)=ppy0matrix(k);
            else
                ry(1)=RMatrix(k,ti);
                pfy(1)=PFMatrix(k,ti);
                pby(1)=PBMatrix(k,ti);
                pry(1)=PRMatrix(k,ti);
                py(1)=PMatrix(k,ti);
                ppy(1)=PPMatrix(k,ti);
            end
            for n=1:N
                if k<=tailbud
                    ry(n+1)=ry(n)+h*(rp+rd*ry(n));
                else
                    ry(n+1)=ry(n)+h*rd*ry(n);
                end
                if ti>trlatedel&&ti>interndel
                    pfy(n+1)=pfy(n)+h*(pfp*RMatrix(k,ti-trlatedel)+pfd*pfy(n)-pb*PFMatrix(k,ti-interndel)+Diff*(PFMatrix(k+1,ti)+PFMatrix(k-1,ti)-2*pfy(n)));
                else
                    pfy(n+1)=pfy(n);
                end
                if ti>interndel
                    pby(n+1)=pby(n)+h*(pb*PFMatrix(k,ti-interndel)+pbd*pby(n));
                else
                    pby(n+1)=pby(n);
                end
                if ti>reprdel
                    pry(n+1)=pry(n)+h*(prp*PPMatrix(k,ti-reprdel)+prd*pry(n));
                else
                    pry(n+1)=pry(n);
                end
                ppy(n+1)=ppy(n)+h*(ppp*pby(n)*py(n)+pd*ppy(n)-ppd*pry(n)*ppy(n));
                py(n+1)=py(n)+h*(pp+ppd*pry(n)*ppy(n)+pd*py(n)-ppp*pby(n)*py(n));
                
            end
            
            ry0=ry(N+1);
            pfy0=pfy(N+1);
            pby0=pby(N+1);
            pry0=pry(N+1);
            py0=py(N+1);
            ppy0=ppy(N+1);
            
            ry0matrix(k)=ry0;
            pfy0matrix(k)=pfy0;
            pby0matrix(k)=pby0;
            pry0matrix(k)=pry0;
            y0matrix(k)=py0;
            ppy0matrix(k)=ppy0;
            
            RMatrix(k,tf)=ry0;
            PFMatrix(k,tf)=pfy0;
            PBMatrix(k,tf)=pby0;
            PRMatrix(k,tf)=pry0;
            PMatrix(k,tf)=py0;
            PPMatrix(k,tf)=ppy0;
        end
        i=ceil(tf/71);
        if mod(tf,floor(1/Vg(i)))==0
            StatMatrix(:,tf+1-floor(1/Vg(i)):tf)=PPMatrix(:,tf+1-floor(1/Vg(i)):tf);
            RMatrix=circshift(RMatrix,[1 0]);
            PFMatrix=circshift(PFMatrix,[1 0]);
            PBMatrix=circshift(PBMatrix,[1 0]);
            PRMatrix=circshift(PRMatrix,[1 0]);
            PMatrix=circshift(PMatrix,[1 0]);
            PPMatrix=circshift(PPMatrix,[1 0]);
            PFMatrix(1,:)=PFMatrix(2,:);
            RMatrix(1,:)=RMatrix(2,:);
            PRMatrix(1,:)=PRMatrix(2,:);
            PBMatrix(1,:)=PBMatrix(2,:);
            PMatrix(1,:)=PMatrix(2,:);
            PPMatrix(1,:)=PPMatrix(2,:);
            %%if ((cond>40&&ti>T-18.5*71)||ti>T-10)&&lock==0 %stepwise drug effect
            if ((cond>40&&ti>T-18.5*71)||ti>T-10)&&lock<10 % gradual drug effect
                
                %%lock=1; %stepwise drug effect
                lock=lock+1; % gradual drug effect
                switch cond
                    case {22,42}
                        %%ppp=ppp/2; %stepwise drug effect
                        ppp=ppp*exp(-log(1.2)*lock/10); % gradual drug effect
                    case {23,43}
                        ppp=ppp/3;
                    case {24,44}
                        ppp=ppp/4;
                    case {25,45}
                        ppp=ppp/5;
                    case {26,46}
                        ppp=ppp/10;
                    case 32
                        ppd=ppd/2;
                    case 33
                        ppd=ppd/3;
                    case 34
                        ppd=ppd/4;
                    case 35
                        ppd=ppd/5;
                    case 36
                        ppd=ppd/10;
                end
            end
            
        end
        PFMatrix(K,tf)=0;
        RMatrix(K,tf)=0;
        PRMatrix(K,tf)=0;
        PBMatrix(K,tf)=0;
        PMatrix(K,tf)=PMatrix(K-1,tf);
        PPMatrix(K:2*K,tf)=0;
        
    end
    
    PPMatrices(c,:,:)=PPMatrix;
end

window=53;
cyc=length(Vg)+1-floor(window/71);
klist=zeros(cyc-1,2);
T=71*length(Vg);
threshold=0.25;

for c=1:2
window=53;
aveklist=zeros(cyc-1,3);
    for w=0:2
        window=window+3;
    zmatrix=zeros(2*K-tailbud,T-window);
    for t=1:T-window
         zmatrix(:,t)=PPMatrices(c,tailbud+1:end,t+window)-0.02*PPMatrices(c,tailbud+1,t+window);
    end
    xnew = linspace(1,2*K-tailbud,8*(2*K-tailbud)+1);
    znew = interp1(zmatrix,xnew,'linear');
    Length=size(znew,1);
    
    DiffpFGFMatrix=zeros(Length,T-window);    
    for t=1:T-window
        for k=12:Length-4
            DiffpFGFMatrix(k,t)=sum(znew(k-11:k-4,t))/sum(znew(k-3:k+4,t))-1; % Normalized Delta Neighbor with Time Delay t0
        end
    end
    
    frontlist=zeros(cyc,1);
    
    aveklist(1,w+1)=48;
    idx =find(abs(DiffpFGFMatrix(12:Length,72)-threshold)<threshold/500); % Index of coordinate in array
    if isempty(idx)
        idx =find(abs(DiffpFGFMatrix(12:Length,72)-threshold)<threshold/100); % Index of coordinate in array
        if isempty(idx)
            idx =find(abs(DiffpFGFMatrix(12:Length,72)-threshold)<threshold/10); % Index of coordinate in array
            if isempty(idx)
                frontlist(4)=frontlist(3);
            else
                frontlist(4)=idx(end);
            end
        else
            frontlist(4)=idx(end);
        end
    else
        frontlist(4)=idx(end);
    end
    
    frontlist(3)=frontlist(4)+aveklist(1,w+1);
    frontlist(2)=frontlist(3)+aveklist(1,w+1);
    frontlist(1)=frontlist(2)+aveklist(1,w+1);
    
    
    for i=3:cyc-1
        idx =find(abs(DiffpFGFMatrix(12:frontlist(i),(i-1)*71+1)-threshold)<threshold/500); % Index of coordinate in array
        if isempty(idx)
            idx =find(abs(DiffpFGFMatrix(12:frontlist(i),(i-1)*71+1)-threshold)<threshold/100); % Index of coordinate in array
            if isempty(idx)
                idx =find(abs(DiffpFGFMatrix(12:frontlist(i),(i-1)*71+1)-threshold)<threshold/10); % Index of coordinate in array
                if isempty(idx)
                    frontlist(i+1)=frontlist(i);
                else
                    frontlist(i+1)=idx(end);
                end
            else
                frontlist(i+1)=idx(end);
            end
        else
            frontlist(i+1)=idx(end);
        end
    end
    
    for i=0:cyc-3
        aveklist(i+2,w+1)=frontlist(i+2)-frontlist(i+3);
        if aveklist(i+2,w+1)<8
            aveklist(i+2,w+1)=8;
        end
        frontlist(i+3)=frontlist(i+2)-aveklist(i+2,w+1);
    end
    
    end
    klist(:,c)=mean(aveklist,2);
    
    
end
klist(1,:)=[];
klist(2,:)=klist(1,:);
plot(1:length(dataY)-2,dataY(3:end),1:cyc-8,klist(1:cyc-8,1),1:cyc-8,klist(1:cyc-8,2));


[h,p]=ttest2(dataY,klist.','Tail','both');
fid1=fopen('parameter.txt','a');
fprintf(fid1,'\r\n');
fprintf(fid1,'%.3f\t',[fit cond h p]);
fclose(fid1);


fid2=fopen('frontlist.txt','a');
fprintf(fid2,'\r\n');
fprintf(fid2,'%d\t',frontlist);
fclose(fid2);
fid3=fopen('klist.txt','a');
fprintf(fid3,'\r\n');
fprintf(fid3,'%d\t',klist);
fclose(fid3);