clear;
dataX=[455 415 375 335 305 285 270 259 249];
dataY=[40 30 20 15 11 10 9 8 9];

pb=0.3;
ppd=0.105;
ppp=0.015;
Diff=1.0;

SomiteSize=5; % Initial Somite Size in Cell #
Vg=0.07; % Tailbud Growth Speed

pbd=-0.5;
rd=-0.01; % Ligand RNA Decay Rate
pfp=2.5; % Free Ligand Protein Translation Rate
pfd=-0.5; % Free Ligand Protein Decay Rate


prp=0.5; % Repressor Protein Translation Rate
prd=-0.5; % Repressor Protein Decay Rate


pd=-0.5; % Inactive Signaling Protein Decay Rate
pp=70; % Inactive Signaling Protein Translation Rate
rp=5; % Ligand RNA Transcription Rate

trlatedel=14; % Ligand Protein Translation Time Delay
actdel=14; % Bound Receptor Cascade Activation Time Delay
reprdel=28; % Repressor Protein Translation Time Delay

ppp_hs=ppp/10;
ppp_norm=ppp;
hstype=2; % type 1 long well, type 2 short well
hs=0;
well1=70;
Ths=300;

switch hs
    case 0
        well2=well1;
    case 1
        well2=well1+13;
    case 2
        well2=well1+4;
end



K=85; % PSM Size in Cell #
tailbud=15; % # of Tailbud Cells
T=2218; % Simulation Duration
window=T-Ths-200;

PSM=1:K; % Position Vector
time=1:T; % Time Vector
for round=1:2
    cut=0;
    r0=0; % Initial Ligand RNA
    pf0=0; % Initial Free Ligand Protein
    pb0=0; % Initial Bound Ligand Protein
    pr0=0; % Initial Repressor Protein
    pp0=0; % Initial Active Signaling Protein
    p0=0; % Initial Inactive Signaling Protein
    
    RMatrix=r0*ones(K,1)*ones(1,T); % 4D Ligand RNA Matrix
    PFMatrix=pf0*ones(K,1)*ones(1,T); % 4D Ligand Protein Matrix
    PBMatrix=pb0*ones(K,1)*ones(1,T); % 4D Bound Ligand Receptor Protein Matrix
    PRMatrix=pr0*ones(K,1)*ones(1,T); % 4D Repressor Protein Matrix
    PMatrix=p0*ones(K,1)*ones(1,T); % 4D Inactive Signaling Protein Matrix
    PPMatrix=pp0*ones(K,1)*ones(1,T); % 4D Active Signaling Protein Matrix
    StatMatrix=PPMatrix;
    PPMatrixShift=pp0*ones(K,1)*ones(1,T-window);
    
    
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
                
                if ti>actdel&&ti>trlatedel
                    pfy(n+1)=pfy(n)+h*(pfp*RMatrix(k,ti-trlatedel)+pfd*pfy(n)-pb*PFMatrix(k,ti-actdel));
                    pby(n+1)=pby(n)+h*(pb*PFMatrix(k,ti-actdel)+pbd*pby(n));
                else
                    pfy(n+1)=pfy(n);
                    pby(n+1)=pby(n);
                end
                if ti>reprdel
                    pry(n+1)=pry(n)+h*(prp*PPMatrix(k,ti-reprdel)+prd*pry(n));
                    ppy(n+1)=ppy(n)+h*(ppp*pby(n)*py(n)+pd*ppy(n)-ppd*pry(n)*ppy(n));
                    py(n+1)=py(n)+h*(pp+ppd*pry(n)*ppy(n)+pd*py(n)-ppp*pby(n)*py(n));
                    
                else
                    pry(n+1)=pry(n);
                    ppy(n+1)=ppy(n);
                    py(n+1)=py(n);
                    
                end
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
        
        past=floor((k-tailbud)/Vg);
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
                
                if ti>actdel&&ti>trlatedel
                    pfy(n+1)=pfy(n)+h*(pfp*RMatrix(k,ti-trlatedel)+pfd*pfy(n)-pb*PFMatrix(k,ti-actdel));
                    pby(n+1)=pby(n)+h*(pb*PFMatrix(k,ti-actdel)+pbd*pby(n));
                else
                    pfy(n+1)=pfy(n);
                    pby(n+1)=pby(n);
                end
                
                if ti>reprdel
                    pry(n+1)=pry(n)+h*(prp*PPMatrix(k,ti-reprdel)+prd*pry(n));
                    ppy(n+1)=ppy(n)+h*(ppp*pby(n)*py(n)+pd*ppy(n)-ppd*pry(n)*ppy(n));
                    py(n+1)=py(n)+h*(pp+ppd*pry(n)*ppy(n)+pd*py(n)-ppp*pby(n)*py(n));
                    
                else
                    pry(n+1)=pry(n);
                    ppy(n+1)=ppy(n);
                    py(n+1)=py(n);
                end
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
            if well1<k+floor((T-ti)*Vg)&&k+floor((T-ti)*Vg)<well2&&ti>T-Ths
                ppp=ppp_hs;
            else
                ppp=ppp_norm;
            end
            
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
                if ti>trlatedel&&ti>actdel
                    pfy(n+1)=pfy(n)+h*(pfp*RMatrix(k,ti-trlatedel)+pfd*pfy(n)-pb*PFMatrix(k,ti-actdel));
                else
                    pfy(n+1)=pfy(n);
                end
                if ti>actdel
                    pby(n+1)=pby(n)+h*(pb*PFMatrix(k,ti-actdel)+pbd*pby(n));
                else
                    pby(n+1)=pby(n);
                end
                if ti>reprdel
                    pry(n+1)=pry(n)+h*(prp*PPMatrix(k,ti-reprdel)+prd*pry(n));
                    ppy(n+1)=ppy(n)+h*(ppp*pby(n)*py(n)+pd*ppy(n)-ppd*pry(n)*ppy(n));
                    py(n+1)=py(n)+h*(pp+ppd*pry(n)*ppy(n)+pd*py(n)-ppp*pby(n)*py(n));
                    
                else
                    pry(n+1)=pry(n);
                    ppy(n+1)=ppy(n);
                    py(n+1)=py(n);
                    
                end
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
            if well1<k+(T-ti-1)/floor(1/Vg)&&k+(T-ti-1)/floor(1/Vg)<well2&&ti>T-Ths
                ppp=ppp_hs;
            else
                ppp=ppp_norm;
            end
            
            
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
                if ti>trlatedel&&ti>actdel
                    pfy(n+1)=pfy(n)+h*(pfp*RMatrix(k,ti-trlatedel)+pfd*pfy(n)-pb*PFMatrix(k,ti-actdel)+Diff*(PFMatrix(k+1,ti)+PFMatrix(k-1,ti)-2*pfy(n)));
                else
                    pfy(n+1)=pfy(n);
                end
                if ti>actdel
                    pby(n+1)=pby(n)+h*(pb*PFMatrix(k,ti-actdel)+pbd*pby(n));
                else
                    pby(n+1)=pby(n);
                end
                if ti>reprdel
                    pry(n+1)=pry(n)+h*(prp*PPMatrix(k,ti-reprdel)+prd*pry(n));
                    ppy(n+1)=ppy(n)+h*(ppp*pby(n)*py(n)+pd*ppy(n)-ppd*pry(n)*ppy(n));
                    py(n+1)=py(n)+h*(pp+ppd*pry(n)*ppy(n)+pd*py(n)-ppp*pby(n)*py(n));
                    
                else
                    pry(n+1)=pry(n);
                    ppy(n+1)=ppy(n);
                    py(n+1)=py(n);
                    
                end
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
        
        for k=K
            if ti>0
                tf=ti+1; % initial and final time
                h=(tf-ti)/N; % Time step
                ry=zeros(N+1,1);
                pfy=zeros(N+1,1);
                pby=zeros(N+1,1);
                pry=zeros(N+1,1);
                py=zeros(N+1,1);
                ppy=zeros(N+1,1);
                
                ry(1)=RMatrix(k,ti);
                pfy(1)=PFMatrix(k,ti);
                pby(1)=PBMatrix(k,ti);
                pry(1)=PRMatrix(k,ti);
                py(1)=PMatrix(k,ti);
                ppy(1)=PPMatrix(k,ti);
                if ti>reprdel
                    for n=1:N
                        ry(n+1)=ry(n)+h*rd*ry(n);
                        pfy(n+1)=pfy(n)+h*(pfp*RMatrix(k,ti-trlatedel)+pfd*pfy(n)-pb*PFMatrix(k,ti-actdel)+Diff*(PFMatrix(k-1,ti)-2*pfy(n)));
                        pby(n+1)=pby(n)+h*(pb*PFMatrix(k,ti-actdel)+pbd*pby(n));
                        pry(n+1)=pry(n)+h*(prp*PPMatrix(k,ti-reprdel)+prd*pry(n));
                        ppy(n+1)=ppy(n)+h*(ppp*pby(n)*py(n)+pd*ppy(n)-ppd*pry(n)*ppy(n));
                        py(n+1)=py(n)+h*(pp+ppd*pry(n)*ppy(n)+pd*py(n)-ppp*pby(n)*py(n));
                    end
                else
                    ry(N+1)=ry(1);
                    pfy(N+1)=pfy(1);
                    pby(N+1)=pby(1);
                    pry(N+1)=pry(1);
                    py(N+1)=py(1);
                    ppy(N+1)=ppy(1);
                end
                ry0=ry(N+1);
                pfy0=pfy(N+1);
                pby0=pby(N+1);
                pry0=pry(N+1);
                py0=py(N+1);
                ppy0=ppy(N+1);
                RMatrix(k,tf)=ry0;
                PFMatrix(k,tf)=pfy0;
                PBMatrix(k,tf)=pby0;
                PRMatrix(k,tf)=pry0;
                PMatrix(k,tf)=py0;
                PPMatrix(k,tf)=ppy0;
            end
        end
        
        
        
        
        if mod(tf,floor(1/Vg))==0
            if tf>window&&cut==0
                cut=1;
                PPMatrixShift(:,1:tf-window)=PPMatrix(:,window+1:tf);
            end
            
            if tf>window+floor(1/Vg)
                PPMatrixShift(1:K-cut,tf-floor(1/Vg)-window+1:tf-window)=PPMatrix(cut+1:K,tf-floor(1/Vg)+1:tf);
                cut=cut+1;
            end
            StatMatrix(:,tf+1-floor(1/Vg):tf)=PPMatrix(:,tf+1-floor(1/Vg):tf);
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
        end               
        
    end
    
    hs=hstype;
    switch hs
    case 0
        well2=well1;
    case 1
        well2=well1+13;
    case 2
        well2=well1+4;
    end
    
    if round==1
        PPMatrixwell1=PPMatrixShift(well1-ceil((T-window)/floor(1/Vg)),:);
        PPMatrixwell2=PPMatrixShift(well2-ceil((T-window)/floor(1/Vg)),:);
    else
        PPMatrixShift(well1-ceil((T-window)/floor(1/Vg)),:)=PPMatrixwell1;
        PPMatrixShift(well2-ceil((T-window)/floor(1/Vg)),:)=PPMatrixwell2;
    end
    
end
for k=tailbud+1:K-1
    PPMatrixShift(k,:)=PPMatrixShift(k,:)-.02*PPMatrixShift(tailbud+1,:);
end

plot(1:K,PPMatrixShift(:,T-window-1),1:K,PPMatrixShift(:,T-window-72),1:K,PPMatrixShift(:,T-window-143),1:K,PPMatrixShift(:,T-window-214),1:K,PPMatrixShift(:,T-window-285),1:K,PPMatrixShift(:,T-window-356),1:K,PPMatrixShift(:,T-window-427),1:K,PPMatrixShift(:,T-window-498))
    
cond=1;
freq=71;
t0=71;
if cond==4
    frontlist=zeros(10,1);
    front=75;
    k0=5;
    window=window+213;
    k=1;
    deltaF=0;
    deltaF0=PPMatrixShift(front,window+1)-PPMatrixShift(front+k0,window+1);
    while deltaF<deltaF0
        k=k+1;
        deltaF=PPMatrixShift(front-k,72)-PPMatrixShift(front,window+72);
    end
    for i=2:floor((T-window)/71)
        if front>2*k+1
            frontlist(i-1,1)=front;
            front=front-k;
            klist(i-1,1)=k;
            delta=PPMatrixShift(front-k,window+i*71+1)-PPMatrixShift(front,window+i*71+1);
            if delta>deltaF0&&front>k+1
                while delta>deltaF0&&front>k+1
                    k=k-1;
                    delta=PPMatrixShift(front-k,window+i*71+1)-PPMatrixShift(front,window+i*71+1);
                    if delta<=deltaF0
                        k=k+1;
                    end
                end
            else
                while delta<deltaF0&&front>k+1
                    k=k+1;
                    delta=PPMatrixShift(front-k,window+i*71+1)-PPMatrixShift(front,window+i*71+1);
                    if delta>=deltaF0
                        k=k-1;
                    end
                end
            end
            if k<1
                k=1;
            end
        end
    end
else
    
    DiffpLigandMatrix=zeros(K-tailbud,T-window);
    for t=1:T-window
        for k=tailbud+1:K-1
            switch cond
                case 1
                    DiffpLigandMatrix(k-tailbud,t)=PPMatrixShift(k,t); % Constant Threshold
                case 2
                    DiffpLigandMatrix(k-tailbud,t)=sum(PPMatrixShift(k,t-freq:t)); % Time Integration with freq
                case 3
                    DiffpLigandMatrix(k-tailbud,t)=PPMatrixShift(k,t-t0)/PPMatrixShift(k,t)-1; % Time Derivative with t0 step
                case 5
                    DiffpLigandMatrix(k-tailbud,t)=abs(PPMatrixShift(k-1,t)/PPMatrixShift(k,t)-1)-abs(PPMatrixShift(k+1,t)/PPMatrixShift(k,t)-1); % Spatial Fold Change
            end
        end
    end
    
    
    detfront=41;        
    threshold=DiffpLigandMatrix(detfront,1);
    DiffState=zeros(floor((T-window)/71),K-tailbud-1);
    
    for i=1:floor((T-window)/71)
        for j=1:K-tailbud-1
            if DiffpLigandMatrix(j,i*71-70)<=threshold&&max(DiffState(1:i,j))==0
                DiffState(i,j)=i;
            end
        end
    end
                
    DiffState=DiffState-3;
    plot(DiffState(1:floor((T-window)/71),:)','LineStyle','none','Marker','s','MarkerSize',7,'LineWidth',2)
    
    axis([10,70,-2,11])
    ax=gca;
    ax.FontSize=12;
    ax.FontWeight='bold';
    ax.XTick=15:10:65;
    ax.YTick=-2:2:10;
    ax.PlotBoxAspectRatio=[5,1,0.1];
    switch cond
        case 1
            ax.Title.String='Constant Threshold';
        case 2
            ax.Title.String='Integration in Time';
        case 3
            ax.Title.String='Derivative in Time';
        case 4
            ax.Title.String='Derivative in Space';
        case 5
            ax.Title.String='Fold Change in Space';
    end
    
    xlabel('PSM Cell','FontSize',12,'FontWeight','bold')
    ylabel('Cycle','FontSize',12,'FontWeight','bold')
end

