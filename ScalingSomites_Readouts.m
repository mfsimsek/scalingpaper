clear;
count=0;
dataX=[455 415 375 335 305 285 270 259 249];
dataY=[40 30 20 15 11 10 9 8 9];

for pb=0.01:0.01:0.1 % Ligand Protein Receptor Binding Rate
    for ppd=0.015:0.015:0.15 % Signaling Protein Inactivation Rate
        for ppp=0.015:0.015:0.15 % Signaling Protein Activation Rate
            for Diff=0.2:0.2:2 % Ligand Diffusion Coefficient
                count=count+1;
                
                SomiteSize=5; % Initial Somite Size in Cell #
                Vg=0.07; % Tailbud Growth Speed
                
                pbd=-0.5;
                rd=-0.01; % Ligand RNA Decay Rate
                pfp=0.5; % Free Ligand Protein Translation Rate
                pfd=-0.1; % Free Ligand Protein Decay Rate
                
                
                prp=0.1; % Repressor Protein Translation Rate
                prd=-0.1; % Repressor Protein Decay Rate
                
                
                pd=-0.1; % Inactive Signaling Protein Decay Rate
                pp=10; % Inactive Signaling Protein Translation Rate
                rp=5; % Ligand RNA Transcription Rate
                
                trlatedel=14; % Ligand Protein Translation Time Delay
                actdel=14; % Bound Receptor Cascade Activation Time Delay
                reprdel=28; % Repressor Protein Translation Time Delay
                
                
                lock=0;
                
                K=85; % PSM Size in Cell #
                tailbud=15; % # of Tailbud Cells
                T=2218; % Simulation Duration
                PSM=1:K; % Position Vector
                time=1:T; % Time Vector
                
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
                                pfy(n+1)=pfy(n)+h*(pfp*RMatrix(k,ti-trlatedel)+pfd*pfy(n)-pb*PFMatrix(k,ti-actdel)+Diff*(PFMatrix(k+1,ti)-2*pfy(n)));
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
                
                skip=15;
                K=K-skip;
                
                window=200;
                freq=71;
                t0=10;
                
                RMatrix(:,1:window)=RMatrix(:,T-window:T-1);
                PFMatrix(:,1:window)=PFMatrix(:,T-window:T-1);
                PRMatrix(:,1:window)=PRMatrix(:,T-window:T-1);
                PBMatrix(:,1:window)=PBMatrix(:,T-window:T-1);
                PMatrix(:,1:window)=PMatrix(:,T-window:T-1);
                PPMatrix(:,1:window)=PPMatrix(:,T-window:T-1);
                
                Told=T;
                T=1000+window;
                tailbud=tailbud-skip;
                
                RMatrix(1:skip,:)=[];
                PFMatrix(1:skip,:)=[];
                PBMatrix(1:skip,:)=[];
                PRMatrix(1:skip,:)=[];
                PMatrix(1:skip,:)=[];
                PPMatrix(1:skip,:)=[];
                
                RMatrix(:,T+1:Told)=[];
                PFMatrix(:,T+1:Told)=[];
                PBMatrix(:,T+1:Told)=[];
                PRMatrix(:,T+1:Told)=[];
                PMatrix(:,T+1:Told)=[];
                PPMatrix(:,T+1:Told)=[];
                
                for ti=window:T-1
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
                        
                        ry(1)=RMatrix(k,ti);
                        pfy(1)=PFMatrix(k,ti);
                        pby(1)=PBMatrix(k,ti);
                        pry(1)=PRMatrix(k,ti);
                        py(1)=PMatrix(k,ti);
                        ppy(1)=PPMatrix(k,ti);
                        
                        for n=1:N
                            ry(n+1)=ry(n)+h*rp+h*rd*ry(n);
                            pfy(n+1)=pfy(n)+h*(pfp*RMatrix(k,ti-trlatedel)+pfd*pfy(n)-pb*PFMatrix(k,ti-actdel)+Diff*(PFMatrix(k+1,ti)-2*pfy(n)));
                            pby(n+1)=pby(n)+h*(pb*PFMatrix(k,ti-actdel)+pbd*pby(n));
                            pry(n+1)=pry(n)+h*(prp*PPMatrix(k,ti-reprdel)+prd*pry(n));
                            ppy(n+1)=ppy(n)+h*(ppp*pby(n)*py(n)+pd*ppy(n)-ppd*pry(n)*ppy(n));
                            py(n+1)=py(n)+h*(pp+ppd*pry(n)*ppy(n)+pd*py(n)-ppp*pby(n)*py(n));
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
                    
                    for k=2:K-1
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
                        
                        for n=1:N
                            if k<=tailbud
                                ry(n+1)=ry(n)+h*rp+h*rd*ry(n);
                            else
                                ry(n+1)=ry(n)+h*rd*ry(n);
                            end
                            pfy(n+1)=pfy(n)+h*(pfp*RMatrix(k,ti-trlatedel)+pfd*pfy(n)-pb*PFMatrix(k,ti-actdel)+Diff*(PFMatrix(k+1,ti)+PFMatrix(k-1,ti)-2*pfy(n)));
                            
                            pby(n+1)=pby(n)+h*(pb*PFMatrix(k,ti-actdel)+pbd*pby(n));
                            pry(n+1)=pry(n)+h*(prp*PPMatrix(k,ti-reprdel)+prd*pry(n));
                            ppy(n+1)=ppy(n)+h*(ppp*pby(n)*py(n)+pd*ppy(n)-ppd*pry(n)*ppy(n));
                            py(n+1)=py(n)+h*(pp+ppd*pry(n)*ppy(n)+pd*py(n)-ppp*pby(n)*py(n));
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
                    
                    for k=K
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
                        
                        for n=1:N
                            ry(n+1)=ry(n)+h*rd*ry(n);
                            pfy(n+1)=pfy(n)+h*(pfp*RMatrix(k,ti-trlatedel)+pfd*pfy(n)-pb*PFMatrix(k,ti-actdel)+Diff*(PFMatrix(k-1,ti)-2*pfy(n)));
                            pby(n+1)=pby(n)+h*(pb*PFMatrix(k,ti-actdel)+pbd*pby(n));
                            pry(n+1)=pry(n)+h*(prp*PPMatrix(k,ti-reprdel)+prd*pry(n));
                            ppy(n+1)=ppy(n)+h*(ppp*pby(n)*py(n)+pd*ppy(n)-ppd*pry(n)*ppy(n));
                            py(n+1)=py(n)+h*(pp+ppd*pry(n)*ppy(n)+pd*py(n)-ppp*pby(n)*py(n));
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
                
                zmatrix=zeros(K-tailbud,T);
                for i=1:K-tailbud
                    for j=1:T
                        zmatrix(i,j+1)=PPMatrix(tailbud+i,j)-0.02*PPMatrix(tailbud+1,j);
                    end
                end
                xnew = linspace(1,K-tailbud,8*(K-tailbud)+1);
                znew = interp1(zmatrix,xnew,'linear');
                
                
                
                for cond=1:5
                    
                    if cond==4
                        klist=zeros(10,1);
                        frontlist=klist;
                        front=71;
                        front0=(front-skip-tailbud)*8;
                        front=front0;
                        k0=40;
                        k=1;
                        deltaF=0;
                        deltaF0=znew(front,window+1)-znew(front+k0,window+1);
                        %deltaF0=5;
                        while deltaF<deltaF0
                            k=k+1;
                            deltaF=znew(front-k,72)-znew(front,window+72);
                        end
                        for i=2:11
                            if front>2*k+1
                            frontlist(i-1,1)=front-120;
                            front=front-k;
                            klist(i-1,1)=k;
                            %deltaF=znew(front,(i-1)*71+1)-znew(front+k0,(i-1)*71+1);
                            delta=znew(front-k,window+i*71+1)-znew(front,window+i*71+1);
                            if delta>deltaF0&&front>k+1
                                while delta>deltaF0&&front>k+1
                                    k=k-1;
                                    delta=znew(front-k,window+i*71+1)-znew(front,window+i*71+1);
                                    if delta<=deltaF0
                                        k=k+1;
                                    end
                                end
                            else
                                while delta<deltaF0&&front>k+1
                                    k=k+1;
                                    delta=znew(front-k,window+i*71+1)-znew(front,window+i*71+1);
                                    if delta>=deltaF0
                                        k=k-1;
                                    end
                                end
                            end
                            if k<8
                                k=8;
                            end
                            end
                        end
                        frontlist=circshift(frontlist,3);
                        for i=1:3
                            frontlist(4-i)=frontlist(5-i)+k0;
                        end
                    else
                        
                        DiffpMatrix=zeros(560,T-window);
                        for t=window+1:T
                            for k=12:556
                                switch cond
                                    case 1
                                        DiffpMatrix(k,t-window)=sum(znew(k-3:k+4,t)); % Constant Threshold
                                    case 2
                                        DiffpMatrix(k,t-window)=sum(sum(znew(k-3:k+4,t-freq:t))); % Time Integration with freq
                                    case 3
                                        DiffpMatrix(k,t-window)=sum(znew(k-3:k+4,t-t0))/sum(znew(k-3:k+4,t))-1; % Time Derivative with t0 step
                                    case 5
                                        DiffpMatrix(k,t-window)=sum(znew(k-11:k-4,t))/sum(znew(k-3:k+4,t))-1; % Spatial Fold Change
                                end
                            end
                        end
                        
                        frontlist=zeros(12,1);
                        klist=zeros(11,1);
                        frontlist(2)=528;
                        klist(1)=40;
                        threshold=DiffpMatrix(frontlist(2)-120,1);
                        
                        idx =find(abs(DiffpMatrix(12:frontlist(2),72)-threshold)<threshold/500); % Index of coordinate in array
                        if isempty(idx)
                            idx =find(abs(DiffpMatrix(12:frontlist(2),72)-threshold)<threshold/100); % Index of coordinate in array
                            if isempty(idx)
                                idx =find(abs(DiffpMatrix(12:frontlist(2),72)-threshold)<threshold/10); % Index of coordinate in array
                                frontlist(3)=frontlist(2);
                            else
                                frontlist(3)=idx(end);
                            end
                        else
                            frontlist(3)=idx(end);
                        end
                        
                        frontlist(2)=frontlist(3)+klist(1);
                        frontlist(1)=frontlist(2)+klist(1);
                        
                        
                        for i=3:11
                            idx =find(abs(DiffpMatrix(12:frontlist(i),(i-1)*71+1)-threshold)<threshold/500); % Index of coordinate in array
                            if isempty(idx)
                                idx =find(abs(DiffpMatrix(12:frontlist(i),(i-1)*71+1)-threshold)<threshold/100); % Index of coordinate in array
                                if isempty(idx)
                                    idx =find(abs(DiffpMatrix(12:frontlist(i),(i-1)*71+1)-threshold)<threshold/10); % Index of coordinate in array
                                    frontlist(i+1)=frontlist(i);
                                else
                                    frontlist(i+1)=idx(end);
                                end
                            else
                                frontlist(i+1)=idx(end);
                            end
                        end
                        
                        for i=0:9
                            frontlist(i+2)=frontlist(i+1)-klist(i+1);
                            klist(i+2)=frontlist(i+2)-frontlist(i+3);
                            if klist(i+2)<8
                                klist(i+2)=8;
                            end
                        end
                        
                        pds=3;
                        frontlist=circshift(frontlist,pds-1);
                        for i=1:pds-1
                            frontlist(pds-i)=frontlist(pds+1-i)+klist(1);
                        end
                        klist(1)=[];
                        frontlist(1)=[];
                        frontlist(11)=[];
                    end
                    group=zeros(1,16);
                    group(:,10:16)=1;
                    X=[dataX frontlist(1:7).'];
                    Y=[dataY klist(1:7).'];
                    
                    [h,a,c,s]=aoctool(log(X),log(Y),group,0.05,...
                        '','','','off','separate lines');
                    c1=multcompare(s,'Alpha',0.2,'Estimate','slope','Display','off');
                    pslope=c1(1,6);
                    c2= multcompare(s,'Alpha',0.2,'Estimate','intercept','Display','off');
                    pint=c2(1,6);
                    fid1=fopen('parameter.txt','a');
                    fprintf(fid1,'\r\n');
                    fprintf(fid1,'%.3f\t',[cond pb ppd ppp Diff pslope pint]);
                    fclose(fid1);
                    
                    fid2=fopen('frontlist.txt','a');
                    fprintf(fid2,'\r\n');
                    fprintf(fid2,'%d\t',frontlist);
                    fclose(fid2);
                    fid3=fopen('klist.txt','a');
                    fprintf(fid3,'\r\n');
                    fprintf(fid3,'%d\t',klist);
                    fclose(fid3);
                end
                
            end
        end
    end
end
