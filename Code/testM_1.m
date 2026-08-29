clear
close all
clc

tic
M=8:5:53; % Number of elements of the array
SNR_db=-5;
SNR=10.^(SNR_db/10);
No_MntCrl=1;
% qv=[6.0,2.0,0.5,2.5,3.0,1.0,5.5,10.0,6.0, 2.0, 0.5, 2.5, 3.0,1.0,5.5,10.0,6.0,2.0,0.5,2.5,3.0,1.0,5.5,10.0,6.0,2.0,0.5,2.5,3.0,1.0,5.5,10.0,...
%     6.0,2.0,0.5,2.5,3.0,1.0,5.5,10.0,6.0,2.0,0.5,2.5,3.0,1.0,5.5,10.0,6.0,2.0,0.5,2.5,3.0];
for l=1:length(M)
%     Qv=diag(qv(1:M(l)));

    %sigmasq=[.1 .1]; %Vecor of source powers. Power of noise is equal to 1.
    %(Note that these values are not in dB. They are absolute values)
    sumS1=0;
    sumS2=0;
    sumS3=0;
    sumS_IMLSE=0;
    sumS_ILSSE=0;
    %************************************************************************
    %************ Initialization ********************************************
    degsais=[-3 6];
    psais=degsais*pi/180;
    %psais=[10 20 60 -70]*pi/180; %Vecor of source direcions
    
    %sigmasq=[.1 .1]; %Vecor of source powers. Power of noise is equal to 1.
    %(Note that these values are not in dB. They are absolute values)
    scnum=length(psais); %Number of sources
    
    ns=100; %Number of snapshots of array output that are generated and used
    %in estimating spatial correlation matrix R
    c=1500; %Wave velocity
    freq=29900; %Frequency that the beampatterns are obtained for it. Note that
    %this frequency must be very close to fc. In fact
    %(M-1)*abs(freq-fc)/fc must be very smaller than one. Otherwise,
    %the computed beampattern will not be correct.
    fc=30000; % Center frequency of the narrowband signal
    d=c/(2*fc); %Distance between adjacent elements of the array
    cors=0; % The correlation coefficient between the signals of the first two
    % sources will be equal to "cors". Note that cors must be in the
    % interval [-1,1].
    
    %*************************************************************************
    %********************* Data Generation ***********************************
    smt=zeros(M(l),ns);% Matrix containing sum of source signals in the output of array
    ss=zeros(scnum,ns);% Matrix containing source signals in the origin of coordinates
    pmt=zeros(scnum,M(l));% Matrix of steering vectors of sources
    for ii=1:scnum
        tmp1=1i*2*pi*fc*d*(sin(psais(ii)))/c;
        tmp2=1i*2*pi*fc*d*(cos(psais(ii)))/c;
        for jj=1:M(l)
            pmt(ii,jj)= exp((jj-1)*tmp1);%Steering vector of the ii'th source
            pmt_dot(ii,jj)=(tmp2)*(jj-1)*exp((jj-1)*tmp1);
        end
    end
    %% CRB_Stochastic
    A=pmt.'; % Steering Matrix
    for cnt=1:No_MntCrl
    Qv=diag(0.5*randi([1,25],1,M(l)));
    suminvQv=sum(1./diag(Qv));
    sigmasq=(M(l)/suminvQv)*SNR*ones(1,2);
    A_tild=(Qv)^(-0.5)*A;
    D=pmt_dot.'; % A_diff
    D_tild=(Qv)^(-0.5)*D;
    % Projection of A_tild
    proj_A_tild=A_tild*(inv(A_tild'*A_tild))*A_tild';
    ortg_proj=eye(M(l))-proj_A_tild;
    P=diag(sigmasq);
    R=A*P*A'+Qv;
    R_tild=(Qv)^(-0.5)*R*(Qv)^(-0.5);
    M_CRB=2*real(((inv(R_tild)*A_tild*P).').*(D_tild'*ortg_proj));
    T_CRB=inv(conj(inv(R_tild)).*inv(R_tild)-((conj(proj_A_tild*inv(R_tild))).*(proj_A_tild*inv(R_tild))));
    CRB_stoc(:,:,l)=(1/ns)*inv(2*real((P*A_tild'*inv(R_tild)*A_tild*P).*...
                    (D_tild'*ortg_proj*inv(R_tild)*D_tild).')-M_CRB*T_CRB*M_CRB.');
    
    
        smt=zeros(M(l),ns);% Matrix containing sum of source signals in the output of array
        ss=zeros(scnum,ns);% Matrix containing source signals in the origin of coordinates
        if (cors==0)||(scnum==1)
            for ii=1:scnum
                ss(ii,:)=sqrt(sigmasq(ii))*sqrt(1/2)*(randn(1,ns)+1i*randn(1,ns));
                % Note that random('Normal',0,1,1,ns) can be used in place of
                % randn(1,ns)
                smt=smt+(pmt(ii,:)).'*ss(ii,:);
            end
        else
            alpha=sqrt(abs(cors));
            beta=sign(cors)*sqrt(abs(cors));
            v1=sqrt(1/2)*(randn(1,ns)+1i*randn(1,ns));
            v2=sqrt(1/2)*(randn(1,ns)+1i*randn(1,ns));
            v3=sqrt(1/2)*(randn(1,ns)+1i*randn(1,ns));
            ss(1,:)=sqrt(sigmasq(1))*(sqrt(1-alpha^2)*v1+alpha*v3);
            smt=smt+(pmt(1,:)).'*ss(1,:);
            ss(2,:)=sqrt(sigmasq(2))*(sqrt(1-beta^2)*v2+beta*v3);
            smt=smt+(pmt(2,:)).'*ss(2,:);
            if scnum>2
                for ii=3:scnum
                    ss(ii,:)=sqrt(sigmasq(ii))*sqrt(1/2)*(randn(1,ns)+1i*randn(1,ns));
                    smt=smt+(pmt(ii,:)).'*ss(ii,:);
                end
            end
        end
        
        
        scm=(ss*ss')/ns;% Estimating signal covariance matrix
        nsemt=zeros(M(l),ns);% Matrix of noise in the outputs of array elements
        
        nsemt=sqrt(1/2)*sqrt(Qv)*(randn(M(l),ns)+1i*randn(M(l),ns)); % In this case the additive
        % noise of sensors will have complex Gaussian distribution
        % Note that random('Normal',0,1,M,ns) can be used in place of randn(M,ns)
        %nsemt=sqrt(1/2)*(random('unif',-1,1,M,ns)+1i*random('unif',-1,1,M,ns));
        xmt=smt+nsemt;% Adding noise to the outputs of array elements (sensors)
        
        R=(xmt*xmt')/ns;% Estimating output covariance matrix
        
        
        %*************************************************************************
        %*********************    SDP solution  **********************************
        u=diag(ones(M(l)^2,1));
        k=0:M(l)-1;v=(M(l)+1)*k+1;
        u(v,:)=[];
        J=u;
        e = 5;
        cvx_begin quiet
        variable R0(M(l),M(l)) hermitian semidefinite
        minimize( trace(R0 ) )
        subject to
        norm(J*reshape(R0-R,M(l)^2,1),2)<= e
        cvx_end
        
        %*************************************************************************
        %********************* Music Algorithm ***********************************
        %% Traditional Method
        [Vn1,Dn1]=eigs(R,M(l)-scnum,'sr');
        
        %% IMLSE Method
        Q_IMLSE=inv(diag(diag(inv(R)))); % initial value for Nonuniform Covariance Matrix
        epst=0.0001; % Error Threshold
        k=0;
        LL=[];
        while (1)
            R_hat_tild_IMLSE=Q_IMLSE^(-0.5)*R*Q_IMLSE^(-0.5); % Improved Array Covariance Matix
            [Uv_ML,Sv_ML]=eigs(R_hat_tild_IMLSE,scnum,'lr'); % Eigen Value Decomposition
            [Un_ML,Sn_ML]=eigs(R_hat_tild_IMLSE,M(l)-scnum,'sr');
            B_new_IMLSE=Q_IMLSE^(0.5)* Uv_ML*(Sv_ML-eye(scnum,scnum))^(0.5);% B Matrix for DOA Estimation
            k=k+1;
            LL(k)=log10(det(Q_IMLSE))+trace(R_hat_tild_IMLSE)+scnum-sum(diag(Sv_ML)-log10(diag(Sv_ML)));% Likelihood Function
            Q_IMLSE=diag(diag(R-B_new_IMLSE*B_new_IMLSE'));
            if k>1
                if abs(LL(k)-LL(k-1))<=epst
                    break;
                end
            end
            
        end
        %% ILSSE Method
        Q_ILSSE=inv(diag(diag(inv(R)))); % initial value for Nonuniform Covariance Matrix
        epst=0.0001; % Error Threshold
        k=0;
        f=[];
        while (1)
            R_hat_tild_ILSSE=R-Q_ILSSE; % Improved Array Covariance Matix
            [Uv_LS,Sv_LS]=eigs(R_hat_tild_ILSSE,scnum,'lr'); % Eigen Value Decomposition
            [Un_LS,Sn_LS]=eigs(R_hat_tild_ILSSE,M(l)-scnum,'sr'); % Eigen Value Decomposition
            B_new_ILSSE=Uv_LS*Sv_LS^(0.5);% B Matrix for DOA Estimation
            Sigma=Sv_LS'*Sv_LS;
            k=k+1;  
            f(k)=trace(R^2+Q_ILSSE^2-2*Q_ILSSE*R)-trace(Sigma^2);
            Q_ILSSE=diag(diag(R-B_new_ILSSE*B_new_ILSSE'));
            if k>1
                if abs(f(k)-f(k-1))<=epst
                    break;
                end
            end    
        end
        %% Matrix Completion
        [Vn2,Dn2]=eigs(R0,M(l)-scnum,'sr');
        
        %% Proposed Method
        u1=diag(ones(M(l)^2,1));
        k1=0:M(l)-1;v1=(M(l)+1)*k1+1;
        u1(v1,:)=0;
        J1=u1;
        Rv=J1*reshape(R,M(l)^2,1);
        R_Pr=reshape(Rv,M(l),M(l))+min(diag(R))*eye(M(l));
        [Vn3,Dn3]=eigs(R_Pr,M(l)-scnum,'sr');

        %% DOA Estimation
        si=-90:0.01:90;
        cnt1=0; cnt2=0; cnt3=0; cnt4=0; cnt5=0;
        for t=1:length(si)
            tmp1=1i*2*pi*fc*d*(sin(si(t)*pi/180))/c;
            for jj=1:M(l)
                stv(jj,1)= exp((jj-1)*tmp1);%Steering vector of the ii'th source
            end
            a_hat_tild_IMLSE=Q_IMLSE^(-0.5)*stv;
            S1(t)=1/norm(Vn1'*stv,2)^2;
            S2(t)=1/norm(Vn2'*stv,2)^2;
            S3(t)=1/norm(Vn3'*stv,2)^2;
            S_IMLSE(t)=norm(a_hat_tild_IMLSE'*Un_ML,2)^(-2);
            S_ILSSE(t)=norm((stv'*Un_LS),2)^(-2);
            if (t>2 && t<length(si))
                if (10*log10(S1(t-1))>10*log10(S1(t-2))&& 10*log10(S1(t-1))>10*log10(S1(t)) && 10*log10(S1(t-1))>-4)
                    cnt1=cnt1+1;
                    logS1_max(cnt1)=10*log10(S1(t-1));
                    IdxS1(cnt1)=find(10*log10(S1)==logS1_max(cnt1));
                    degS1(cnt1)=-90+(IdxS1(cnt1)/100);
                end
                if (10*log10(S2(t-1))>10*log10(S2(t-2))&& 10*log10(S2(t-1))>10*log10(S2(t)) && 10*log10(S2(t-1))>0)
                    cnt2=cnt2+1;
                    logS2_max(cnt2)=10*log10(S2(t-1));
                    IdxS2(cnt2)=find(10*log10(S2)==logS2_max(cnt2));
                    degS2(cnt2)=-90+(IdxS2(cnt2)/100);
                end
                if (10*log10(S3(t-1))>10*log10(S3(t-2))&& 10*log10(S3(t-1))>10*log10(S3(t)) && 10*log10(S3(t-1))>0)
                    cnt3=cnt3+1;
                    logS3_max(cnt3)=10*log10(S3(t-1));
                    IdxS3(cnt3)=find(10*log10(S3)==logS3_max(cnt3));
                    degS3(cnt3)=-90+(IdxS3(cnt3)/100);
                end
                if (10*log10(S_IMLSE(t-1))>10*log10(S_IMLSE(t-2))&& 10*log10(S_IMLSE(t-1))>10*log10(S_IMLSE(t)) && 10*log10(S_IMLSE(t-1))>0)% 
                    cnt4=cnt4+1;
                    logS_IMLSE_max(cnt4)=10*log10(S_IMLSE(t-1));
                    IdxS_IMLSE(cnt4)=find(10*log10(S_IMLSE)==logS_IMLSE_max(cnt4));
                    degS_IMLSE(cnt4)=-90+(IdxS_IMLSE(cnt4)/100);
                end
                if (10*log10(S_ILSSE(t-1))>10*log10(S_ILSSE(t-2))&& 10*log10(S_ILSSE(t-1))>10*log10(S_ILSSE(t)) && 10*log10(S_ILSSE(t-1))>0)% 
                    cnt5=cnt5+1;
                    logS_ILSSE_max(cnt5)=10*log10(S_ILSSE(t-1));
                    IdxS_ILSSE(cnt5)=find(10*log10(S_ILSSE)==logS_ILSSE_max(cnt5));
                    degS_ILSSE(cnt5)=-90+(IdxS_ILSSE(cnt5)/100);
                end
            end
        end
        plot(si,S1);figure; plot(si,S2);figure ;plot(si,S3);figure; plot(si,S_IMLSE); figure; plot(si,S_ILSSE)
        if (length(degsais)>length(degS1))
            sumS1=(sum((degsais-[degS1 degS1(end)]).^2)/scnum)+sumS1;
        elseif (length(degsais)<length(degS1))
            [logS1_sort,ind_S1]=sort(logS1_max,'descend');
            degS1_sc=-90+(IdxS1(ind_S1(1:2))/100);
            sumS1=(sum((degsais-degS1_sc).^2)/scnum)+sumS1;
        else
            sumS1=(sum((degsais-degS1).^2)/scnum)+sumS1;
        end
        
        if (length(degsais)>length(degS2))
            sumS2=(sum((degsais-[degS2 degS2(end)]).^2)/scnum)+sumS2;
        elseif (length(degsais)<length(degS2))
            [logS2_sort,ind_S2]=sort(logS2_max,'descend');
            degS2_sc=-90+(IdxS2(ind_S2(1:2))/100);
            sumS2=(sum((degsais-degS2_sc).^2)/scnum)+sumS2;
        else
            sumS2=(sum((degsais-degS2).^2)/scnum)+sumS2;
        end
        
        if (length(degsais)>length(degS3))
            sumS3=(sum((degsais-[degS3 degS3(end)]).^2)/scnum)+sumS3;
        elseif (length(degsais)<length(degS3))
            [logS3_sort,ind_S3]=sort(logS3_max,'descend');
            degS3_sc=-90+(IdxS3(ind_S3(1:2))/100);
            sumS3=(sum((degsais-degS3_sc).^2)/scnum)+sumS3;
        else
            sumS3=(sum((degsais-degS3).^2)/scnum)+sumS3;
        end
        
        if (length(degsais)>length(degS_IMLSE))
            sumS_IMLSE=(sum((degsais-[degS_IMLSE degS_IMLSE(end)]).^2)/scnum)+sumS_IMLSE;
        elseif (length(degsais)<length(degS_IMLSE))
            [logS_IMLSE_sort,ind_IMLSE]=sort(logS_IMLSE_max,'descend');
            degS_IMLSE_sc=-90+(IdxS_IMLSE(ind_IMLSE(1:2))/100);
            sumS_IMLSE=(sum((degsais-degS_IMLSE_sc).^2)/scnum)+sumS_IMLSE;
        else
            sumS_IMLSE=(sum((degsais-degS_IMLSE).^2)/scnum)+sumS_IMLSE;
        end
        if (length(degsais)>length(degS_ILSSE))
            sumS_ILSSE=(sum((degsais-[degS_ILSSE degS_ILSSE(end)]).^2)/scnum)+sumS_ILSSE;
        elseif (length(degsais)<length(degS_ILSSE))
            [logS_ILSSE_sort,ind_ILSSE]=sort(logS_ILSSE_max,'descend');
            degS_ILSSE_sc=-90+(IdxS_ILSSE(ind_ILSSE(1:2))/100);
            sumS_ILSSE=(sum((degsais-degS_ILSSE_sc).^2)/scnum)+sumS_ILSSE;        
        else
            sumS_ILSSE=(sum((degsais-degS_ILSSE).^2)/scnum)+sumS_ILSSE;
        end
        cnt
    end
    RMSES1(l)=sqrt(sumS1/No_MntCrl);
    RMSES2(l)=sqrt(sumS2/No_MntCrl);
    RMSES3(l)=sqrt(sumS3/No_MntCrl);
    RMSES_IMLSE(l)=sqrt(sumS_IMLSE/No_MntCrl);
    RMSES_ILSSE(l)=sqrt(sumS_ILSSE/No_MntCrl);
    l
end
for t=1:length(M)
    CRB_stoc1(t)=(CRB_stoc(1,1,t));
end
for t=1:length(M)
    CRB_stoc2(t)=(CRB_stoc(2,2,t));
end
CRB_stocf=sqrt((CRB_stoc1+CRB_stoc2)/scnum)*(180/pi);
semilogy(M,CRB_stocf,'Color',[0 1 0.1],'LineWidth',1.25)
hold on
plot(M,RMSES1,'-.sr')
plot(M,RMSES2,'--*b')
plot(M,RMSES3,'-om')
plot(M,RMSES_IMLSE,'LineStyle','-','Marker','>','Color',[0.9 0.3 0.1],'LineWidth',1.25)
plot(M,RMSES_ILSSE,'LineStyle','-.','Marker','h','Color',[0.1 0.5 0.4],'LineWidth',1.25)
% axis([-5 15 0 10])
title({'RMSEs(deg.^{\bullet}) of DOA estimation versus SNR';'No of Iteration:500,SNR=5dB '...
      ;'Correlation Coefficient=0, Sourse Degree:[-3 6]'})
legend('CRB_{stoc}','Traditional MUSIC','Matrix Completion Method','Proposed Method','IMLSE','ILSSE')
xlabel('Number of Sensors'); ylabel('DOA RMSE(Deg.)');
grid on
toc
