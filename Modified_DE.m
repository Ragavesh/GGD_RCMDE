%finalised for two diagrams 20 06 2020
% this is used for plot 1/f, awgn, noise vs rcmde, mde, de
%this is used to get plot of simulated pd signal vs mde
%imported some portion from mde_signals_noisev1.m
clear all;
close all;

%this is used to get plot of simulated pd signal vs mde
p_noise_sel=[1 1 1 111 10 100];  %[1 1 1 1 10 for pink noise
p_nlevel=[0 0 20 0 100 100]; %[-5 0 5 10

% this is used for plot 1/f, awgn, noise vs rcmde, mde, de
%p_noise_sel=[10 1 111 1 1 1 1 1 101 101 101 101 101 111 111 111 111 111 10 100];  %[1 1 1 1
%p_nlevel=[100 -5 -5 0 5 10 20 -5 0 5 10 20 -5 0 5 10 20 100 100 100]; %[-5 0 5 10

Sig=1;  %
p_method=1;
skip_noise=0;
add_noise=1; %1 means adds noise to PD signals
add_noise_level=10;
show_1st_only=1;
file_print=0;
Num_Iteration = 1;
sig_size=0;

%% parameter settings

signal_size_bits=16;  % 16 used for noise and pd data 12 for emd
max_samples=1024;
signal_size=2^signal_size_bits;
segment_row=256;
segment_col=signal_size/segment_row;
maxpulses=30;

mult=1:0.01:1.27;
a=4.5e5;
b=16e6;
w=2*pi*550e3;
fs=3e3;


col1={'bo','g*','rd','r*','rx','yo','k*','bd','b*','bx','g*','rd','ko','m*','yd','ko','b*','gd'};
col2={'b','g','r','k','m','y','k','b','g','b','g','r','k','m','y','k','b','g'};

testFunctions = {'Doppler','Noise','Doppler','HeaviSine','Blocks','Ramp', 'Cusp', 'Sing', 'HiSine',...
    'LoSine', 'LinChirp', 'TwoChirp', 'QuadChirp', 'MishMash', 'WernerSorrows', ...
    'Leopold', 'Piece-Regular','Riemann','HypChirps','LinChirps', 'Chirps', 'Gabor',...
    'sineoneoverx','Cusp2','SmoothCusp','Piece-Polynomial' };
%signal_size = 2^10; %linspace(0,1,2^10);
Fs = 2*signal_size;
segment_col=signal_size/segment_row;
for p_type=1:size(testFunctions,2)
    
    i_sig=MakeSignal(testFunctions{p_type},signal_size);
    figure(p_type);
    
    plot(i_sig,'LineWidth',1.1)
    xlim([1 2^signal_size_bits]);
    %ylim([min(i_sig)-0.1*max(i_sig) max(i_sig)+0.1*max(i_sig)]);
    %xlabel('Number of Samples');
    title(testFunctions{p_type})
    S1=testFunctions{p_type};
    v_mx=max(i_sig);
    v_mn=min(i_sig);
    Len= size(i_sig,2);
    T = 1/Fs;
    t=(0:Len-1)*T;
    
    %% noise selection and addition
    if skip_noise==1
        N_Sel_run=1;
        Num_Iteration=1;
    else
        N_Sel_run = length(p_noise_sel);
        if(strcmp(S1,'Noise'))
            N_Sel_run=3;
        end
    end
    for N_Sel=1:N_Sel_run
        Awg=' ';
        dsi=' ';
        col=' ';
        if(~skip_noise)
            p_noise = p_noise_sel(N_Sel);
            if p_noise
                if mod(p_noise,10) % check for awgn
                    n_sig=awgn(i_sig,p_nlevel(N_Sel));
                    Awg = p_nlevel(N_Sel);
                    
                    if(Awg== -5)
                        row=1;
                        lam=1.5;
                    elseif (Awg == 0)
                        row=2;
                        lam=0.55;
                    elseif (Awg == 5)
                        row=3;
                        lam=0.55;
                    elseif (Awg == 10)
                        row=4;
                        lam=0.25;
                    else
                        row =5;
                        lam=0.25;
                    end
                else
                    n_sig=i_sig;
                end
                offset=0;
                if (mod(p_noise,100)>9 && mod(p_noise,100)<12) % check for color noise
                    %pN = dsp.ColoredNoise(1,4e3,1);
                    pN = pinknoise(1,size(i_sig,2))';
                    if (max(i_sig)==0)
                        Tn=pN';
                    elseif (max(i_sig)<.2)
                        Tn = pN()'/40;  % it was 40 changed on 30/09 as 20, 20 as 10, tem as 5
                    elseif (max(i_sig)>.2 && (max(i_sig)<2))
                        Tn = pN()'/20;
                    else
                        Tn = pN()'/10;
                    end
                    
                    e_len=length(n_sig)-length(Tn);
                    if(e_len>0)
                        Tn = wextend(1,'sym', Tn, e_len,'r');
                    else
                        Tn= Tn(1:length(n_sig));
                    end
                    n_sig=Tn+n_sig;
                    col=' COL';
                    offset=4;
                end
                if (mod(p_noise,1000)>99 && mod(p_noise,1000)<112) % check for AM noise
                    m=0.4;
                    fa=1000; % Frequency of modulating signal
                    Ta=1/fa;
                    fc=0.4e6:0.1e6:0.8e6; % Frequency of carrier signal
                    t=0:1e-6:3e-3;
                    am=0;
                    if (max(i_sig)<.2)
                        am_noise_level=0.002;
                    elseif (max(i_sig)>.2 && (max(i_sig)<2))
                        am_noise_level=0.02;
                    else
                        am_noise_level=0.02;
                    end
                    
                    for i=1:5
                        am=am+(1+m*sin(2*pi*fa*t)).*cos(2*pi*fc(i)*t);
                    end
                    Tn = am*am_noise_level;
                    e_len=length(n_sig)-length(Tn);
                    if(e_len>0)
                        Tn = wextend(1,'sym', Tn, e_len,'r');
                    else
                        Tn= Tn(1:length(n_sig));
                    end
                    n_sig=Tn+n_sig;
                    dsi=' DSI';
                    offset=4;
                end
            end
            if(~strcmp(Awg,' ')&& strcmp(dsi, ' DSI')&& strcmp(col, ' COL'))
                S4=sprintf(' AWGN %ddB, %s, %s',Awg, dsi,'1/f');
            elseif(strcmp(Awg,' ')&& strcmp(dsi, ' ')&& strcmp(col, ' COL'))
                S4=sprintf(' %s', '1/f noise');
                %             elseif(strcmp(Awg,' ')&& strcmp(dsi, ' ')&& strcmp(col, ' '))
                %                 S4=sprintf(' AWGN %ddB, %s',Awg);
            else
                S4=sprintf(' AWGN %ddB',Awg);
            end
        else
            if(~add_noise)
                n_sig=i_sig;
            end
            % n_sig=i_sig;
            lam=0.15;
            S4=sprintf('PD %s', S1);
        end
        %n_sig = -1 + 2.*(n_sig - min(n_sig))./(max(n_sig) - min(n_sig));
        if(show_1st_only)
            %print the created signals
            figure(1);
            clf;
            subplot(1,2,1)
            plot(i_sig); xlim([1 Len]); ylabel('Amplitude');
            xlabel('Sample points');
            subplot(1,2,2)
            plot(n_sig); xlim([1 Len]); ylabel('Amplitude');
            xlabel('Sample points');
            set(gcf, 'Position',  [800, 500, 500, 150])
        end
%% select PD data
if(N_Sel==1)
        n_sig=importdata('E:/Dropbox/PhD_GCU/PSA_V1/PD_Data1.txt'); %PD100m500ns
        n_sig=n_sig(1:signal_size)';
else
        n_sig=importdata('E:/Dropbox/PhD_GCU/PSA_V1/PD100m500ns.txt'); %PD100m500ns
        n_sig=n_sig(1:signal_size)';
end
 % to plot rcmde for wgn, pink, a block noisy signal
 % select noise/signal type from the signal, 10db awgn, pink, 2^12 selected
 
plt_val=[];   % first row mde, then rcmde and mod_rcmde b=1, 2,3
plt_val=[plt_val; MDE(n_sig,3,6,1,15)];
plt_val=[plt_val; RCMDE(n_sig,3,6,1,15)]; 
for beta=1:3
      plt_val=[plt_val; mod_RCMDE(n_sig,3,6,1,15,beta)];
end
figure(2);
set(gcf, 'Position',  [1100 330 720 580])
plt_noise=["AWGN","Pink"];
%plt_noise=["PD Data","Field Noise"];
plt_label=["MDE","NCDF",strcat('\beta=', num2str(1)), strcat('\beta=', num2str(2)), strcat('\beta=', num2str(3)) ];
for n=1:5
    plot(1:15, plt_val(n,:),strcat(col1{(N_Sel-1)*5+n},'-'),'LineWidth',1.15,'MarkerSize',5,'DisplayName',strcat(plt_noise{N_Sel},"-",plt_label{n}));
    hold on;
end
    xlabel('Scale factor');
    ylabel('Entropy Value');
    h =legend('show','location','best');
    set(h,'FontSize',8);
    set(h, 'NumColumns',2)
   
    x=5;
    
%end of "to plot rcmde for wgn, pink, a block noisy signal

%% old code only for plots
 %% reshape and segment the signal and find de and rcmde
        re_n_sig=reshape(n_sig, segment_row, segment_col);
        DE=[];
        DE1=[];
        o_MDE=[];
        o_RCMDE=[];
        o_MSE=[];
        for num_it=1:segment_col
            DE=[DE DisEn_NCDF(re_n_sig(:,num_it)',2,3,1)]; %DisEn_NCDF(x,m,nc,tau)
            DE1=[DE1 DisEn_GCDF(re_n_sig(:,num_it)',2,3,1,1)]; %DisEn_NCDF(x,m,nc,tau)
            tmp=MDE(re_n_sig(:,num_it),2,3,1,12);
            o_MDE=[o_MDE tmp(1,numel(tmp))];  %MDE(x,m,c,tau,Scale)
            tmp=RCMDE(re_n_sig(:,num_it),2,3,1,12);
            o_RCMDE=[o_RCMDE tmp(1,numel(tmp))];  %MDE(x,m,c,tau,Scale)
            o_MSE= [o_MSE multiscaleSampleEntropy(re_n_sig(:,num_it), 2, 0.15, 12)]; %signal, m, r(0.15??), tau (scale)
        end
        
        %% end of emd
        figure(10);
        clf;
        subplot(2,2,1);
        plot(DE,o_MSE, 'k+','DisplayName','MSE','MarkerSize',4);
        hold on;
        plot(DE,o_RCMDE, 'b*','DisplayName','RCMDE','MarkerSize',4)
        %plot(DE,o_MDE, 'ro')
        xlabel('DE');
        ylabel('Entropy value');
        legend('show','location','northeast');
        xlim([0 4]);
        ylim([0 4]);
        
        %figure(11);
        %clf;
        subplot(2,2,2);
        plot(o_MSE,DE, 'k+','DisplayName','DE','MarkerSize',4);
        hold on;
        plot(o_MSE,o_RCMDE, 'b*','DisplayName','RCMDE','MarkerSize',4)
        %plot(o_MSE,o_MDE, 'ro')
        xlabel('MSE');
        ylabel('Entropy value');
        legend('show','location','northeast');
        xlim([0 4]);
        ylim([0 4]);
        
        %figure(12);
        %clf;
        subplot(2,2,3);
        plot(o_RCMDE,DE, 'b*','DisplayName','DE','MarkerSize',4)
        hold on;
        %plot(o_RCMDE,o_MDE, 'k+','DisplayName','MDE','MarkerSize',4)
        plot(o_RCMDE,o_MSE, 'ro','DisplayName','MSE','MarkerSize',4)
        xlabel('RCMDE');
        ylabel('Entropy value');
        legend('show','location','northeast');
        xlim([0 4]);
        ylim([0 4]);
        
        %figure(13);
        %clf;
        subplot(2,2,4);
        plot(o_MDE,DE, 'b*','DisplayName','DE','MarkerSize',4)
        hold on;
        %plot(o_MDE,o_MSE, 'k+','DisplayName','MSE'); %
        plot(o_MDE,o_RCMDE, 'ro','DisplayName','RCMDE','MarkerSize',4)
        xlabel('MDE');
        ylabel('Entropy value');
        legend('show','location','northeast');
        xlim([0 4]);
        ylim([0 4]);
        
        %nt=sprintf('Signal:  %s',S4);
        %title(nt);
        figure(34);
        sig_size=sig_size+1;
        xax=1:1:segment_col;
        name=S4;
        
        plot(xax,DE,strcat(col1{sig_size}),'MarkerSize',4,'DisplayName',name);
        xlabel('Sample Segment');
        ylabel('DE');
        h=legend('Location','best','Orientation','horizontal');
        legend('show','location','best');
        set(h,'FontSize',10);
        set(h, 'NumColumns',2)
        ylim([1 3.7]);
        
        hold on;
        figure(35);
        plot(o_RCMDE,DE,strcat(col1{sig_size}),'MarkerSize',4,'DisplayName',name);
        hold on;
        xlabel('RCMDE');
        ylabel('DE');
        h=legend('Location','best','Orientation','horizontal');
        %legend('show','location','best');
        set(h,'FontSize',10);
        set(h, 'NumColumns',2)
        xlim([0.25 3.5]);
        ylim([1 3.6]);
        
        %% patch DE
        figure(5);
        plot(DE);
        DE=[DE DisEn_NCDF(re_n_sig(:,num_it)',2,5,1)]; %DisEn_NCDF(x,m,nc,tau)
        n1=1;
        for beta=1.0:0.5:3
            DE1=[];
            for num_it=1:segment_col
                DE1=[DE1 DisEn_GCDF(re_n_sig(:,num_it)',2,5,1,beta)]; %DisEn_NCDF(x,m,nc,tau)
            end
            hold on;
            plot(DE1,strcat(col2{n1}));
            n1=n1+1;
        end
        %% apply emd and test each imf
        
        IMF = emd(n_sig);
figure(4);
for beta=1:3
            x=IMF(beta,:);
          %  beta=2;
        s=std(x);
m_x=mean(x);
    rho=sqrt(gamma(1/beta)/gamma(3/beta))*s;
    % %pdf
    Fx=(beta/(2*rho*gamma(1/beta)))*exp(-(abs(x - m_x)/rho).^beta);
    % %cdf
    %c = 0.5*sign(x-m_x);
    %Fx = (c+0.5) - c.*(gammainc((1.0/beta),(abs(x-m_x)/rho)).^beta); %/(2*gamma(1/beta));
    
    plot(x,Fx,'.')
    hold on;
    pd=makedist('Normal',m_x,s); % distribution, mu and sigma value.
    Fx=pdf(pd,x);
    plot(x,Fx,'.')
end
        decom_size=size(IMF,1);
        limit1=size(IMF,1);
        if(limit1>12)  % this is for plotting first 12 imfs only
            limit1=12;
        end
        if(mod(limit1,2)~=0)
            limit1=limit1-1; % ignore last residue to make plot uniform ??
        end
        if(limit1<=6)
            lx=3;
            ly=2;
        elseif (limit1>6 && limit1<=12)
            lx= 3;
            ly= 4;
        else
            lx=4;
            ly=4;
        end
        figure(7);
        clf;  % create only one plot
        for n=1:limit1
            re_n_sig=reshape(IMF(n,:), segment_row, segment_col);
            m_rcdme=[];
            for num_it=1:segment_col
                tmp= RCMDE(re_n_sig(:,num_it),2,6,1,12); %DisEn_NCDF(x,m,nc,tau) ;
                m_rcdme = [m_rcdme tmp(1,numel(tmp))];
            end
            subplot_tight(lx,ly,n);
%            subplot(lx,ly,n)
            plot(m_rcdme,'LineWidth',1.25,'DisplayName','NCDF');
            hold on;
            ylim([1.25 3.5]);
            n1=1;
            for beta=1.0:1:3
                m_rcdme=[];
                for num_it=1:segment_col
                    tmp= mod_RCMDE(re_n_sig(:,num_it),2,6,1,12,beta); %DisEn_NCDF(x,m,nc,tau) ;
                    m_rcdme = [m_rcdme tmp(1,numel(tmp))];
                end
                hold on;
                plot(m_rcdme,strcat(col2{n1},'-'),'LineWidth',1.25,'DisplayName',strcat('\beta=', num2str(beta)));
                n1=n1+1;
            end
            h =legend('show','location','best');
            set(h,'FontSize',8);
            set(h, 'NumColumns',1)
            % set(gcf, 'Position',  [1400 250 520 760])
            if(n > (lx-1)*ly)
            xlabel('Sample Segment');
            end
            if(mod(n,ly)==1)
            ylabel('Entropy Value');
            end
            Lgnd = legend('show');
            Lgnd.Position(1) = 0.84;
            Lgnd.Position(2) = 0.2;             
        end
        
    end
end
