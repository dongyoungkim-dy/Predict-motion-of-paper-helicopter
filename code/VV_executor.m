function [area_metric, tha] = 
VV_executor(model,h,clip_no,time_data,calib_data_cd,cd_model,simul_kind)
%%% Function to execute MCMC simulations and predictive distribution area
%%% metric comparison
 
Heli = struct('model',model,'calib_data',calib_data_cd,'time_data',time_data,...
    'clip_no',clip_no,'height',h,'cd_model',cd_model,'simul',simul_kind);
% Store the input data into struct variable named 'Heli'
 
simul_kind = Heli.simul;
% Select Validation area metric or P-box approach 
 
n= length(Heli.calib_data); mx = mean(Heli.calib_data); sx = std(Heli.calib_data); N=1e4*2; w=.005;
tha = MCMC(n,mx,sx,N,w);
%%% Executing MCMC simulations with the drag coefficient data of calibration conditions
 
    switch simul_kind
    % Depending on simul_kind
        case 'vali'
        % Run Validation area metric
            area_metric = validationAreaMetric(N,tha, Heli);
        case 'pbox'
        % Run P-box approach (not perfect, reliable)
            M=400;
            area_metric = pbox_simulation(N,M,tha,Heli);   
    end
end
function tha = MCMC(n,mx,sx,N,w)
%%% The reference code is from the professor Joo-ho choi,
%%% from School of Aerospace and Mechanical Engineering, Korea Aerospace University
 
    % joint posterior distribution
    s2=sx^2; % my example
    dx=0.005; dy=0.002; xx=mx-0.2:dx:mx+0.2;  yy=dy:dy:0.1;
    [x,y]=meshgrid(xx,yy); % x is mew, y is sig2.
    pdf=y.^(-n/2-1).*exp(-1/2*y.^(-1).*((n-1)*s2+n*(mx-x).^2));
    figure; contour(x,y,pdf); grid on figure; mesh(x,y,pdf);
 
    % MH algorithm for two parameters (Gelman p290)   
    th=[mx sx^2]; tha=zeros(N,2); tha(1,:)=th; % initialize th. 
    for i=2:N;
        ths=th+w*(2*rand(1,2)-1); % sample ths from proposal N(th,0.2I)
        pths= ths(2).^(-n/2-1).*exp(-1/2*ths(2).^(-1).*((n-1)*s2+n*(mx-ths(1)).^2));
        pth= th(2).^(-n/2-1).*exp(-1/2*th(2).^(-1).*((n-1)*s2+n*(mx-th(1)).^2));
        if ths(2) >0 & rand<pths/pth; th=ths; end
        tha(i,:)=th;
    end 
    tha(:,2) = sqrt(tha(:,2));
    figure(2); plot(tha(:,1),tha(:,2),'.'); figure;
    hist3([tha(:,1),tha(:,2)],[60,60]); % 3D histogram. 
    close all
    
end
function area_metrics = validationAreaMetric(N,tha, Heli)
    %%% Generate data randomly from each N(mx,sx)
    model = Heli.model; % Current Heli model (S/M/L/H)
    exper_t = Heli.time_data; % Current time data for validation
    clip_no = Heli.clip_no;
    h = Heli.height;
    cd_model = Heli.cd_model;
    
    %%% With N groups of mean and std of drag coefficient, generate one sample cd from each group
    M = 1; 
    extracted_data = zeros(N,M);   
    for i=1:1:N
        extracted_data(i,:) = normrnd(tha(i,1),tha(i,2),[M 1]); 
        % Generate single sample Cd 
        while extracted_data(i,:) < 0 
        %%% Case when negative value returned
            extracted_data(i,:) = normrnd(tha(i,1),tha(i,2),[M 1]); 
            % Generate single sample Cd 
        end
    end
    
    %%% Convert N Cd of calibration conditions to fall time with the selected other helicopter condition
    [rho, g, w, a] = modelSelector(model, clip_no);
    falltime = zeros(N,1);
    for j=1:1:N
        falltime(j,:) = fallTimeCalculator(cd_model,w,rho,extracted_data(j,:),a,h,g);
        % Calculate the fall time of selected heli. model with the calibrated Cd
    end
    
    fig = figure; 
    % Genereate a figure
    mx = mean(falltime); sx = std(falltime); dx = 0.0001;
    % Calculate the mean, std and decide the differential
    
    % Calculate the min and max value for selecting the range of the CDF
    if mx-1 < 0
        lx = min(mx-0.5,min(exper_t));
        ux = max(mx+0.5,max(exper_t));
    else     
        lx = min(mx-1,min(exper_t));
        ux = max(mx+1,max(exper_t));
    end

    xx = transpose(lx:dx:ux); 
    yy = normcdf(xx,mx,sx);
    p_data = [xx yy];
    
    h = plot(xx,yy,'m'); hold on grid on set(h,'Color','r');
    % Plotting the predictive CDF 
    h1 = cdfplot(exper_t); set(h1,'Color','k');
    % Plotting the empirical CDF
        em_data= [transpose(get(h1,'XData')), transpose(get(h1,'YData'))];
    % Get the data of cdfplot
    
    %%% Calculating validation area metrics
    ff = 1; ranges = zeros(1,2);
    for kk=2:2:length(em_data)
        indices = find(em_data(kk-1,1) < p_data(:,1) & p_data(:,1) < em_data(kk,1));
        if ~isempty(indices)
            p_data(em_data(kk-1,1) < p_data(:,1) & p_data(:,1) < em_data(kk,1),3) = em_data(kk,2);
            ranges(ff,:) = [indices(1) indices(end)];
            ff = ff+1;    
        end
    end
          
    temp = abs(p_data(:,2)-p_data(:,3));
    area_metrics = 0;
    for kk=2:1:length(temp)-2
        dx = p_data(kk+1,1)-p_data(kk,1);
        area_metrics = area_metrics+dx*(temp(kk+1)+temp(kk))/2;
    end
    
    %%% Fill the area metric with green color
    for i=1:1:length(ranges)
        index1 = ranges(i,1); index2 = ranges(i,2);
        index = [index1 index1:index2 index2];
        
        fill(p_data(index,1),...
            [p_data(index1,3); p_data(index1:index2,2);
 p_data(index1,3)],'g','EdgeColor','none');
        drawnow;
    end
    %%% Redrawing the plot
    delete([h h1]);
    h = plot(xx,yy,'m'); hold on grid on set(h,'Color','r'); 
    h1 = cdfplot(exper_t); set(h1,'Color','k');  
    xlabel('Fall time (sec)'); ylabel('CDF of fall time');
end
 
function area_metrics = pbox_simulation(N,M,tha, Heli)
    %%% Generate data randomly from each N(mx,sx)
    model = Heli.model; % Current Heli model (S/M/L/H)
    exper_t = Heli.time_data;
    clip_no = Heli.clip_no;
    h = Heli.height;
    cd_model = Heli.cd_model;
    
    %%% With N groups of mean and std of drag coefficient, generate M sample cd from each group
    extracted_data = zeros(N,M);
    pbox_cdf = zeros(N,2); 
    
    for i=1:1:N
        extracted_data(i,:) = normrnd(tha(i,1),tha(i,2),[M 1]);
        % Generate M sample Cd 
    end
       [rho, g, w, a] = modelSelector(model, clip_no);
    parfor j=1:N
        falltime(j,:) = fallTimeCalculator(cd_model,w,rho,extracted_data(j,:),a,h,g);
    end
    
    % With 2.5 and 97.5 percentile index of M, plotting CDF 
    for i=1:N
        temp = sort(falltime(i,:));
        pbox_cdf(i,:) = [temp(M*0.025) temp(M*0.975)];
    end
    
    figure; 
    [xx1,cdf_data1,h] = cdf_plot(pbox_cdf(:,1),1.5);
    hold on
    [xx2,cdf_data2,h1] = cdf_plot(pbox_cdf(:,2),1.5);
    h2 = cdfplot(exper_t);
    
    %%% Lower bound
    p_data = [xx1' cdf_data1'];
    em_data= [transpose(get(h2,'XData')), transpose(get(h2,'YData'))];  
    
    %%% Calculating validation area metrics (not perfect)
    ff = 1; ranges = zeros(1,2);
    for kk=2:2:length(em_data)
        indices = find(em_data(kk-1,1) < p_data(:,1) & p_data(:,1) < em_data(kk,1));
        if ~isempty(indices)
            p_data(em_data(kk-1,1) < p_data(:,1) & p_data(:,1) < em_data(kk,1),3) = em_data(kk,2);
            ranges(ff,:) = [indices(1) indices(end)];
            ff = ff+1;    
        end
    end
          
    temp = p_data(:,3)-p_data(:,2);
    index = find(temp>0);
    if ~isempty(index)
        if index(end) == length(p_data)
            index(end) = [];
        end
 
        area_metrics1 = 0;
        for kk=1:1:length(index)
            dx = abs(p_data(index(kk)+1,1)-p_data(index(kk),1));
            area_metrics1 = area_metrics1+dx*abs((temp(index(kk)+1)+temp(index(kk))))/2;
 
        end

        %%% Lower
        for i=1:1:length(ranges)
            index1 = ranges(i,1); index2 = ranges(i,2);
            temp = p_data(index1:index2,3)-p_data(index1:index2,2);
            parseindex = find(temp>0);
            if ~isempty(parseindex)
                indices = index1:index2;
                indices = indices(parseindex); %#ok<FNDSB>
                index = [indices(1) indices indices(end)]; 
                h = fill(p_data(index,1),...
                    [p_data(indices(1),3); p_data(indices,2); p_data(indices(1),3)],'g','EdgeColor','none');
 
                drawnow;
            end
        end
    else
        area_metrics1 = 0;
    end
    %%% Upper bound
    p_data = [xx2' cdf_data2'];
    em_data= [transpose(get(h2,'XData')), transpose(get(h2,'YData'))];  
    
    ff = 1; ranges = zeros(1,2);
    for kk=2:2:length(em_data)
        indices = find(em_data(kk-1,1) < p_data(:,1) & p_data(:,1) < em_data(kk,1));
        if ~isempty(indices)
            p_data(em_data(kk-1,1) < p_data(:,1) & p_data(:,1) < em_data(kk,1),3) = em_data(kk,2);
            ranges(ff,:) = [indices(1) indices(end)];
            ff = ff+1;    
        end
    end
          
    temp = p_data(:,3)-p_data(:,2);
    index = find(temp<0);
    if ~isempty(index)
        if index(end) == length(p_data)
            index(end) = [];
        end
        area_metrics2 = 0;
        for kk=1:1:length(index)
            dx = abs(p_data(index(kk)+1,1)-p_data(index(kk),1));
            area_metrics2 = area_metrics2+dx*abs((temp(index(kk)+1)+temp(index(kk))))/2;
        end
 
        %%% Upper
        for i=1:1:length(ranges)
            index1 = ranges(i,1); index2 = ranges(i,2);
            temp = p_data(index1:index2,3)-p_data(index1:index2,2);
            parseindex = find(temp<0);
            if ~isempty(parseindex)
                indices = index1:index2;
                indices = indices(parseindex); %#ok<FNDSB>
                index = [indices(1) indices indices(end)]; 
                fill(p_data(index,1),...
                    [p_data(indices(1),3); p_data(indices,2);
 p_data(indices(1),3)],'g','EdgeColor','none');
 
                drawnow;
            end
        end
    else
        area_metrics2 = 0;
    end
    
    %%% Redrawing
    delete([h h1 h2]);
    [~,~,h] = cdf_plot(pbox_cdf(:,1),0.5);
    hold on grid on
    [~,~,h1] = cdf_plot(pbox_cdf(:,2),0.5);
    h2 = cdfplot(exper_t);
    set(h,'Color','r');
    set(h1,'Color','b');
    set(h2,'Color','k');
    xlabel('Fall time (sec)'); ylabel('CDF of fall time');
    
    area_metrics = area_metrics1+area_metrics2;
    % Calculate the total area metric
end