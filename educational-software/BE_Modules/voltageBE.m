function voltageBE(action)
    
    
    % Variable Clamp Voltage taken from the GUI 
    clampVoltage_handle=findobj(gcbf,'Tag','clampVoltage');
    V_clamp1=str2num(get(clampVoltage_handle,'String'));


    EditHandle1a=findobj(gcbf,'Tag','EditText1a');
    use_modified = get(EditHandle1a,'Visible');
    
    % Variable Time delay 
    % is the initial delay before the voltage of the resting stage changes 
    t_delay0_handle=findobj(gcbf,'Tag','t_delay0');
    tdel=str2num(get(t_delay0_handle,'String'));

    % Variable time end
    % The delay time at which the clamp voltage is switched
    t_delay1_handle=findobj(gcbf,'Tag','t_delay1');
    tend=str2num(get(t_delay1_handle,'String'));
    
    if strcmp(use_modified,'on')
        V_clamp0=str2num(get(EditHandle1a,'String'));
        EditHandle2a=findobj(gcbf,'Tag','EditText2a');
        tdel2=str2num(get(EditHandle2a,'String'));
    else 
        tdel2 = tend;
        V_clamp0=V_clamp1;
    end
    
    y_min_handle=findobj(gcbf,'Tag','y_min');
    ymin=str2double(get(y_min_handle,'String'));
    y_max_handle=findobj(gcbf,'Tag','y_max');
    ymax=str2num(get(y_max_handle,'String'));

    na_out_handle=findobj(gcbf,'Tag','na_out'); % LAB_A_Task_3
    Nae=str2num(get(na_out_handle,'String'));

    E_NA=25*log(Nae/45.0)+60;
    eNa_handle=findobj(gcbf,'Tag','eNa');
    e_na_str=num2str(E_NA);
    set(eNa_handle,'String',e_na_str);
    
    PopUpHandle1=findobj(gcbf,'Tag','PopUp1');
    plotstrcell=get(PopUpHandle1,'String');
    plotflag=get(PopUpHandle1,'Value');
    PopUpHandle2=findobj(gcbf,'Tag','PopUp2');
    displaycell=get(PopUpHandle2,'String');
    displayflag=get(PopUpHandle2,'Value');
    
    PopUpHandle3=findobj(gcbf,'Tag','PopUp3');
    ionchannelblockercell=get(PopUpHandle3,'String');
    ionflag=get(PopUpHandle3,'Value');
    
    switch(action)
    case 'run',
    
    %%----------------------------------------------------------------------
    % BME445/545 - Quantitative Neural Function
    %
    % HODGKIN-HUXLEY: 
    % Numerical Simulation of the Response of Squid Axon to Current Injection
    % Voltage Clamping Experiment.
    %----------------------------------------------------------------------
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % INITIALIZATION
    %%%%%%%%%%%%%%%%%%%%%%%
    DT = 0.05;
    T_MAX = tend;
    t = 0:DT:T_MAX;
    
    V_REST = 0.0;
    
    G_NA = 120.0;	% mS/cm^2
    G_K = 25.0;	    % mS/cm^2
    G_LEAK = 0.3;	% mS/cm^2 % Leakage Conductance will not be covered in 445/545.
    
    E_K = -8.0;         % mV
    E_LEAK = 10.613;	% mV  % Leakage Voltage will not be covered in 445/545.
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % MAIN SIMULATION LOOP
    %%%%%%%%%%%%%%%%%%%%%%%
    % This part of code is covered in Lectures 4-6
    clear v  m  h  n			% clear old varibles
    v = V_REST;					% initial membrane voltage
    m = alpha_m(v)/(alpha_m(v)+beta_m(v));	% initial (steady-state) m
    h = alpha_h(v)/(alpha_h(v)+beta_h(v));	% initial (steady-state) h
    n = alpha_n(v)/(alpha_n(v)+beta_n(v));	% initial (steady-state) n
    idel=floor(tdel1/.05);
    idel2=floor(tdel2/.05);
    
    g_K = [];
    g_Na = [];
    
    for i=1:length(t)
        if i<=idel
            V=V_REST;
        elseif i >= idel2
            V=V_clamp1;
        else
            V=V_clamp0;
        end
    
        m_inf = alpha_m(V)/(alpha_m(V)+beta_m(V));
        h_inf = alpha_h(V)/(alpha_h(V)+beta_h(V));
        n_inf = alpha_n(V)/(alpha_n(V)+beta_n(V));
        tau_m = 1/(alpha_m(V)+beta_m(V));
        tau_h = 1/(alpha_h(V)+beta_h(V));
        tau_n = 1/(alpha_n(V)+beta_n(V));
        M     = m_inf-(m_inf-m)*exp(-(i-(idel+1))*DT/tau_m);
        H     = h_inf-(h_inf-h)*exp(-(i-(idel+1))*DT/tau_h);
        N     = n_inf-(n_inf-n)*exp(-(i-(idel+1))*DT/tau_n);
        gNa   = G_NA * M^3 * H;
        gK    = G_K  * N^4;
    
        if i >= idel2
            V=V_clamp1;
        else
            V=V_clamp0;
        end
        
        Im(i) = (gNa*(V-E_NA) + gK*(V-E_K));
        %++BME445/545 *** USEFUL CODE COMMENTED OUT (Used in LAB_B)
        %Im(i) = (gNa*(V-ENA) + gK*(V-EK) + GLEAK*(V-ELEAK));
        %Im(i) = gNa*(V-ENA) + gK*(V-EK)+GLEAK*(V-ELEAK)+(v(i)-v(i-1))/DT;
        %Im(i) = gNa*(V-ENA) + gK*(V-EK)+(v(i)-v(i-1))/DT;
        %--2022.01.24
        INa(i)= gNa*(V-E_NA);
        IK(i) = gK*(V-E_K);
        g_Na(i)= gNa;
        g_K(i)= gK;
    end
    
    
    %% BME 445/545 - LAB_A Voltage Clamp Experiment Plotting
    figure(445)
    
    % This line runs: 'hold on' or 'hold off' based on what is set up in
    % the popup GUI option.
    eval(sprintf('%s',plotstrcell{plotflag}));
    
    % Actual plotting.
    if displayflag == 1
    
        plot(t,Im,'b-');
        axis([0 T_MAX ymin ymax]);
        title('Measured Membrane Current vs Time Plot');
        xlabel('time (mSec)');
        ylabel('Im (mA/sqcm)');
        
    elseif displayflag == 2
        
        plot(t,g_K,'b-');
        axis([0 T_MAX ymin ymax]);
        title('Potassium Conductance vs Time Plot');
        xlabel('time (mSec)');
        ylabel('g_K ');
        
    end
    
    %++BME445/545 LAB_A_CodingTask_I - Programming the Output Metadata - 
    % metadata are parameters associated with this experiment that can be
    % programmatically saved into variable 'meta' as shown below: 
    
    % Task_1b. - Create your metadata variable, 'meta' here - 
    meta = [];
    meta.Im = Im; 
    meta.E_NA = E_NA;
    meta.V_clamp = V_clamp1;
    
    % Task_1c. - Saving your Experiment - 
    % The "eval( --> sprintf( )" below is a powerful way to add flexibile
    % string manipulation to your code and any outputs.   
    %   In the ORIGINAL CODE form, each experiment will overwrite your results.
    %   How will you index different experiments without overwriting (or 
    %   not losing) your repeated voltage clamp experiments? 
    % ORIGINAL CODE:  
    save im Im  % be sure to eventually comment this out.
    
    % ADD MODIFIED CODE HERE :  
    % save ...
    % Syntax Hint: >> eval(sprintf('save im_%s Im meta', ['JP']));
    %--2022.01.24
    
    end
    
    
    % Support functions are provided below:
    %*********************
    function rate = alpha_h(v)
        rate = zeros(size(v));
        rate = 0.07*exp(-v/20.0);
        
    %***********************
    function rate = beta_h(v)
        rate = zeros(size(v));
        rate =  1.0 ./ (exp((-v+30.0)/10.0) + 1.0);
        
    %***********************
    function rate = alpha_m(v)
        rate = zeros(size(v));		% DEFAULT RATE TO ZERO
        rate = 0.1.*(25.-v)./(exp((25.-v)./10)-1);
    
    %**********************
    function rate =  beta_m(v)
        rate = zeros(size(v));
        rate = 4.0*exp(-v/18.0);
    
    %***********************
    function rate = alpha_n(v)
        rate = zeros(size(v));
        rate = 0.01.*(10.-v)./(exp((10.-v)/10)-1);
      
    %**********************
    function rate = beta_n(v)
        rate = zeros(size(v));
        rate = 0.125*exp(-v/80.0);
    