function voltageBEG(action)


     %% ----------------------------
     % 1. Clamp Voltage Variables
     %% ----------------------------

     % Main voltage clamp value (V_clamp1)
     clampVoltage_handle = findobj(gcbf,'Tag','clampVoltage');
     V_clamp1 = str2num(get(clampVoltage_handle,'String'));

     %% ----------------------------
     % 2. Plot Range Settings
     %% ----------------------------

     % Minimum Y-axis value for plots
     y_min_handle = findobj(gcbf,'Tag','y_min');
     ymin = str2double(get(y_min_handle,'String'));

     % Maximum Y-axis value for plots
     y_max_handle = findobj(gcbf,'Tag','y_max');
     ymax = str2num(get(y_max_handle,'String'));

     %% ----------------------------
     % 3. Potential Variables
     %% ----------------------------

     % External sodium concentration
     na_out_handle = findobj(gcbf,'Tag','na_out');
     Nae = str2num(get(na_out_handle,'String'));

     E_NA = 25 * log(Nae / 45.0) + 60;
     eNa_handle = findobj(gcbf,'Tag','eNa');
     e_na_str = num2str(E_NA);
     set(eNa_handle, 'String', e_na_str);

     %% ----------------------------
     % 4. GUI Options (Display & Plot Flags)
     %% ----------------------------

     PopUpHandle1 = findobj(gcbf,'Tag','PopUp1');
     plotstrcell = get(PopUpHandle1,'String');
     plotflag = get(PopUpHandle1,'Value');


     PopUpHandle2 = findobj(gcbf,'Tag','PopUp2');
     displaycell = get(PopUpHandle2,'String');
     displayflag = get(PopUpHandle2,'Value');

     PopUpHandle3 = findobj(gcbf,'Tag','PopUp3');
     ionchannelblockercell = get(PopUpHandle3,'String');
     ionflag = get(PopUpHandle3,'Value'); 

     %% ----------------------------
     % 5. Time Variables
     %% ----------------------------

     % Initial delay before voltage clamp starts (t_delay0)
     t_delay0_handle = findobj(gcbf,'Tag','t_delay0');
     tdel1 = str2num(get(t_delay0_handle,'String'));

     % Final delay / switch time (t_delay1)
     t_delay1_handle = findobj(gcbf,'Tag','t_delay1');
     tend = str2num(get(t_delay1_handle,'String'));

     % Optional delay and clamp value if modification is visible
     EditHandle1a = findobj(gcbf,'Tag','EditText1a');
     use_modified_str = get(EditHandle1a,'Visible');
     use_modified_flag = strcmp(use_modified_str, 'on');

     V_clamp0 = V_clamp1; % Default if not using modified
     tdel2 = tend;       % Default if not using modified

     if use_modified_flag
         V_clamp0 = str2num(get(EditHandle1a,'String'));
         EditHandle2a = findobj(gcbf,'Tag','EditText2a');
         tdel2 = str2num(get(EditHandle2a,'String'));
     end
     %% ----------------------------

     switch(action)
     case 'run'

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


     % Constants used by the CUDA kernel (hardcoded in .cu file or can be passed)
     V_REST = 0.0;
     G_NA = 120.0;	
     G_K = 25.0;	    
     E_K = -8.0;         

     
    %%%%%%%%%%%%%%%%%%%%%%%
    % MAIN SIMULATION LOOP
    %%%%%%%%%%%%%%%%%%%%%%%
    
    clear v  m  h  n			% clear old varibles
    v = V_REST;					% initial membrane voltage
    m = alpha_m(v)/(alpha_m(v)+beta_m(v));	% initial (steady-state) m
    h = alpha_h(v)/(alpha_h(v)+beta_h(v));	% initial (steady-state) h
    n = alpha_n(v)/(alpha_n(v)+beta_n(v));	% initial (steady-state) n
    idel=floor(tdel1/.05);
    idel2=floor(tdel2/.05);
    
    g_K = [];
    g_Na = [];
    tic;

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

        INa(i)= gNa*(V-E_NA);
        IK(i) = gK*(V-E_K);
        g_Na(i)= gNa;
        g_K(i)= gK;
    end
    noncuda_time = toc;

    fprintf('Simulation Time using NON-CUDA %.4f \n', noncuda_time);

     %%%%%%%%%%%%%%%%%%%%%%%
     % MAIN SIMULATION - USING CUDA MEX
     %%%%%%%%%%%%%%%%%%%%%%%

     % The original for loop and associated initializations are replaced by a call to the MEX function.
     % The MEX function performs the time-stepping calculations on the GPU.

     % Call the CUDA MEX function
     % Inputs: V_clamp1, Nae, tdel1, tend, use_modified_flag (0 or 1), [V_clamp0, tdel2 if use_modified_flag is 1]
     % Outputs: Im, INa, IK, g_Na, g_K, E_NA
 
     
     tic;
         if use_modified_flag
              [Im, INa, IK, g_Na, g_K, E_NA_mex] = voltage_kernel(V_clamp1, Nae, tdel1, tend, 1, V_clamp0, tdel2);
         else
              [Im, INa, IK, g_Na, g_K, E_NA_mex] = voltage_kernel(V_clamp1, Nae, tdel1, tend, 0);
         end

         % Update the ENa GUI handle string with the E_NA value returned by the MEX
         eNa_handle = findobj(gcbf,'Tag','eNa');
         set(eNa_handle, 'String', num2str(E_NA_mex));

     cuda_time = toc;
     fprintf('Simulation Time using CUDA %.4f \n', cuda_time);




     %% BME 445/545 - LAB_A Voltage Clamp Experiment Plotting
     figure(445)

     % This line runs: 'hold on' or 'hold off' based on what is set up in
     % the popup GUI option.
     eval(sprintf('%s',plotstrcell{plotflag}));

     % Actual plotting.
     % Ensure the output arrays from the MEX function are not empty before plotting
     if ~isempty(t) && ~isempty(Im) && ~isempty(g_K)
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
     else
         % Handle case where simulation failed and arrays are empty
         title('Simulation Failed');
         xlabel('');
         ylabel('');
         cla; % Clear the current axes
     end


     end % End of case 'run'
     % Add other cases here if action can be something else (e.g., 'stop', 'reset')

 end % End of function