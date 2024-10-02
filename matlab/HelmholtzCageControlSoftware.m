classdef HelmholtzCageControlSoftware_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                       matlab.ui.Figure
        Image2                         matlab.ui.control.Image
        TigerSatsHelmholtzCageSimulatorLabel  matlab.ui.control.Label
        TabGroup                       matlab.ui.container.TabGroup
        OrbitalsimulatorTab            matlab.ui.container.Tab
        GridLayout                     matlab.ui.container.GridLayout
        ComputemagneticfieldButton     matlab.ui.control.Button
        SavemagneticfieldButton        matlab.ui.control.Button
        PlotmagneticfieldButton        matlab.ui.control.Button
        DevdigvijaySingh2023Label      matlab.ui.control.Label
        StatusMessageLabel             matlab.ui.control.Label
        PropagateorbitButton           matlab.ui.control.Button
        StartsecondEditField           matlab.ui.control.NumericEditField
        StartsecondEditFieldLabel      matlab.ui.control.Label
        StartminuteEditField           matlab.ui.control.NumericEditField
        StartminuteEditFieldLabel      matlab.ui.control.Label
        StarthourEditField             matlab.ui.control.NumericEditField
        StarthourEditFieldLabel        matlab.ui.control.Label
        SetISSorbitButton              matlab.ui.control.Button
        VisualizeorbitButton           matlab.ui.control.Button
        SetpolarorbitButton            matlab.ui.control.Button
        MeananomalydegEditField        matlab.ui.control.NumericEditField
        MeananomalydegEditFieldLabel   matlab.ui.control.Label
        OrbitalparametersLabel         matlab.ui.control.Label
        MissionparametersLabel         matlab.ui.control.Label
        PerigeeargdegEditField         matlab.ui.control.NumericEditField
        PerigeeargdegEditFieldLabel    matlab.ui.control.Label
        Eccentricity0circEditField     matlab.ui.control.NumericEditField
        Eccentricity0circEditFieldLabel  matlab.ui.control.Label
        RAANdegEditField               matlab.ui.control.NumericEditField
        RAANdegLabel                   matlab.ui.control.Label
        InclinationidegEditField       matlab.ui.control.NumericEditField
        InclinationidegEditFieldLabel  matlab.ui.control.Label
        ClosestapproachaltkmEditField  matlab.ui.control.NumericEditField
        ClosestapproachaltkmEditFieldLabel  matlab.ui.control.Label
        DurationdaysEditField          matlab.ui.control.NumericEditField
        DurationdaysEditFieldLabel     matlab.ui.control.Label
        MissiontitleEditField          matlab.ui.control.EditField
        MissiontitleEditFieldLabel     matlab.ui.control.Label
        StartdateDatePicker            matlab.ui.control.DatePicker
        StartdateLabel                 matlab.ui.control.Label
        Image                          matlab.ui.control.Image
        LivecoilcontrolsTab            matlab.ui.container.Tab
        GridLayout2                    matlab.ui.container.GridLayout
        GroundtrackPanel               matlab.ui.container.Panel
        CurrentloadedfilereadoutLabel  matlab.ui.control.Label
        CurrentloadedfileLabel         matlab.ui.control.Label
        ResetButton                    matlab.ui.control.Button
        StopButton                     matlab.ui.control.Button
        StartButton                    matlab.ui.control.Button
        LoadtransientfileButton        matlab.ui.control.Button
        LastcommunicationtimereadoutLabel  matlab.ui.control.Label
        LastcommunicationtimeLabel     matlab.ui.control.Label
        LastcommunicationreadoutLabel  matlab.ui.control.Label
        LastcommunicationLabel         matlab.ui.control.Label
        ControlmodeButtonGroup         matlab.ui.container.ButtonGroup
        TransientButton                matlab.ui.control.RadioButton
        StaticButton                   matlab.ui.control.RadioButton
        DisableButton_2                matlab.ui.control.RadioButton
        TransientcontrolLabel          matlab.ui.control.Label
        RestartcontrolsButton          matlab.ui.control.Button
        GridLayout5_2                  matlab.ui.container.GridLayout
        BzEditField_2                  matlab.ui.control.NumericEditField
        BzEditField_2Label             matlab.ui.control.Label
        ByEditField_2                  matlab.ui.control.NumericEditField
        ByEditField_2Label             matlab.ui.control.Label
        BxEditField_2                  matlab.ui.control.NumericEditField
        BxEditField_2Label             matlab.ui.control.Label
        SetstaticButton                matlab.ui.control.Button
        StaticcontrolTLabel            matlab.ui.control.Label
        BasiccommandsLabel             matlab.ui.control.Label
        GridLayout5                    matlab.ui.container.GridLayout
        BzEditField                    matlab.ui.control.NumericEditField
        BzEditFieldLabel               matlab.ui.control.Label
        ByEditField                    matlab.ui.control.NumericEditField
        ByEditFieldLabel               matlab.ui.control.Label
        BxEditField                    matlab.ui.control.NumericEditField
        BxEditFieldLabel               matlab.ui.control.Label
        SetambientButton               matlab.ui.control.Button
        GridLayout4                    matlab.ui.container.GridLayout
        TemperatureReadout             matlab.ui.control.Label
        TemperatureLabel               matlab.ui.control.Label
        CoolingfanLabel                matlab.ui.control.Label
        CoolingfanSwitch               matlab.ui.control.Switch
        AmbientfieldTLabel             matlab.ui.control.Label
        SystemstatusLabel              matlab.ui.control.Label
        SerialcommunicationsLabel      matlab.ui.control.Label
        UpdateCOMButton                matlab.ui.control.Button
        DropDown                       matlab.ui.control.DropDown
        GridLayout3                    matlab.ui.container.GridLayout
        Lamp_2                         matlab.ui.control.Lamp
        ConnectButton                  matlab.ui.control.Button
        UIAxes                         matlab.ui.control.UIAxes
    end


    properties (Access = private)
        % ---------- orbital simulator ---------- %
        satSDP4 % Propagated satellite orbit
        sc      % Satellite scenario
        B_time  % Time series corresponding to magnetic field
        B_field % Magnetic field components in ECEF coordinate system
        B_loc   % Latitude, longitude, and elevation time series

        % ---------- live coils control ---------- %
        BAUD_RATE = 9600;
        CONNECTION_TIMEOUT = 2;
        coils_control % Serial port connection
        valid_file_loaded
        num_fail_cmds
        MAX_FAIL_CMDS = 5;

        % timers
        TEMP_INTERVAL = 1;
        temp_timer
        TRANSIENT_INTERVAL = 0.5;
        transient_timer

        % transient profile
        transient_time
        transient_time_s
        transient_start_time
        transient_Bx
        transient_By
        transient_Bz
        transient_lat
        transient_lon
        transient_h

        % plots
        plt1
        sc1
        yl1
        plt2
        sc2
        current_time
    end

    methods (Access = private)
        function [] = propagate_orbit(app)
            % physical constants
            const_G = 6.67430e-11;
            const_Re = 6.371e6;
            const_Me = 5.97219e24;

            % propagation constants
            const_sample_time = 60;

            % ---- user input parameters ---- %
            param_title = app.MissiontitleEditField.Value;

            % line 1 parameters
            param_start_time = app.StartdateDatePicker.Value;
            param_duration = days(app.DurationdaysEditField.Value);
            % line 2 parameters
            param_closest_alt = app.ClosestapproachaltkmEditField.Value*1e3;  % closest approach altitude, in meters
            param_inclination = app.InclinationidegEditField.Value;     % inclination, in degrees
            param_RAAN = app.RAANdegEditField.Value;             % right ascension of the ascending node, in degrees
            param_ecc = app.Eccentricity0circEditField.Value;              % eccentricity (0: circular)
            param_perigee_arg = app.PerigeeargdegEditField.Value;      % angle from RAAN of closest approach, in degrees
            param_mean_anamoly = app.MeananomalydegEditField.Value;   % mean anomaly - fraction of orbit completed, in degrees

            % ---- depreciated parameters ---- %
            % line 1 constants
            const_catalog_num = -1;     % satellite catalog number
            const_class = 'U';          % launch classification (U - unclassified)
            const_int_des_yr = 20;      % international designation, last two digits of year
            const_int_des_num = -1;     % international designation, launch number of year
            const_int_des_piece = 'A';  % international designation, piece of the launch
            param_mean_motion_dir1 = 0; % 1st derivative of mean motion - ballistic coeff.
            % line 2 constants
            const_epoch_rev_num = 0;    % revolution number at epoch

            % ---- calculated parameters ---- %
            % line 1
            param_stop_date = param_start_time + param_duration;
            % epoch is the start date in this context
            param_epoch_yr = round(mod(year(param_start_time)/100,1)*100);        % last two digits of epoch year (when TLE is true)
            param_epoch_day = day(param_start_time, 'dayofyear') + (hour(param_start_time) + (minute(param_start_time) + second(param_start_time)/60)/60)/24;     % days into year of epoch (when TLE is true)
            % line 2
            param_mean_motion = sqrt(const_G*const_Me/(const_Re + param_closest_alt)^3) * 86400 / (2*pi);      % mean motion, in revolutions per days

            % ---- format TLE ---- %
            tle0 = sprintf('%s', param_title);
            tle1 = sprintf('1 %05d%c %02d%03d%s %02d%012.8f %c.%08d %s %s %s %s', ...
                const_catalog_num, const_class, const_int_des_yr, const_int_des_num, ...
                pad(const_int_des_piece, 3), param_epoch_yr, param_epoch_day, ...
                '+', 0, ...
                '+00000+0', '+00000+0', '0', '00000');
            tle2 = sprintf('2 %05d %08.4f %08.4f %07d %08.4f %08.4f %011.8f %05d', ...
                const_catalog_num, param_inclination, param_RAAN, floor(1e7*param_ecc), ...
                param_perigee_arg, param_mean_anamoly, param_mean_motion, const_epoch_rev_num);
            tle_full = sprintf('%s\n%s\n%s', tle0, tle1, tle2);

            % write to temp file
            tleFile = "temp.tle";
            fileID = fopen(tleFile,'w');
            fprintf(fileID, '%s', tle_full);
            fclose(fileID);

            % ---- simulate orbit ---- %
            app.sc = satelliteScenario(param_start_time,param_stop_date,const_sample_time);

            try
                app.satSDP4 = satellite(app.sc,tleFile, ...
                    "Name",app.MissiontitleEditField.Value, ...
                    "OrbitPropagator","sgp4");
                app.StatusMessageLabel.Text = "Successfully simulated orbit.";

            catch ME
                app.StatusMessageLabel.Text = "Invalid orbit parameters.";
                return;
            end

            enable_compute_magnetic_field(app);
            enable_save_button(app);

        end

        function [] = disable_save_button(app)
            app.VisualizeorbitButton.Enable = false;
            app.ComputemagneticfieldButton.Enable = false;
            app.PlotmagneticfieldButton.Enable = false;
            app.SavemagneticfieldButton.Enable = false;
            app.StatusMessageLabel.Text = "Parameter change detected.";
        end

        function [] = enable_save_button(app)
            app.VisualizeorbitButton.Enable = true;
        end

        function [] = enable_compute_magnetic_field(app)
            app.ComputemagneticfieldButton.Enable = true;
        end

        function [] = enable_save_magnetic_field(app)
            app.PlotmagneticfieldButton.Enable = true;
            app.SavemagneticfieldButton.Enable = true;
        end

        function [] = set_polar_orbit(app)
            app.InclinationidegEditField.Value = 90;

        end

        function [] = set_ISS_orbit(app)
            url = 'https://www.heavens-above.com/orbit.aspx?satid=25544';
            code = webread(url);
            code = code(~isspace(code));

            ecc = str2double(cell2mat(extractBetween(code,'lblEccentricity">','</span>')));
            per_alt = str2double(cell2mat(extractBetween(code,'lblPerigee">','km</span>')));
            RAAN = str2double(cell2mat(extractBetween(code,'lblNode">','</span>')));
            inc = str2double(cell2mat(extractBetween(code,'lblInclination">','</span>')));
            per_arg = str2double(cell2mat(extractBetween(code,'lblArgP">','</span>')));
            mean_anom = str2double(cell2mat(extractBetween(code,'lblMA">','</span>')));

            app.Eccentricity0circEditField.Value = ecc;
            app.ClosestapproachaltkmEditField.Value = per_alt;
            app.RAANdegEditField.Value = RAAN;
            app.InclinationidegEditField.Value = inc;
            app.PerigeeargdegEditField.Value = per_arg;
            app.MeananomalydegEditField.Value = mean_anom;
        end

        function [] = compute_magnetic_field(app)
            try
                [pos,~,time] = states(app.satSDP4);
                pos = pos';
                time = time';

                wgs84 = wgs84Ellipsoid('meter');
                [lat,lon,h] = ecef2geodetic(wgs84,pos(:,1),pos(:,2),pos(:,3));

                B = NaN(length(lat), 3);

                for i = 1:length(lat)
                    [B_temp,~,~,~,~] = wrldmagm(h(i),lat(i),lon(i),decyear(year(time(i)), month(time(i)), day(time(i))));
                    [B(i,1), B(i,2), B(i,3)] = ned2ecefv(B_temp(1),B_temp(2),B_temp(3),lat(i),lon(i));
                end

                app.B_time = time;
                app.B_field = B/1e3;
                app.B_loc = [lat lon h];
            catch ME
                app.StatusMessageLabel.Text = "Magnetic field computation failed.";
                return;
            end

            enable_save_magnetic_field(app);
        end

        function [] = plot_magnetic_field(app)
            figure(1)
            clf;
            plot(app.B_time, app.B_field);
            xlabel('time')
            ylabel('B (\muT)')
            lgd = legend('X', 'Y', 'Z');
            title('magnetic field profile in ECEF coordinate system')
        end

        function [] = save_magnetic_field(app)
            [file,location] = uiputfile(sprintf('%s.csv', app.MissiontitleEditField.Value));
            time = datestr(app.B_time);
            B_x = app.B_field(:,1);
            B_y = app.B_field(:,2);
            B_z = app.B_field(:,3);
            lat = app.B_loc(:,1);
            lon = app.B_loc(:,2);
            h = app.B_loc(:,3);
            data = table(time, B_x, B_y, B_z, lat, lon, h);
            writetable(data, fullfile(location, file));
        end



        function [] = update_COM_ports(app)
            app.ConnectButton.Enable = 'off';
            app.DropDown.Enable = 'off';
            app.DropDown.Items = serialportlist();
            app.ConnectButton.Enable = 'on';
            app.DropDown.Enable = 'on';
        end

        function [] = connect_coils_control(app)
            f = msgbox(sprintf("Attempting connection to %s...", app.DropDown.Value),"Connecting");
            if (isa(app.coils_control,'internal.Serialport'))
                app.coils_control.delete;
            end
            try
                app.DropDown.Value
                app.coils_control = serialport(app.DropDown.Value, app.BAUD_RATE);
                tic;
                t = toc;
                while (t < app.CONNECTION_TIMEOUT && app.coils_control.NumBytesAvailable == 0)
                    pause(0.1)
                    t = toc;
                end
                if (app.coils_control.NumBytesAvailable == 0)
                    errID = 'myComponent:inputError';
                    msgtext = 'Connection error.';
                    ME = MException(errID,msgtext);
                    throw(ME)
                else
                    resp = app.coils_control.readline;
                    if (~strcmp(strtrim(resp), sprintf("<confirmed>")))
                        errID = 'myComponent:inputError';
                        msgtext = 'Connection error.';
                        ME = MException(errID,msgtext);
                        throw(ME)
                    end
                end
            catch ME
                if (exist('f', 'var'))
                    delete(f)
                end
                f = errordlg("Connection failed.");
                end_connection(app);
                return;
            end
            if (exist('f', 'var'))
                delete(f)
            end
            msgbox("Connection successful!","Success");
            start_connection(app);
        end

        function [] = end_connection(app)
            % disable elements

            app.Lamp_2.Enable = 'off';
            app.TemperatureReadout.Enable = 'off';
            app.CoolingfanSwitch.Enable = 'off';
            app.BxEditField.Enable = 'off';
            app.ByEditField.Enable = 'off';
            app.BzEditField.Enable = 'off';
            app.SetambientButton.Enable = 'off';
            app.StartButton.Enable = 'off';
            app.StopButton.Enable = 'off';
            app.ResetButton.Enable = 'off';
            app.ControlmodeButtonGroup.Enable = 'off';
            app.RestartcontrolsButton.Enable = 'off';

            % disable all control modes
            disable_static(app);
            disable_transient(app);

            % stop and delete all timers
            delete(timerfind)
        end

        function [] = start_connection(app)
            % enable all elements

            app.Lamp_2.Enable = 'on';
            app.TemperatureReadout.Enable = 'on';
            app.CoolingfanSwitch.Enable = 'on';
            app.BxEditField.Enable = 'on';
            app.ByEditField.Enable = 'on';
            app.BzEditField.Enable = 'on';
            app.SetambientButton.Enable = 'on';
            app.StartButton.Enable = 'on';
            app.StopButton.Enable = 'on';
            app.ResetButton.Enable = 'on';
            app.ControlmodeButtonGroup.Enable = 'on';
            app.RestartcontrolsButton.Enable = 'on';

            update_control_mode_GUI(app);

            initialize_temperature_timer(app);
        end

        function [] = update_control_mode_GUI(app)
            if (app.ControlmodeButtonGroup.SelectedObject.Text(1) == 'T')
                % transient
                if (app.valid_file_loaded)
                    send_command(app,'<sc1>');
                    disable_static(app);
                    enable_transient(app);
                end
            elseif (app.ControlmodeButtonGroup.SelectedObject.Text(1) == 'S')
                % static
                send_command(app,'<sc1>');
                enable_static(app);
                disable_transient(app);
            else
                % disable
                send_command(app,'<sc0>');
                disable_static(app);
                disable_transient(app);
            end
        end

        function [] = enable_static(app)
            app.BxEditField_2.Enable = 'on';
            app.ByEditField_2.Enable = 'on';
            app.BzEditField_2.Enable = 'on';
            app.SetstaticButton.Enable = 'on';
        end

        function [] = disable_static(app)
            app.BxEditField_2.Enable = 'off';
            app.ByEditField_2.Enable = 'off';
            app.BzEditField_2.Enable = 'off';
            app.SetstaticButton.Enable = 'off';
        end

        function [] = enable_transient(app)
            % check whether a valid profile is loaded prior to calling fcn
            app.StartButton.Enable = 'on';
            app.StopButton.Enable = 'on';
            app.ResetButton.Enable = 'on';
        end

        function [] = disable_transient(app)
            app.StartButton.Enable = 'off';
            app.StopButton.Enable = 'off';
            app.ResetButton.Enable = 'off';
            if (app.valid_file_loaded)
                reset_transient(app);
            end
        end

        function [] = initialize_temperature_timer(app)
            app.temp_timer = timer;
            app.temp_timer.Period = app.TEMP_INTERVAL;
            app.temp_timer.ExecutionMode = 'fixedRate';
            app.temp_timer.TimerFcn = @(~,thisEvent) process_temp(app);
            start(app.temp_timer);
        end

        function [] = process_temp(app)
            resp = send_command(app, '<gt>');
            temp_str = strrep(resp,'<','');
            temp_str = strrep(temp_str, '>', '');
            app.TemperatureReadout.Text = sprintf('%s C', temp_str);
        end

        function [resp] = send_command(app, cmd)
            try
                app.coils_control.flush();
                resp = app.coils_control.writeread(cmd);
                app.num_fail_cmds = 0;
            catch ME
                app.num_fail_cmds = app.num_fail_cmds + 1;
                check_connections(app);
                return;
            end
            app.LastcommunicationreadoutLabel.Text = resp;
            app.LastcommunicationtimereadoutLabel.Text = datestr(datetime('now'));
        end

        function [] = check_connections(app)
            if (app.num_fail_cmds >= app.MAX_FAIL_CMDS)
                end_connection(app);
            end
        end

        function [] = create_plot(app)
            ax = app.UIAxes;
            app.plt1(1) = plot(ax, app.transient_time, app.transient_Bx);
            hold(ax, 'on')
            app.plt1(2) = plot(ax, app.transient_time, app.transient_By);
            app.plt1(3) = plot(ax, app.transient_time, app.transient_Bz);
            app.yl1 = ax.YLim;
            app.sc1 = plot(ax, [0,0], app.yl1, 'k');
            lgd = legend(ax, 'B_x', 'B_y', 'B_z');
            xlabel(ax, 'time')
            ylabel(ax, 'magnetic field (\muT)')
            hold(ax, 'off')

            ax2 = geoaxes(app.GroundtrackPanel);
            app.plt2 = geoplot(ax2, app.transient_lat, app.transient_lon, 'r');
            hold(ax2,'on')
            app.sc2 = geoplot(ax2, app.transient_lat(1), app.transient_lat(1), 'ko');
            hold(ax2,'off')
        end

        function [] = update_plot(app, current_time_s, lat, lon)
            app.sc1.XData = min(app.transient_time) + seconds([current_time_s, current_time_s]);
            app.sc2.XData = lat;
            app.sc2.YData = lon;
        end

        function [] = stop_transient(app)
            if (isa(app.transient_timer, 'timer'))
                delete(app.transient_timer)
            end
        end

        function [] = reset_transient(app)
            if (isa(app.transient_timer, 'timer'))
                delete(app.transient_timer)
            end
            app.sc1.XData = min(app.transient_time) + seconds([0 0]);
            app.sc2.XData = app.transient_lat(1);
            app.sc2.YData = app.transient_lon(1);
            send_command(app, sprintf('<sb%.3f,%.3f,%.3f>', app.transient_Bx(1), app.transient_By(1), app.transient_Bz(1)));
        end

        function [] = update_transient(app)
            ct = seconds(datetime('now') - app.transient_start_time);

            if (ct > max(app.transient_time_s))
                reset_transient(app);
                return;
            end

            Bx_i = interp1(app.transient_time_s, app.transient_Bx, ct);
            By_i = interp1(app.transient_time_s, app.transient_By, ct);
            Bz_i = interp1(app.transient_time_s, app.transient_Bz, ct);

            lat_i = interp1(app.transient_time_s, app.transient_lat, ct);
            lon_i = interp1(app.transient_time_s, app.transient_lon, ct);

            send_command(app, sprintf('<sb%.3f,%.3f,%.3f>', Bx_i, By_i, Bz_i));
            update_plot(app, ct, lat_i, lon_i)
        end

        function [] = start_transient(app)
            if (isa(app.transient_timer, 'timer'))
                delete(app.transient_timer)
            end
            app.transient_timer = timer;
            app.transient_timer.Period = app.TRANSIENT_INTERVAL;
            app.transient_timer.ExecutionMode = 'fixedRate';
            app.transient_timer.TimerFcn = @(~,thisEvent) update_transient(app);
            app.transient_start_time = datetime('now');
            start(app.transient_timer);
            app.StopButton.Enable = 'on';
            app.ResetButton.Enable = 'on';
        end
    end


    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: VisualizeorbitButton
        function VisualizeorbitButtonPushed(app, event)
            v = satelliteScenarioViewer(app.sc);
            play(app.sc)
            camtarget(v,app.satSDP4);
        end

        % Value changed function: StartsecondEditField
        function StartsecondEditFieldValueChanged(app, event)
            disable_save_button(app);
        end

        % Value changed function: StartminuteEditField
        function StartminuteEditFieldValueChanged(app, event)
            disable_save_button(app);
        end

        % Value changed function: StarthourEditField
        function StarthourEditFieldValueChanged(app, event)
            disable_save_button(app);
        end

        % Value changed function: MeananomalydegEditField
        function MeananomalydegEditFieldValueChanged(app, event)
            disable_save_button(app);
        end

        % Value changed function: PerigeeargdegEditField
        function PerigeeargdegEditFieldValueChanged(app, event)
            disable_save_button(app);
        end

        % Value changed function: Eccentricity0circEditField
        function Eccentricity0circEditFieldValueChanged(app, event)
            disable_save_button(app);
        end

        % Value changed function: RAANdegEditField
        function RAANdegEditFieldValueChanged(app, event)
            disable_save_button(app);
        end

        % Value changed function: InclinationidegEditField
        function InclinationidegEditFieldValueChanged(app, event)
            disable_save_button(app);
        end

        % Value changed function: ClosestapproachaltkmEditField
        function ClosestapproachaltkmEditFieldValueChanged(app, event)
            disable_save_button(app);
        end

        % Value changed function: DurationdaysEditField
        function DurationdaysEditFieldValueChanged(app, event)
            disable_save_button(app);
        end

        % Value changed function: MissiontitleEditField
        function MissiontitleEditFieldValueChanged(app, event)
            disable_save_button(app);
        end

        % Value changed function: StartdateDatePicker
        function StartdateDatePickerValueChanged(app, event)
            disable_save_button(app);
        end

        % Button pushed function: PropagateorbitButton
        function PropagateorbitButtonPushed(app, event)
            propagate_orbit(app);
        end

        % Button pushed function: PlotmagneticfieldButton
        function PlotmagneticfieldButtonPushed(app, event)
            plot_magnetic_field(app);
        end

        % Button pushed function: ComputemagneticfieldButton
        function ComputemagneticfieldButtonPushed(app, event)
            compute_magnetic_field(app);
        end

        % Button pushed function: SavemagneticfieldButton
        function SavemagneticfieldButtonPushed(app, event)
            save_magnetic_field(app);
        end

        % Button down function: LivecoilcontrolsTab
        function LivecoilcontrolsTabButtonDown(app, event)
            update_COM_ports(app);
            app.Lamp_2.Enable = 'off';
        end

        % Button pushed function: ConnectButton
        function ConnectButtonPushed2(app, event)
            connect_coils_control(app);
        end

        % Button pushed function: SetISSorbitButton
        function SetISSorbitButtonPushed(app, event)
            set_ISS_orbit(app);
        end

        % Button pushed function: SetpolarorbitButton
        function SetpolarorbitButtonPushed(app, event)
            set_polar_orbit(app);
        end

        % Button pushed function: UpdateCOMButton
        function UpdateCOMButtonPushed2(app, event)
            update_COM_ports(app);
        end

        % Selection changed function: ControlmodeButtonGroup
        function ControlmodeButtonGroupSelectionChanged(app, event)
            selectedButton = app.ControlmodeButtonGroup.SelectedObject;
            update_control_mode_GUI(app);
        end

        % Button pushed function: SetambientButton
        function SetambientButtonPushed(app, event)
            Bx = app.BxEditField.Value;
            By = app.ByEditField.Value;
            Bz = app.BzEditField.Value;

            cmd = sprintf('<sn%.3f,%.3f,%.3f>', Bx,By,Bz);
            send_command(app,cmd);
        end

        % Button pushed function: SetstaticButton
        function SetstaticButtonPushed(app, event)
            Bx = app.BxEditField_2.Value;
            By = app.ByEditField_2.Value;
            Bz = app.BzEditField_2.Value;

            cmd = sprintf('<sb%.3f,%.3f,%.3f>', Bx,By,Bz);
            send_command(app,cmd);
        end

        % Value changed function: CoolingfanSwitch
        function CoolingfanSwitchValueChanged(app, event)
            value = app.CoolingfanSwitch.Value;
            if (value(2) == 'f')
                send_command(app,'<df>')
            else
                send_command(app,'<af>')
            end
        end

        % Button pushed function: RestartcontrolsButton
        function RestartcontrolsButtonPushed(app, event)
            send_command(app,'<rs>')
        end

        % Button pushed function: LoadtransientfileButton
        function LoadtransientfileButtonPushed(app, event)
            try
                [file,location] = uigetfile('data.csv');
                data = readtable(fullfile(location, file));
                app.transient_time = datetime(data.time);
                app.transient_time_s = seconds(app.transient_time - min(app.transient_time));
                app.transient_Bx = data.B_x;
                app.transient_By = data.B_y;
                app.transient_Bz = data.B_z;
                app.transient_lat = data.lat;
                app.transient_lon = data.lon;
                app.transient_h = data.h;

                create_plot(app);
            catch ME
                errordlg("Invalid file input.");
                return;
            end
            app.valid_file_loaded = true;
            app.CurrentloadedfilereadoutLabel.Text = file;
            update_control_mode_GUI(app);
        end

        % Close request function: UIFigure
        function UIFigureCloseRequest(app, event)
            delete(timerfind)
            delete(app)
        end

        % Button pushed function: StartButton
        function StartButtonPushed(app, event)
            start_transient(app);
        end

        % Button pushed function: StopButton
        function StopButtonPushed(app, event)
            stop_transient(app);
        end

        % Button pushed function: ResetButton
        function ResetButtonPushed(app, event)
            reset_transient(app);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Color = [0.94 0.94 0.94];
            app.UIFigure.Position = [100 100 800 640];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.CloseRequestFcn = createCallbackFcn(app, @UIFigureCloseRequest, true);

            % Create TabGroup
            app.TabGroup = uitabgroup(app.UIFigure);
            app.TabGroup.Position = [1 1 800 601];

            % Create OrbitalsimulatorTab
            app.OrbitalsimulatorTab = uitab(app.TabGroup);
            app.OrbitalsimulatorTab.Title = 'Orbital simulator';

            % Create GridLayout
            app.GridLayout = uigridlayout(app.OrbitalsimulatorTab);
            app.GridLayout.ColumnWidth = {'1x', '1x', '1x', '1x'};
            app.GridLayout.RowHeight = {'1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x'};
            app.GridLayout.BackgroundColor = [1 1 1];

            % Create Image
            app.Image = uiimage(app.GridLayout);
            app.Image.Layout.Row = [1 10];
            app.Image.Layout.Column = [3 4];
            app.Image.ImageSource = 'orbital_params.png';

            % Create StartdateLabel
            app.StartdateLabel = uilabel(app.GridLayout);
            app.StartdateLabel.HorizontalAlignment = 'right';
            app.StartdateLabel.Layout.Row = 3;
            app.StartdateLabel.Layout.Column = 1;
            app.StartdateLabel.Text = 'Start date';

            % Create StartdateDatePicker
            app.StartdateDatePicker = uidatepicker(app.GridLayout);
            app.StartdateDatePicker.ValueChangedFcn = createCallbackFcn(app, @StartdateDatePickerValueChanged, true);
            app.StartdateDatePicker.Layout.Row = 3;
            app.StartdateDatePicker.Layout.Column = 2;

            % Create MissiontitleEditFieldLabel
            app.MissiontitleEditFieldLabel = uilabel(app.GridLayout);
            app.MissiontitleEditFieldLabel.HorizontalAlignment = 'right';
            app.MissiontitleEditFieldLabel.Layout.Row = 2;
            app.MissiontitleEditFieldLabel.Layout.Column = 1;
            app.MissiontitleEditFieldLabel.Text = 'Mission title';

            % Create MissiontitleEditField
            app.MissiontitleEditField = uieditfield(app.GridLayout, 'text');
            app.MissiontitleEditField.ValueChangedFcn = createCallbackFcn(app, @MissiontitleEditFieldValueChanged, true);
            app.MissiontitleEditField.Layout.Row = 2;
            app.MissiontitleEditField.Layout.Column = 2;

            % Create DurationdaysEditFieldLabel
            app.DurationdaysEditFieldLabel = uilabel(app.GridLayout);
            app.DurationdaysEditFieldLabel.HorizontalAlignment = 'right';
            app.DurationdaysEditFieldLabel.Layout.Row = 7;
            app.DurationdaysEditFieldLabel.Layout.Column = 1;
            app.DurationdaysEditFieldLabel.Text = 'Duration (days)';

            % Create DurationdaysEditField
            app.DurationdaysEditField = uieditfield(app.GridLayout, 'numeric');
            app.DurationdaysEditField.ValueChangedFcn = createCallbackFcn(app, @DurationdaysEditFieldValueChanged, true);
            app.DurationdaysEditField.Layout.Row = 7;
            app.DurationdaysEditField.Layout.Column = 2;

            % Create ClosestapproachaltkmEditFieldLabel
            app.ClosestapproachaltkmEditFieldLabel = uilabel(app.GridLayout);
            app.ClosestapproachaltkmEditFieldLabel.HorizontalAlignment = 'right';
            app.ClosestapproachaltkmEditFieldLabel.Layout.Row = 10;
            app.ClosestapproachaltkmEditFieldLabel.Layout.Column = 1;
            app.ClosestapproachaltkmEditFieldLabel.Text = 'Closest approach alt. (km)';

            % Create ClosestapproachaltkmEditField
            app.ClosestapproachaltkmEditField = uieditfield(app.GridLayout, 'numeric');
            app.ClosestapproachaltkmEditField.Limits = [0 10000];
            app.ClosestapproachaltkmEditField.ValueChangedFcn = createCallbackFcn(app, @ClosestapproachaltkmEditFieldValueChanged, true);
            app.ClosestapproachaltkmEditField.Layout.Row = 10;
            app.ClosestapproachaltkmEditField.Layout.Column = 2;

            % Create InclinationidegEditFieldLabel
            app.InclinationidegEditFieldLabel = uilabel(app.GridLayout);
            app.InclinationidegEditFieldLabel.HorizontalAlignment = 'right';
            app.InclinationidegEditFieldLabel.Layout.Row = 11;
            app.InclinationidegEditFieldLabel.Layout.Column = 1;
            app.InclinationidegEditFieldLabel.Text = 'Inclination (i, deg.)';

            % Create InclinationidegEditField
            app.InclinationidegEditField = uieditfield(app.GridLayout, 'numeric');
            app.InclinationidegEditField.Limits = [0 360];
            app.InclinationidegEditField.ValueChangedFcn = createCallbackFcn(app, @InclinationidegEditFieldValueChanged, true);
            app.InclinationidegEditField.Layout.Row = 11;
            app.InclinationidegEditField.Layout.Column = 2;

            % Create RAANdegLabel
            app.RAANdegLabel = uilabel(app.GridLayout);
            app.RAANdegLabel.HorizontalAlignment = 'right';
            app.RAANdegLabel.Layout.Row = 12;
            app.RAANdegLabel.Layout.Column = 1;
            app.RAANdegLabel.Text = 'RAAN (Ω, deg.)';

            % Create RAANdegEditField
            app.RAANdegEditField = uieditfield(app.GridLayout, 'numeric');
            app.RAANdegEditField.Limits = [0 360];
            app.RAANdegEditField.ValueChangedFcn = createCallbackFcn(app, @RAANdegEditFieldValueChanged, true);
            app.RAANdegEditField.Layout.Row = 12;
            app.RAANdegEditField.Layout.Column = 2;

            % Create Eccentricity0circEditFieldLabel
            app.Eccentricity0circEditFieldLabel = uilabel(app.GridLayout);
            app.Eccentricity0circEditFieldLabel.HorizontalAlignment = 'right';
            app.Eccentricity0circEditFieldLabel.Layout.Row = 13;
            app.Eccentricity0circEditFieldLabel.Layout.Column = 1;
            app.Eccentricity0circEditFieldLabel.Text = 'Eccentricity (0 = circ.)';

            % Create Eccentricity0circEditField
            app.Eccentricity0circEditField = uieditfield(app.GridLayout, 'numeric');
            app.Eccentricity0circEditField.Limits = [0 1];
            app.Eccentricity0circEditField.ValueChangedFcn = createCallbackFcn(app, @Eccentricity0circEditFieldValueChanged, true);
            app.Eccentricity0circEditField.Layout.Row = 13;
            app.Eccentricity0circEditField.Layout.Column = 2;

            % Create PerigeeargdegEditFieldLabel
            app.PerigeeargdegEditFieldLabel = uilabel(app.GridLayout);
            app.PerigeeargdegEditFieldLabel.HorizontalAlignment = 'right';
            app.PerigeeargdegEditFieldLabel.Layout.Row = 14;
            app.PerigeeargdegEditFieldLabel.Layout.Column = 1;
            app.PerigeeargdegEditFieldLabel.Text = 'Perigee arg. (ω, deg.)';

            % Create PerigeeargdegEditField
            app.PerigeeargdegEditField = uieditfield(app.GridLayout, 'numeric');
            app.PerigeeargdegEditField.Limits = [0 360];
            app.PerigeeargdegEditField.ValueChangedFcn = createCallbackFcn(app, @PerigeeargdegEditFieldValueChanged, true);
            app.PerigeeargdegEditField.Layout.Row = 14;
            app.PerigeeargdegEditField.Layout.Column = 2;

            % Create MissionparametersLabel
            app.MissionparametersLabel = uilabel(app.GridLayout);
            app.MissionparametersLabel.HorizontalAlignment = 'center';
            app.MissionparametersLabel.FontName = 'Arial';
            app.MissionparametersLabel.FontSize = 16;
            app.MissionparametersLabel.FontWeight = 'bold';
            app.MissionparametersLabel.Layout.Row = 1;
            app.MissionparametersLabel.Layout.Column = [1 2];
            app.MissionparametersLabel.Text = 'Mission parameters';

            % Create OrbitalparametersLabel
            app.OrbitalparametersLabel = uilabel(app.GridLayout);
            app.OrbitalparametersLabel.HorizontalAlignment = 'center';
            app.OrbitalparametersLabel.FontName = 'Arial';
            app.OrbitalparametersLabel.FontSize = 16;
            app.OrbitalparametersLabel.FontWeight = 'bold';
            app.OrbitalparametersLabel.Layout.Row = 8;
            app.OrbitalparametersLabel.Layout.Column = [1 2];
            app.OrbitalparametersLabel.Text = 'Orbital parameters';

            % Create MeananomalydegEditFieldLabel
            app.MeananomalydegEditFieldLabel = uilabel(app.GridLayout);
            app.MeananomalydegEditFieldLabel.HorizontalAlignment = 'right';
            app.MeananomalydegEditFieldLabel.Layout.Row = 15;
            app.MeananomalydegEditFieldLabel.Layout.Column = 1;
            app.MeananomalydegEditFieldLabel.Text = 'Mean anomaly (deg.)';

            % Create MeananomalydegEditField
            app.MeananomalydegEditField = uieditfield(app.GridLayout, 'numeric');
            app.MeananomalydegEditField.Limits = [0 360];
            app.MeananomalydegEditField.ValueChangedFcn = createCallbackFcn(app, @MeananomalydegEditFieldValueChanged, true);
            app.MeananomalydegEditField.Layout.Row = 15;
            app.MeananomalydegEditField.Layout.Column = 2;

            % Create SetpolarorbitButton
            app.SetpolarorbitButton = uibutton(app.GridLayout, 'push');
            app.SetpolarorbitButton.ButtonPushedFcn = createCallbackFcn(app, @SetpolarorbitButtonPushed, true);
            app.SetpolarorbitButton.Layout.Row = 9;
            app.SetpolarorbitButton.Layout.Column = 1;
            app.SetpolarorbitButton.Text = 'Set polar orbit';

            % Create VisualizeorbitButton
            app.VisualizeorbitButton = uibutton(app.GridLayout, 'push');
            app.VisualizeorbitButton.ButtonPushedFcn = createCallbackFcn(app, @VisualizeorbitButtonPushed, true);
            app.VisualizeorbitButton.Enable = 'off';
            app.VisualizeorbitButton.Layout.Row = 14;
            app.VisualizeorbitButton.Layout.Column = 3;
            app.VisualizeorbitButton.Text = 'Visualize orbit';

            % Create SetISSorbitButton
            app.SetISSorbitButton = uibutton(app.GridLayout, 'push');
            app.SetISSorbitButton.ButtonPushedFcn = createCallbackFcn(app, @SetISSorbitButtonPushed, true);
            app.SetISSorbitButton.Layout.Row = 9;
            app.SetISSorbitButton.Layout.Column = 2;
            app.SetISSorbitButton.Text = 'Set ISS orbit';

            % Create StarthourEditFieldLabel
            app.StarthourEditFieldLabel = uilabel(app.GridLayout);
            app.StarthourEditFieldLabel.HorizontalAlignment = 'right';
            app.StarthourEditFieldLabel.Layout.Row = 4;
            app.StarthourEditFieldLabel.Layout.Column = 1;
            app.StarthourEditFieldLabel.Text = 'Start hour';

            % Create StarthourEditField
            app.StarthourEditField = uieditfield(app.GridLayout, 'numeric');
            app.StarthourEditField.Limits = [0 24];
            app.StarthourEditField.ValueChangedFcn = createCallbackFcn(app, @StarthourEditFieldValueChanged, true);
            app.StarthourEditField.Layout.Row = 4;
            app.StarthourEditField.Layout.Column = 2;

            % Create StartminuteEditFieldLabel
            app.StartminuteEditFieldLabel = uilabel(app.GridLayout);
            app.StartminuteEditFieldLabel.HorizontalAlignment = 'right';
            app.StartminuteEditFieldLabel.Layout.Row = 5;
            app.StartminuteEditFieldLabel.Layout.Column = 1;
            app.StartminuteEditFieldLabel.Text = 'Start minute';

            % Create StartminuteEditField
            app.StartminuteEditField = uieditfield(app.GridLayout, 'numeric');
            app.StartminuteEditField.Limits = [0 60];
            app.StartminuteEditField.ValueChangedFcn = createCallbackFcn(app, @StartminuteEditFieldValueChanged, true);
            app.StartminuteEditField.Layout.Row = 5;
            app.StartminuteEditField.Layout.Column = 2;

            % Create StartsecondEditFieldLabel
            app.StartsecondEditFieldLabel = uilabel(app.GridLayout);
            app.StartsecondEditFieldLabel.HorizontalAlignment = 'right';
            app.StartsecondEditFieldLabel.Layout.Row = 6;
            app.StartsecondEditFieldLabel.Layout.Column = 1;
            app.StartsecondEditFieldLabel.Text = 'Start second';

            % Create StartsecondEditField
            app.StartsecondEditField = uieditfield(app.GridLayout, 'numeric');
            app.StartsecondEditField.Limits = [0 60];
            app.StartsecondEditField.ValueChangedFcn = createCallbackFcn(app, @StartsecondEditFieldValueChanged, true);
            app.StartsecondEditField.Layout.Row = 6;
            app.StartsecondEditField.Layout.Column = 2;

            % Create PropagateorbitButton
            app.PropagateorbitButton = uibutton(app.GridLayout, 'push');
            app.PropagateorbitButton.ButtonPushedFcn = createCallbackFcn(app, @PropagateorbitButtonPushed, true);
            app.PropagateorbitButton.Layout.Row = 13;
            app.PropagateorbitButton.Layout.Column = 3;
            app.PropagateorbitButton.Text = 'Propagate orbit';

            % Create StatusMessageLabel
            app.StatusMessageLabel = uilabel(app.GridLayout);
            app.StatusMessageLabel.HorizontalAlignment = 'center';
            app.StatusMessageLabel.Layout.Row = 12;
            app.StatusMessageLabel.Layout.Column = [3 4];
            app.StatusMessageLabel.Text = '';

            % Create DevdigvijaySingh2023Label
            app.DevdigvijaySingh2023Label = uilabel(app.GridLayout);
            app.DevdigvijaySingh2023Label.HorizontalAlignment = 'center';
            app.DevdigvijaySingh2023Label.FontSize = 8;
            app.DevdigvijaySingh2023Label.Layout.Row = 16;
            app.DevdigvijaySingh2023Label.Layout.Column = [1 4];
            app.DevdigvijaySingh2023Label.Text = '© Devdigvijay Singh 2023';

            % Create PlotmagneticfieldButton
            app.PlotmagneticfieldButton = uibutton(app.GridLayout, 'push');
            app.PlotmagneticfieldButton.ButtonPushedFcn = createCallbackFcn(app, @PlotmagneticfieldButtonPushed, true);
            app.PlotmagneticfieldButton.Enable = 'off';
            app.PlotmagneticfieldButton.Layout.Row = 14;
            app.PlotmagneticfieldButton.Layout.Column = 4;
            app.PlotmagneticfieldButton.Text = 'Plot magnetic field';

            % Create SavemagneticfieldButton
            app.SavemagneticfieldButton = uibutton(app.GridLayout, 'push');
            app.SavemagneticfieldButton.ButtonPushedFcn = createCallbackFcn(app, @SavemagneticfieldButtonPushed, true);
            app.SavemagneticfieldButton.Enable = 'off';
            app.SavemagneticfieldButton.Layout.Row = 15;
            app.SavemagneticfieldButton.Layout.Column = 4;
            app.SavemagneticfieldButton.Text = 'Save magnetic field';

            % Create ComputemagneticfieldButton
            app.ComputemagneticfieldButton = uibutton(app.GridLayout, 'push');
            app.ComputemagneticfieldButton.ButtonPushedFcn = createCallbackFcn(app, @ComputemagneticfieldButtonPushed, true);
            app.ComputemagneticfieldButton.Enable = 'off';
            app.ComputemagneticfieldButton.Layout.Row = 13;
            app.ComputemagneticfieldButton.Layout.Column = 4;
            app.ComputemagneticfieldButton.Text = 'Compute magnetic field';

            % Create LivecoilcontrolsTab
            app.LivecoilcontrolsTab = uitab(app.TabGroup);
            app.LivecoilcontrolsTab.Title = 'Live coil controls';
            app.LivecoilcontrolsTab.ButtonDownFcn = createCallbackFcn(app, @LivecoilcontrolsTabButtonDown, true);

            % Create GridLayout2
            app.GridLayout2 = uigridlayout(app.LivecoilcontrolsTab);
            app.GridLayout2.ColumnWidth = {'1x', '1x', '1x', '1x'};
            app.GridLayout2.RowHeight = {'1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x'};
            app.GridLayout2.BackgroundColor = [1 1 1];

            % Create UIAxes
            app.UIAxes = uiaxes(app.GridLayout2);
            title(app.UIAxes, 'Magnetic field')
            xlabel(app.UIAxes, 'X')
            ylabel(app.UIAxes, 'Y')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.FontName = 'Helvetica';
            app.UIAxes.GridLineStyle = '-';
            app.UIAxes.YGrid = 'off';
            app.UIAxes.ZGrid = 'off';
            app.UIAxes.ColorOrder = [0 0.447 0.741;0.85 0.325 0.098;0.929 0.694 0.125;0.494 0.184 0.556;0.466 0.674 0.188;0.301 0.745 0.933;0.635 0.078 0.184];
            app.UIAxes.TitleFontSizeMultiplier = 1.1;
            app.UIAxes.Layout.Row = [10 16];
            app.UIAxes.Layout.Column = [1 2];
            colormap(app.UIAxes, 'parula')
            app.UIAxes.Visible = 'off';

            % Create GridLayout3
            app.GridLayout3 = uigridlayout(app.GridLayout2);
            app.GridLayout3.RowHeight = {'1x'};
            app.GridLayout3.Padding = [0 0 0 0];
            app.GridLayout3.Layout.Row = 1;
            app.GridLayout3.Layout.Column = 4;

            % Create ConnectButton
            app.ConnectButton = uibutton(app.GridLayout3, 'push');
            app.ConnectButton.ButtonPushedFcn = createCallbackFcn(app, @ConnectButtonPushed2, true);
            app.ConnectButton.Layout.Row = 1;
            app.ConnectButton.Layout.Column = 1;
            app.ConnectButton.Text = 'Connect';

            % Create Lamp_2
            app.Lamp_2 = uilamp(app.GridLayout3);
            app.Lamp_2.Enable = 'off';
            app.Lamp_2.Layout.Row = 1;
            app.Lamp_2.Layout.Column = 2;

            % Create DropDown
            app.DropDown = uidropdown(app.GridLayout2);
            app.DropDown.Layout.Row = 1;
            app.DropDown.Layout.Column = 3;

            % Create UpdateCOMButton
            app.UpdateCOMButton = uibutton(app.GridLayout2, 'push');
            app.UpdateCOMButton.ButtonPushedFcn = createCallbackFcn(app, @UpdateCOMButtonPushed2, true);
            app.UpdateCOMButton.Layout.Row = 1;
            app.UpdateCOMButton.Layout.Column = 2;
            app.UpdateCOMButton.Text = 'Update COM';

            % Create SerialcommunicationsLabel
            app.SerialcommunicationsLabel = uilabel(app.GridLayout2);
            app.SerialcommunicationsLabel.HorizontalAlignment = 'right';
            app.SerialcommunicationsLabel.FontSize = 14;
            app.SerialcommunicationsLabel.FontWeight = 'bold';
            app.SerialcommunicationsLabel.Layout.Row = 1;
            app.SerialcommunicationsLabel.Layout.Column = 1;
            app.SerialcommunicationsLabel.Text = 'Serial communications';

            % Create SystemstatusLabel
            app.SystemstatusLabel = uilabel(app.GridLayout2);
            app.SystemstatusLabel.HorizontalAlignment = 'right';
            app.SystemstatusLabel.FontSize = 14;
            app.SystemstatusLabel.FontWeight = 'bold';
            app.SystemstatusLabel.Layout.Row = 2;
            app.SystemstatusLabel.Layout.Column = 1;
            app.SystemstatusLabel.Text = 'System status';

            % Create AmbientfieldTLabel
            app.AmbientfieldTLabel = uilabel(app.GridLayout2);
            app.AmbientfieldTLabel.HorizontalAlignment = 'right';
            app.AmbientfieldTLabel.FontSize = 14;
            app.AmbientfieldTLabel.FontWeight = 'bold';
            app.AmbientfieldTLabel.Layout.Row = 3;
            app.AmbientfieldTLabel.Layout.Column = 1;
            app.AmbientfieldTLabel.Text = 'Ambient field (μT)';

            % Create GridLayout4
            app.GridLayout4 = uigridlayout(app.GridLayout2);
            app.GridLayout4.ColumnWidth = {'1x', '1x', '1x', '1x', '1x'};
            app.GridLayout4.RowHeight = {'1x'};
            app.GridLayout4.Padding = [0 0 0 0];
            app.GridLayout4.Layout.Row = 2;
            app.GridLayout4.Layout.Column = [2 4];

            % Create CoolingfanSwitch
            app.CoolingfanSwitch = uiswitch(app.GridLayout4, 'slider');
            app.CoolingfanSwitch.ValueChangedFcn = createCallbackFcn(app, @CoolingfanSwitchValueChanged, true);
            app.CoolingfanSwitch.Enable = 'off';
            app.CoolingfanSwitch.Layout.Row = 1;
            app.CoolingfanSwitch.Layout.Column = 5;
            app.CoolingfanSwitch.Value = 'On';

            % Create CoolingfanLabel
            app.CoolingfanLabel = uilabel(app.GridLayout4);
            app.CoolingfanLabel.HorizontalAlignment = 'right';
            app.CoolingfanLabel.FontWeight = 'bold';
            app.CoolingfanLabel.Layout.Row = 1;
            app.CoolingfanLabel.Layout.Column = 4;
            app.CoolingfanLabel.Text = 'Cooling fan:';

            % Create TemperatureLabel
            app.TemperatureLabel = uilabel(app.GridLayout4);
            app.TemperatureLabel.HorizontalAlignment = 'right';
            app.TemperatureLabel.FontWeight = 'bold';
            app.TemperatureLabel.Layout.Row = 1;
            app.TemperatureLabel.Layout.Column = 2;
            app.TemperatureLabel.Text = 'Temperature:';

            % Create TemperatureReadout
            app.TemperatureReadout = uilabel(app.GridLayout4);
            app.TemperatureReadout.Enable = 'off';
            app.TemperatureReadout.Layout.Row = 1;
            app.TemperatureReadout.Layout.Column = 3;
            app.TemperatureReadout.Text = '0';

            % Create GridLayout5
            app.GridLayout5 = uigridlayout(app.GridLayout2);
            app.GridLayout5.ColumnWidth = {'1x', '1x', '1x', '1x', '1x', '1x', '1x'};
            app.GridLayout5.RowHeight = {'1x'};
            app.GridLayout5.Padding = [0 0 0 0];
            app.GridLayout5.Layout.Row = 3;
            app.GridLayout5.Layout.Column = [2 4];

            % Create SetambientButton
            app.SetambientButton = uibutton(app.GridLayout5, 'push');
            app.SetambientButton.ButtonPushedFcn = createCallbackFcn(app, @SetambientButtonPushed, true);
            app.SetambientButton.Enable = 'off';
            app.SetambientButton.Layout.Row = 1;
            app.SetambientButton.Layout.Column = 7;
            app.SetambientButton.Text = 'Set';

            % Create BxEditFieldLabel
            app.BxEditFieldLabel = uilabel(app.GridLayout5);
            app.BxEditFieldLabel.HorizontalAlignment = 'right';
            app.BxEditFieldLabel.Layout.Row = 1;
            app.BxEditFieldLabel.Layout.Column = 1;
            app.BxEditFieldLabel.Text = 'Bx:';

            % Create BxEditField
            app.BxEditField = uieditfield(app.GridLayout5, 'numeric');
            app.BxEditField.Enable = 'off';
            app.BxEditField.Layout.Row = 1;
            app.BxEditField.Layout.Column = 2;

            % Create ByEditFieldLabel
            app.ByEditFieldLabel = uilabel(app.GridLayout5);
            app.ByEditFieldLabel.HorizontalAlignment = 'right';
            app.ByEditFieldLabel.Layout.Row = 1;
            app.ByEditFieldLabel.Layout.Column = 3;
            app.ByEditFieldLabel.Text = 'By:';

            % Create ByEditField
            app.ByEditField = uieditfield(app.GridLayout5, 'numeric');
            app.ByEditField.Enable = 'off';
            app.ByEditField.Layout.Row = 1;
            app.ByEditField.Layout.Column = 4;

            % Create BzEditFieldLabel
            app.BzEditFieldLabel = uilabel(app.GridLayout5);
            app.BzEditFieldLabel.HorizontalAlignment = 'right';
            app.BzEditFieldLabel.Layout.Row = 1;
            app.BzEditFieldLabel.Layout.Column = 5;
            app.BzEditFieldLabel.Text = 'Bz:';

            % Create BzEditField
            app.BzEditField = uieditfield(app.GridLayout5, 'numeric');
            app.BzEditField.Enable = 'off';
            app.BzEditField.Layout.Row = 1;
            app.BzEditField.Layout.Column = 6;

            % Create BasiccommandsLabel
            app.BasiccommandsLabel = uilabel(app.GridLayout2);
            app.BasiccommandsLabel.HorizontalAlignment = 'right';
            app.BasiccommandsLabel.FontSize = 14;
            app.BasiccommandsLabel.FontWeight = 'bold';
            app.BasiccommandsLabel.Layout.Row = 4;
            app.BasiccommandsLabel.Layout.Column = 1;
            app.BasiccommandsLabel.Text = 'Basic commands';

            % Create StaticcontrolTLabel
            app.StaticcontrolTLabel = uilabel(app.GridLayout2);
            app.StaticcontrolTLabel.HorizontalAlignment = 'right';
            app.StaticcontrolTLabel.FontSize = 14;
            app.StaticcontrolTLabel.FontWeight = 'bold';
            app.StaticcontrolTLabel.Layout.Row = 7;
            app.StaticcontrolTLabel.Layout.Column = 1;
            app.StaticcontrolTLabel.Text = 'Static control  (μT)';

            % Create GridLayout5_2
            app.GridLayout5_2 = uigridlayout(app.GridLayout2);
            app.GridLayout5_2.ColumnWidth = {'1x', '1x', '1x', '1x', '1x', '1x', '1x'};
            app.GridLayout5_2.RowHeight = {'1x'};
            app.GridLayout5_2.Padding = [0 0 0 0];
            app.GridLayout5_2.Layout.Row = 7;
            app.GridLayout5_2.Layout.Column = [2 4];

            % Create SetstaticButton
            app.SetstaticButton = uibutton(app.GridLayout5_2, 'push');
            app.SetstaticButton.ButtonPushedFcn = createCallbackFcn(app, @SetstaticButtonPushed, true);
            app.SetstaticButton.Enable = 'off';
            app.SetstaticButton.Layout.Row = 1;
            app.SetstaticButton.Layout.Column = 7;
            app.SetstaticButton.Text = 'Set';

            % Create BxEditField_2Label
            app.BxEditField_2Label = uilabel(app.GridLayout5_2);
            app.BxEditField_2Label.HorizontalAlignment = 'right';
            app.BxEditField_2Label.Layout.Row = 1;
            app.BxEditField_2Label.Layout.Column = 1;
            app.BxEditField_2Label.Text = 'Bx:';

            % Create BxEditField_2
            app.BxEditField_2 = uieditfield(app.GridLayout5_2, 'numeric');
            app.BxEditField_2.Enable = 'off';
            app.BxEditField_2.Layout.Row = 1;
            app.BxEditField_2.Layout.Column = 2;

            % Create ByEditField_2Label
            app.ByEditField_2Label = uilabel(app.GridLayout5_2);
            app.ByEditField_2Label.HorizontalAlignment = 'right';
            app.ByEditField_2Label.Layout.Row = 1;
            app.ByEditField_2Label.Layout.Column = 3;
            app.ByEditField_2Label.Text = 'By:';

            % Create ByEditField_2
            app.ByEditField_2 = uieditfield(app.GridLayout5_2, 'numeric');
            app.ByEditField_2.Enable = 'off';
            app.ByEditField_2.Layout.Row = 1;
            app.ByEditField_2.Layout.Column = 4;

            % Create BzEditField_2Label
            app.BzEditField_2Label = uilabel(app.GridLayout5_2);
            app.BzEditField_2Label.HorizontalAlignment = 'right';
            app.BzEditField_2Label.Layout.Row = 1;
            app.BzEditField_2Label.Layout.Column = 5;
            app.BzEditField_2Label.Text = 'Bz:';

            % Create BzEditField_2
            app.BzEditField_2 = uieditfield(app.GridLayout5_2, 'numeric');
            app.BzEditField_2.Enable = 'off';
            app.BzEditField_2.Layout.Row = 1;
            app.BzEditField_2.Layout.Column = 6;

            % Create RestartcontrolsButton
            app.RestartcontrolsButton = uibutton(app.GridLayout2, 'push');
            app.RestartcontrolsButton.ButtonPushedFcn = createCallbackFcn(app, @RestartcontrolsButtonPushed, true);
            app.RestartcontrolsButton.Enable = 'off';
            app.RestartcontrolsButton.Layout.Row = 6;
            app.RestartcontrolsButton.Layout.Column = 4;
            app.RestartcontrolsButton.Text = 'Restart controls';

            % Create TransientcontrolLabel
            app.TransientcontrolLabel = uilabel(app.GridLayout2);
            app.TransientcontrolLabel.HorizontalAlignment = 'right';
            app.TransientcontrolLabel.FontSize = 14;
            app.TransientcontrolLabel.FontWeight = 'bold';
            app.TransientcontrolLabel.Layout.Row = 8;
            app.TransientcontrolLabel.Layout.Column = 1;
            app.TransientcontrolLabel.Text = 'Transient control';

            % Create ControlmodeButtonGroup
            app.ControlmodeButtonGroup = uibuttongroup(app.GridLayout2);
            app.ControlmodeButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @ControlmodeButtonGroupSelectionChanged, true);
            app.ControlmodeButtonGroup.Enable = 'off';
            app.ControlmodeButtonGroup.Title = 'Control mode';
            app.ControlmodeButtonGroup.Layout.Row = [4 6];
            app.ControlmodeButtonGroup.Layout.Column = 2;

            % Create DisableButton_2
            app.DisableButton_2 = uiradiobutton(app.ControlmodeButtonGroup);
            app.DisableButton_2.Text = 'Disable';
            app.DisableButton_2.Position = [11 50 62 22];
            app.DisableButton_2.Value = true;

            % Create StaticButton
            app.StaticButton = uiradiobutton(app.ControlmodeButtonGroup);
            app.StaticButton.Text = 'Static';
            app.StaticButton.Position = [11 28 65 22];

            % Create TransientButton
            app.TransientButton = uiradiobutton(app.ControlmodeButtonGroup);
            app.TransientButton.Text = 'Transient';
            app.TransientButton.Position = [11 6 72 22];

            % Create LastcommunicationLabel
            app.LastcommunicationLabel = uilabel(app.GridLayout2);
            app.LastcommunicationLabel.HorizontalAlignment = 'right';
            app.LastcommunicationLabel.FontWeight = 'bold';
            app.LastcommunicationLabel.Layout.Row = 4;
            app.LastcommunicationLabel.Layout.Column = 3;
            app.LastcommunicationLabel.Text = 'Last communication:';

            % Create LastcommunicationreadoutLabel
            app.LastcommunicationreadoutLabel = uilabel(app.GridLayout2);
            app.LastcommunicationreadoutLabel.Layout.Row = 4;
            app.LastcommunicationreadoutLabel.Layout.Column = 4;
            app.LastcommunicationreadoutLabel.Text = 'Disconnected';

            % Create LastcommunicationtimeLabel
            app.LastcommunicationtimeLabel = uilabel(app.GridLayout2);
            app.LastcommunicationtimeLabel.HorizontalAlignment = 'right';
            app.LastcommunicationtimeLabel.FontWeight = 'bold';
            app.LastcommunicationtimeLabel.Layout.Row = 5;
            app.LastcommunicationtimeLabel.Layout.Column = 3;
            app.LastcommunicationtimeLabel.Text = 'Last communication time:';

            % Create LastcommunicationtimereadoutLabel
            app.LastcommunicationtimereadoutLabel = uilabel(app.GridLayout2);
            app.LastcommunicationtimereadoutLabel.Layout.Row = 5;
            app.LastcommunicationtimereadoutLabel.Layout.Column = 4;
            app.LastcommunicationtimereadoutLabel.Text = 'Disconnected';

            % Create LoadtransientfileButton
            app.LoadtransientfileButton = uibutton(app.GridLayout2, 'push');
            app.LoadtransientfileButton.ButtonPushedFcn = createCallbackFcn(app, @LoadtransientfileButtonPushed, true);
            app.LoadtransientfileButton.Layout.Row = 8;
            app.LoadtransientfileButton.Layout.Column = 2;
            app.LoadtransientfileButton.Text = 'Load transient file';

            % Create StartButton
            app.StartButton = uibutton(app.GridLayout2, 'push');
            app.StartButton.ButtonPushedFcn = createCallbackFcn(app, @StartButtonPushed, true);
            app.StartButton.Enable = 'off';
            app.StartButton.Layout.Row = 9;
            app.StartButton.Layout.Column = 2;
            app.StartButton.Text = 'Start';

            % Create StopButton
            app.StopButton = uibutton(app.GridLayout2, 'push');
            app.StopButton.ButtonPushedFcn = createCallbackFcn(app, @StopButtonPushed, true);
            app.StopButton.Enable = 'off';
            app.StopButton.Layout.Row = 9;
            app.StopButton.Layout.Column = 3;
            app.StopButton.Text = 'Stop';

            % Create ResetButton
            app.ResetButton = uibutton(app.GridLayout2, 'push');
            app.ResetButton.ButtonPushedFcn = createCallbackFcn(app, @ResetButtonPushed, true);
            app.ResetButton.Enable = 'off';
            app.ResetButton.Layout.Row = 9;
            app.ResetButton.Layout.Column = 4;
            app.ResetButton.Text = 'Reset';

            % Create CurrentloadedfileLabel
            app.CurrentloadedfileLabel = uilabel(app.GridLayout2);
            app.CurrentloadedfileLabel.HorizontalAlignment = 'right';
            app.CurrentloadedfileLabel.FontWeight = 'bold';
            app.CurrentloadedfileLabel.Layout.Row = 8;
            app.CurrentloadedfileLabel.Layout.Column = 3;
            app.CurrentloadedfileLabel.Text = 'Current loaded file:';

            % Create CurrentloadedfilereadoutLabel
            app.CurrentloadedfilereadoutLabel = uilabel(app.GridLayout2);
            app.CurrentloadedfilereadoutLabel.Layout.Row = 8;
            app.CurrentloadedfilereadoutLabel.Layout.Column = 4;
            app.CurrentloadedfilereadoutLabel.Text = 'no file loaded';

            % Create GroundtrackPanel
            app.GroundtrackPanel = uipanel(app.GridLayout2);
            app.GroundtrackPanel.TitlePosition = 'centertop';
            app.GroundtrackPanel.Title = 'Ground track';
            app.GroundtrackPanel.Layout.Row = [10 16];
            app.GroundtrackPanel.Layout.Column = [3 4];
            app.GroundtrackPanel.FontWeight = 'bold';

            % Create TigerSatsHelmholtzCageSimulatorLabel
            app.TigerSatsHelmholtzCageSimulatorLabel = uilabel(app.UIFigure);
            app.TigerSatsHelmholtzCageSimulatorLabel.HorizontalAlignment = 'center';
            app.TigerSatsHelmholtzCageSimulatorLabel.FontName = 'Segoe UI';
            app.TigerSatsHelmholtzCageSimulatorLabel.FontSize = 20;
            app.TigerSatsHelmholtzCageSimulatorLabel.Position = [4 606 796 35];
            app.TigerSatsHelmholtzCageSimulatorLabel.Text = 'TigerSats Helmholtz Cage Simulator';

            % Create Image2
            app.Image2 = uiimage(app.UIFigure);
            app.Image2.Position = [762 604 37 37];
            app.Image2.ImageSource = 'logo.png';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = HelmholtzCageControlSoftware_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end