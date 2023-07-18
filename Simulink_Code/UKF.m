classdef UKF < matlab.System
    properties                  %Sources:    Main Algorithm: https://kodlab.seas.upenn.edu/uploads/Arun/UKFpaper.pdf
                                % assistance understanding points + weights
                                % https://www.cs.unc.edu/~welch/kalman/media/pdf/Julier1997_SPIE_KF.pdf
                                % how sensor data works and process noise
                                % https://ahrs.readthedocs.io/en/latest/filters/aqua.html
        kappa = 0 
        alfa = 0.001
        beta = 6.0
        
        x_0 = 0.001*rand(7,1);
        Q_ = 0.1;
        R_ = 0.1;
    end

    properties (Access = protected)
        lambda_
        gamma
        W
        W0m
        W0c
        matrix_W
        Wi_prime
        x_apriori
        x_aposteriori
        v
        z_apriori
        y
        P_apriori
        P_aposteriori
        sP_aposteriori
        z_sigma
        y_sigma
        x_sigma
        P_xz
        P_vv
        P_zz
        K
        Q
        R
        prev_time
        time_int
        z
        prev_z
    end
    

    methods (Access = protected)
        function setupImpl(obj)
            obj.lambda_ = (6 + obj.kappa) * obj.alfa^2 - 6;
            %obj.gamma = sqrt(6 * obj.lambda_);
            obj.W0m = obj.lambda_ / (6 + obj.lambda_);
            obj.W0c = obj.lambda_ / (6 + obj.lambda_) + (1 - obj.alfa^2 + obj.beta);
            obj.W = 1 / (2 * (6 + obj.lambda_));

            
            obj.x_apriori = [0;0;0;0;0;0];
            obj.x_aposteriori = [0;0;0;0;0;0];

            obj.x_sigma = zeros(7,13);
            obj.y_sigma = zeros(7,13);
            obj.z_sigma = zeros(3,13);
            
            obj.P_apriori = zeros(6);
            obj.P_aposteriori = zeros(6);
            obj.sP_aposteriori = zeros(6);
            obj.matrix_W = zeros(6,12);
            obj.Wi_prime = zeros(6,13);
            obj.P_vv = zeros(3,3);
            obj.K = zeros(6,3);

            obj.prev_time = -0.001;
%             obj.z_apriori = [5;0;0;0;0;0];
%             obj.prev_z = [5;0;0;0;0;0];
%             obj.z = [0;0;0;0;0;0];
%             obj.v = [0;0;0;0;0;0];

            obj.R = 0.003*eye(3); %, zeros(3); zeros(3), 0.00033*eye(3)]; %0.033 degrees/second is the sensor noise from the datasheete
            
        end

        function [out, P] = getOutputSizeImpl(obj)
            % Return size for each output port
            out = [6,1];
            P = [3,1];
        end

        function [out, P] = getOutputDataTypeImpl(obj)
            out = 'double';
            P = 'double';
        end


        function [out, P] = isOutputComplexImpl(obj)
            out = 1;
            P = 1;
        end

        
        function [out, P] = isOutputFixedSizeImpl(obj)
            % Return true for each output port with fixed size
            
            out = 1;
            P = 1;
        end

 

        

%         function resetImpl(obj, Q, R, x_0)
%             % Q - filter process noise covraiance
%             % R - measurement noise covariance,
%             % P - init covariance noise
% 
%             % init of all vectors and matrices where the first dim := n
%             obj.y = zeros(6.00,1);
%             obj.y_P = zeros(6.00,1);
%             obj.P_y = zeros(6.00);
%             obj.P_y_P = zeros(6.00);
%             obj.P_xy = zeros(7.0, 6.00);
%             obj.P_xyP = zeros(7.0, 6.00);
% 
%             obj.K = zeros(7.0, 6.00);
%             obj.K_0 = zeros(7.0, 6.00);
%             obj.K_UKF_T = zeros(6.00, 7.0);
%             obj.y_sigma = zeros(6.00, 2*7.0+1);
%             obj.x_sigma = zeros(7.0, 2*7.0+1);
%             obj.x_sigma_f = zeros(7.0, 2*7.0+1);
% 
%             obj.P_apriori = zeros(7.0);
%             obj.P_aprioriP = zeros(7.0);
%             obj.P_aposteriori = zeros(7.0);
% 
%             obj.x_apriori = x_0(:, 1);
%             obj.x_aposteriori = x_0(:, 1);
%             obj.x_P = zeros(7.0, 1);
% 
%             for i=1:7.0
%                 obj.P_apriori(i,i) = Q;
%                 obj.P_aposteriori(i,i) = Q;
%             end
% 
%             obj.setCovariances(Q, R);
%         end

        function setProcessCovariance(obj,time,vect_X)
%             ARW = 0.003; %need to calculate these in the future
%             BI = 0.0012; 
%             sigma_drift_dot = ((2*pi)/log(2))*((BI^2)/ARW);
%             rate_bias_process_noise = (sigma_drift_dot*(time-prev_time))^2;
%             
            rate_sensornoise = 0.033; %degrees/second equivalent 0.033 degrees/second - root mean squared from datasheet
            attitude_noise = (0.5*time*rate_sensornoise)^2;

            obj.Q = [(vect_X(1:3)*attitude_noise+0.3).*eye(3), zeros(3); zeros(3), 0.003*eye(3)];
        end

%         function measurement_process(obj,z,vect_X,vect_Z,time)
% %             if z(3) >= 0
% %                 quat_accel = [sqrt((z(3)+1)/2) -z(2)/(sqrt(2*(z(3)))) z(1)/(sqrt(2*(z(3)))) 0];
% %             else
% %                 quat_accel = [sqrt((z(3)+1)/2) -z(2)/(sqrt(2*(z(3)))) z(1)/(sqrt(2*(z(3)))) 0];
% %             end
% %             rot_z = rotvecd(quaternion(quat_acc)).';
%            z_grad = compact(quaternion(time*z,'rotvecd'));
%            z_curr = vect_X(1:4,1).';
%            
%            obj.z = [vect_X(1:3)+time*z(1:3); z(1:3)];
%            obj.z = z(1:3);
%         end


        function measurement_model(obj,vect_X)
%              vect_quat = [0, 1, 0, 0]; %vector quaternion for acceleration
%              z_acc = quatmultiply(quatmultiply(vect_X(1:4,:).',vect_quat),quatinv(vect_X(1:4,:).'));
           % obj.z_sigma = [vect_X]
%           z_rot = rotvecd(quaternion(vect_X(1:4,:).')).';
%           obj.z_sigma = [z_rot; vect_X(5:7,:)];
            
            obj.z_sigma = vect_X(5:7,:);
        end

        function meas_meancov(obj,matrix_Z)
            obj.z_apriori = zeros(3,1);
            obj.z_apriori = obj.z_apriori + obj.W0m * matrix_Z(:,1);
            for i = 2:13
                obj.z_apriori = obj.z_apriori + obj.W * matrix_Z(:,i);
            end
            obj.P_zz = zeros(3);
            obj.P_zz = obj.P_zz + obj.W0c * ((matrix_Z(:,1)-obj.z_apriori) .* (matrix_Z(:,1)-obj.z_apriori).');
            for i = 2:13
                obj.P_zz = obj.P_zz + obj.W * ((matrix_Z(:,i)-obj.z_apriori) .* (matrix_Z(:,i)-obj.z_apriori).');
            end
            obj.P_vv = obj.P_zz + obj.R;
        end


        function ap_meancov(obj,vect_Y,vect_guess)
            q_t = compact(quaternion(vect_guess(1:3,:).','rotvecd'));
            q_i = (vect_Y(1:4,:)).';

            ei = quatmultiply(q_i,quatconj(q_t));
            
            rot_ei = rotvecd(quaternion(ei));
               
            rot_e = (1/13) * sum(rot_ei,1);

            e = compact(quaternion(rot_e,'rotvecd'));
            
            q_t = quatmultiply(e,q_t); %iterative way of calculating mean

            while not(all(abs(rot_e) < 0.1))
                ei = quatmultiply(q_i,quatinv(q_t));

                rot_ei = rotvecd(quaternion(ei));
               
                rot_e = (1/13) * sum(rot_ei,1);

                e = compact(quaternion(rot_e,'rotvecd'));
            
                q_t = quatmultiply(e,q_t); %iterative way of calculating mean
            end
            
            rot_mean = rotvecd(quaternion(q_t)).';

            ang_vel_mean = obj.W0m*vect_Y(5:7,1);

            for i = 2:13
                ang_vel_mean = ang_vel_mean + obj.W*vect_Y(5:7,i);
            end

            obj.x_apriori = [rot_mean;ang_vel_mean];

            obj.Wi_prime = [rot_ei';vect_Y(5:7,:)-obj.x_apriori(4:6)];
            
            obj.P_apriori = zeros(6);
            obj.P_apriori = obj.P_apriori + obj.W0c * (obj.Wi_prime(:,1) .* obj.Wi_prime(:,1).');
            
            for i = 2:13
                obj.P_apriori = obj.P_apriori + obj.W * (obj.Wi_prime(:,i) .* obj.Wi_prime(:,i).');
            end
        end

        function cross_correlation(obj)
            obj.P_xz = zeros(6,3);
            obj.P_xz = obj.P_xz + obj.W0c*(obj.Wi_prime(:,1) .* (obj.z_sigma(:,1)-obj.z_apriori).');
            for i = 2:13
                obj.P_xz = obj.P_xz + obj.W * (obj.Wi_prime(:,i) .* (obj.z_sigma(:,i)-obj.z_apriori).');
            end
        end


        function process_model(obj,vect_X,matrix_W)
            
            q_k = compact(quaternion(vect_X(1:3,:).','rotvecd')).';
            
            q_w = compact(quaternion(matrix_W(1:3,:).','rotvecd')).';
            
            w_W = [matrix_W(4,:);matrix_W(5,:);matrix_W(6,:)];

            obj.x_sigma = [[q_k,quatmultiply(q_k',q_w').'];[vect_X(4:6,:),vect_X(4:6,:) + w_W]];

        end

        function time_step(obj,vect_X)

            q_k = vect_X(1:4,:);
            
            q_grad = compact(quaternion(obj.time_int*vect_X(5:7,:).','rotvecd')).';
            
            obj.y_sigma = [quatmultiply(q_k',q_grad')';vect_X(5:7,:)];       
        end

        function [y,P] = stepImpl(obj,time,z)
            
            obj.time_int = time-obj.prev_time;
            
            obj.setProcessCovariance(time,obj.x_aposteriori)

            obj.sP_aposteriori = chol((6 + obj.lambda_)*(obj.P_aposteriori + obj.Q),'upper');

            

            for n = 1:6
                obj.matrix_W(:,n) = obj.sP_aposteriori(:,n);
                obj.matrix_W(:,n+6) = -obj.sP_aposteriori(:,n);
            end
            
            obj.process_model(obj.x_aposteriori,obj.matrix_W); %calculate sigma points

            obj.time_step(obj.x_sigma); %time step forward

            obj.ap_meancov(obj.y_sigma,obj.x_aposteriori); %calculate the mean of y sigma and covariance

            obj.measurement_model(obj.y_sigma); %calculate the predicted measurement based on y_sigma

            obj.meas_meancov(obj.z_sigma); %calculate the mean and covariance of Z_sigma (measurement matrix)
            %also calculates innovation covariance

%             obj.measurement_process(z,obj.x_sigma,obj.prev_z,time-obj.prev_time);
            
            obj.v = z - obj.z_apriori;
%            obj.v = [obj.x_aposteriori(1:3)+(time-obj.prev_time)*0.5*(z+obj.prev_z(4:6)); z] - obj.z_apriori; %calculate innovation
            
            

            obj.cross_correlation(); %calculates cross correlation between y_sigma and z_sigma

            obj.K = obj.P_xz/obj.P_vv; %update the Kalman gains

            obj.x_aposteriori = obj.x_apriori + obj.K*obj.v; %update the state

            

            obj.P_aposteriori = obj.P_apriori - obj.K * obj.P_vv * obj.K.'; %update the covariance state
            
%             obj.prev_z = [obj.x_aposteriori(1:3)+(time-obj.prev_time)*0.5*(z+obj.prev_z(4:6)); z];
            
            obj.prev_time = time;

            y = obj.x_aposteriori;
            P = obj.v;
           
        end
    end
end



       
