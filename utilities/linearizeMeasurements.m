function [ corrected_measurements ] = linearizeMeasurements( measurementsS, measurementsW, measurementsB, varargin )
%linearizeMeasurement Estimate linearization function and correction factors and apply them to
%obtained measurements
%   Detailed explanation goes here

p=inputParser;

p.addRequired('measurementsS');
p.addRequired('measurementsW');
p.addRequired('measurementsB');
p.addParameter('Plot', 0);

p.parse(measurementsS, measurementsW, measurementsB, varargin{:});

for i=1:size(measurementsS,1)
    for j=1:size(measurementsS,2)
        
        for measurementNo=1:size(measurementsS, 3)
            
%             if(i==1 || j==1)
%                 m=1;
%                 n=1;
%             else
%                 %                 m=i-1;
%                 %                 n=j-1;
%             end
            
            m=1;
            n=1;
            
            mNo=measurementNo;
            
            
            
            if(mod(i,2)==1 & mod(j,2)==1)
                corrected_measurements(i, j, measurementNo) = (measurementsS(i,j,measurementNo)) ...
                    ./(measurementsW(m,n,mNo));
            end
            
            if(mod(i,2)==1 & mod(j,2)==0)
                corrected_measurements(i, j, measurementNo) = (measurementsS(i,j,measurementNo)) ...
                    ./(measurementsW(m,n+1,mNo));
            end
            
            if(mod(i,2)==0 & mod(j,2)==1)
                corrected_measurements(i, j, measurementNo) = (measurementsS(i,j,measurementNo)) ...
                    ./(measurementsW(m+1,n,mNo));
            end
            
            if(mod(i,2)==0 & mod(j,2)==0)
                corrected_measurements(i, j, measurementNo) = (measurementsS(i,j,measurementNo)) ...
                    ./(measurementsW(m+1,n+1,mNo));
            end
            
        end
    end
end

end

