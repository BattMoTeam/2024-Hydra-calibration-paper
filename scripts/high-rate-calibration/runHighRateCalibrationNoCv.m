function runHighRateCalibrationNoCv()
%RUNHIGHRATECALIBRATIONNOCV Calibrate with cutoff termination and no CV.

    useCVswitch = false; %#ok<NASGU>
    calibrationSuffix = '-no-cv'; %#ok<NASGU>

    run(fullfile(fileparts(mfilename('fullpath')), 'runHighRateCalibration.m'));

end
