clear all
close all

mrstDebug(0)

data = load('evalObjectiveBattmo.mat');

setupNew = data.setupNew;
opt = data.opt;
extra = data.extra;

schedule = setupNew.schedule;
schedule.step.val = schedule.step.val(1:22);
schedule.step.control = schedule.step.control(1:22);


[wellSols, states] = simulateScheduleAD(setupNew.state0, setupNew.model, schedule, ...
                                            'NonLinearSolver', opt.NonlinearSolver, ...
                                            'OutputMinisteps', false              , ...
                                            'Verbose'        , true        , ...
                                            extra{:});
