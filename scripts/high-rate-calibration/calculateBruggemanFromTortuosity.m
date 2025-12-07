function bman = calculateBruggemanFromTortuosity(model, jsonstructEC)

    am    = 'ActiveMaterial';
    itf   = 'Interface';
    pe    = 'PositiveElectrode';
    ne    = 'NegativeElectrode';
    co    = 'Coating';
    sd    = 'SolidDiffusion';
    ctrl  = 'Control';
    elyte = 'Electrolyte';
    sep   = 'Separator';

    printer = @(s) disp(jsonencode(s, 'PrettyPrint', true));

    bruggeman = @(vf, tau) -log(tau)/log(vf);

    % Set Bruggeman coeffs from tortuosity
    tortuosityRef = struct(pe, 3.46, ...
                           ne, 3, ...
                           sep, 4.2);
    fprintf('Reference tortuosities:\n');
    printer(tortuosityRef);

    bman = struct(pe, bruggeman(jsonstructEC.(pe).(co).volumeFraction, tortuosityRef.(pe)), ...
                  ne, bruggeman(jsonstructEC.(ne).(co).volumeFraction, tortuosityRef.(ne)), ...
                  sep, bruggeman(1 - model.(sep).porosity, tortuosityRef.(sep)));

    fprintf('Bruggeman coefficients calculated from reference tortuosities:\n');
    printer(bman);

end
