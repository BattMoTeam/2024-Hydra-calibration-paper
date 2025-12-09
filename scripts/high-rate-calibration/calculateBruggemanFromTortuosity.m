function bman = calculateBruggemanFromTortuosity(model, jsonstructEC, tortuosityRef)

    am    = 'ActiveMaterial';
    itf   = 'Interface';
    pe    = 'PositiveElectrode';
    ne    = 'NegativeElectrode';
    co    = 'Coating';
    sd    = 'SolidDiffusion';
    ctrl  = 'Control';
    elyte = 'Electrolyte';
    sep   = 'Separator';

    bruggeman = @(vf, tau) -log(tau)/log(vf);

    bman = struct(pe, bruggeman(jsonstructEC.(pe).(co).volumeFraction, tortuosityRef.(pe)), ...
                  ne, bruggeman(jsonstructEC.(ne).(co).volumeFraction, tortuosityRef.(ne)), ...
                  sep, bruggeman(1 - model.(sep).porosity, tortuosityRef.(sep)));

end
