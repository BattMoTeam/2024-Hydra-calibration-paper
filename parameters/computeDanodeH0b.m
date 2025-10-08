function D = computeDanodeH0b(soc)
    
    graphite = [0   1.92E-13
                0.1 1.16E-13
                0.2 7.12E-14
                0.3 4.50E-14
                0.4 3.00E-14
                0.5 2.20E-14
                0.6 1.82E-14
                0.7 1.73E-14
                0.8 1.85E-14
                0.9 2.13E-14
                1   2.57E-14
               ];

    D = interpTable(graphite(:,1), graphite(:,2), soc);

end
