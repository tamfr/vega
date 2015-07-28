# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 22:43:21 2015

@author: Mott
"""

import telnetlib

def ephemeris(body, begin, end):
    t = telnetlib.Telnet()
    t.open('horizons.jpl.nasa.gov', 6775)
    
    expect = ( ( r'Horizons>', body + '\n' ),
               ( r'Select.*E.phemeris.*:', 'E\n'),
               ( r'Observe.*:', 'e\n' ),
               ( r'Coordinate system center.*:', '0\n' ),
               ( r'Reference plane.*:', 'eclip\n'),
               ( r'Starting *CT.* :', begin + '\n' ),
               ( r'Ending *CT.* :', end + '\n' ),
               ( r'Output interval.*:', '1d\n' ),
               ( r'Accept default output.*:', 'y\n')
    )
     
    with open('results.txt', 'w') as fp:
        while True:
            try:
                answer = t.expect(list(i[0] for i in expect), 10)
                print answer
            except EOFError:
                break
            fp.write(answer[2])
            fp.flush()
            t.write(expect[answer[0]][1])