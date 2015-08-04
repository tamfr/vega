# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 22:43:21 2015

@author: Mott
"""

import telnetlib#, re

body="399"
begin="2030-JAN-01 10:12"
end="2030-JAN-04 10:12"

#def ephemeris(body, begin, end):
t = telnetlib.Telnet()
t.open('horizons.jpl.nasa.gov', 6775)

expect = ( ( r'Horizons>', body + '\n' ),
           ( r'Select.*E.phemeris.*:', 'E\n'),
           ( r'Observe.*:', 'e\n' ),
           ( r'Coordinate system center.*:', '10\n' ),
           ( r'Reference plane.*:', 'eclip\n'),
           ( r'Starting *CT.* :', begin + '\n' ),
           ( r'Ending *CT.* :', end + '\n' ),
           ( r'Output interval.*:', '1d\n' ),
           ( r'Accept default output.*:', 'y\n'),
           ( r'Select\.\.\. .A.gain.* :', 'X\n' )
)
n=0 
with open('results.txt', 'w') as fp:
    while True:
        try:
            answer = t.expect(list(i[0] for i in expect), 10)
            n+=1
            if n==10:            
                ephemeri = answer                
        except EOFError:
            break
        fp.write(answer[2])
        fp.flush()
        t.write(expect[answer[0]][1])

# Loop finds all instances of 'EC=' as a generator and averages the values that follow 'EC='.
#tot = 0
#i = 0
#for m in re.finditer('EC=', ephemeri[2]): 
#    i+=1    
#    tot = (tot + float(ephemeri[2][m.start() + 4 : m.start() + 25]))  
#
#print tot/i
