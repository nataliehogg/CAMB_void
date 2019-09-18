import getdist.mcsamples as g
import os
from collections import deque
from datetime import datetime

chain = 'chains/plik18_LSS_fixfid'

file = open('gelman_rubin.txt', 'a+')

g = g.loadMCSamples(chain)
r = g.getConvergeTests(test_confidence=0.95, what=['GelmanRubin'], feedback=True)
file.write('The R-1 of ')
file.write(chain)
file.write(' at ')
file.write(datetime.now().strftime('%H:%M:%S'))
file.write(' on ')
file.write(datetime.now().strftime('%d/%m/%y'))
file.write(' is ')
file.writelines(deque(r, 10))
file.close()
