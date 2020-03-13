from Geometry import PBC
import sys,os

file = sys.argv[-1]
dirname = os.path.dirname(file)
filename = os.path.basename(file).split('.')[0]

pbc = PBC.parseGJF(open(file).read())
pbc.name=filename
open(os.path.join(dirname,'POSCAR'),'w+').write(pbc.toVASP())