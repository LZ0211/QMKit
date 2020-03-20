from Geometry import PBC
import sys,os


if __name__ == "__main__":
    file = sys.argv[-1]
    cmd = sys.argv[-2]
    dirname = os.path.dirname(file)
    filename = os.path.basename(file).split('.')[0]
    if cmd == 'gjf2vasp':
        pbc = PBC.parseGJF(open(file).read())
        pbc.name=filename
        open(os.path.join(dirname,'POSCAR'),'w+').write(pbc.toVASP())
    elif cmd == 'vasp2gjf':
        pbc = PBC.parseVASP(open(file).read())
        open(os.path.join(dirname,pbc.name),'w+').write(pbc.toGJF())


