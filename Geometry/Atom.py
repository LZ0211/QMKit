import uuid
from .Cartesian import Cartesian


class Atom(Cartesian):
    periodic_table=((1,"H",1.008,((1))),(2,"He",4.0026,((2))),(3,"Li",6.94,((2),(1))),(4,"Be",9.0122,((2),(2))),(5,"B",10.81,((2),(2,1))),(6,"C",12.011,((2),(2,2))),(7,"N",14.007,((2),(2,3))),(8,"O",15.999,((2),(2,4))),(9,"F",18.998,((2),(2,5))),(10,"Ne",20.18,((2),(2,6))),(11,"Na",22.99,((2),(2,6),(1))),(12,"Mg",24.305,((2),(2,6),(2))),(13,"Al",26.982,((2),(2,6),(2,1))),(14,"Si",28.085,((2),(2,6),(2,2))),(15,"P",30.974,((2),(2,6),(2,3))),(16,"S",32.06,((2),(2,6),(2,4))),(17,"Cl",35.45,((2),(2,6),(2,5))),(18,"Ar",39.95,((2),(2,6),(2,6))),(19,"K",39.098,((2),(2,6),(2,6),(1))),(20,"Ca",40.078,((2),(2,6),(2,6),(2))),(21,"Sc",44.956,((2),(2,6),(2,6,1),(2))),(22,"Ti",47.867,((2),(2,6),(2,6,2),(2))),(23,"V",50.942,((2),(2,6),(2,6,3),(2))),(24,"Cr",51.996,((2),(2,6),(2,6,5),(1))),(25,"Mn",54.938,((2),(2,6),(2,6,5),(2))),(26,"Fe",55.845,((2),(2,6),(2,6,6),(2))),(27,"Co",58.933,((2),(2,6),(2,6,7),(2))),(28,"Ni",58.693,((2),(2,6),(2,6,8),(2))),(29,"Cu",63.546,((2),(2,6),(2,6,10),(1))),(30,"Zn",65.38,((2),(2,6),(2,6,10),(2))),(31,"Ga",69.723,((2),(2,6),(2,6,10),(2,1))),(32,"Ge",72.63,((2),(2,6),(2,6,10),(2,2))),(33,"As",74.922,((2),(2,6),(2,6,10),(2,3))),(34,"Se",78.971,((2),(2,6),(2,6,10),(2,4))),(35,"Br",79.901,((2),(2,6),(2,6,10),(2,5))),(36,"Kr",83.798,((2),(2,6),(2,6,10),(2,6))),(37,"Rb",85.468,((2),(2,6),(2,6,10),(2,6),(1))),(38,"Sr",87.62,((2),(2,6),(2,6,10),(2,6),(2))),(39,"Y",88.906,((2),(2,6),(2,6,10),(2,6,1),(2))),(40,"Zr",91.224,((2),(2,6),(2,6,10),(2,6,2),(2))),(41,"Nb",92.906,((2),(2,6),(2,6,10),(2,6,3),(2))),(42,"Mo",95.95,((2),(2,6),(2,6,10),(2,6,5),(1))),(43,"Tc",98,((2),(2,6),(2,6,10),(2,6,5),(2))),(44,"Ru",101.07,((2),(2,6),(2,6,10),(2,6,7),(1))),(45,"Rh",102.91,((2),(2,6),(2,6,10),(2,6,8),(1))),(46,"Pd",106.42,((2),(2,6),(2,6,10),(2,6,10))),(47,"Ag",107.87,((2),(2,6),(2,6,10),(2,6,10),(1))),(48,"Cd",112.41,((2),(2,6),(2,6,10),(2,6,10),(2))),(49,"In",114.82,((2),(2,6),(2,6,10),(2,6,10),(2,1))),(50,"Sn",118.71,((2),(2,6),(2,6,10),(2,6,10),(2,2))),(51,"Sb",121.76,((2),(2,6),(2,6,10),(2,6,10),(2,3))),(52,"Te",127.6,((2),(2,6),(2,6,10),(2,6,10),(2,4))),(53,"I",126.9,((2),(2,6),(2,6,10),(2,6,10),(2,5))),(54,"Xe",131.29,((2),(2,6),(2,6,10),(2,6,10),(2,6))),(55,"Cs",132.91,((2),(2,6),(2,6,10),(2,6,10),(2,6),(1))),(56,"Ba",137.33,((2),(2,6),(2,6,10),(2,6,10),(2,6),(2))),(57,"La",138.91,((2),(2,6),(2,6,10),(2,6,10),(2,6,1),(2))),(58,"Ce",140.12,((2),(2,6),(2,6,10),(2,6,10),(2,6),(1),(1),(2))),(59,"Pr",140.91,((2),(2,6),(2,6,10),(2,6,10),(2,6),(3),(2))),(60,"Nd",144.24,((2),(2,6),(2,6,10),(2,6,10),(2,6),(4),(2))),(61,"Pm",145,((2),(2,6),(2,6,10),(2,6,10),(2,6),(5),(2))),(62,"Sm",150.36,((2),(2,6),(2,6,10),(2,6,10),(2,6),(6),(2))),(63,"Eu",151.96,((2),(2,6),(2,6,10),(2,6,10),(2,6),(7),(2))),(64,"Gd",157.25,((2),(2,6),(2,6,10),(2,6,10),(2,6),(7),(1),(2))),(65,"Tb",158.93,((2),(2,6),(2,6,10),(2,6,10),(2,6),(9),(2))),(66,"Dy",162.5,((2),(2,6),(2,6,10),(2,6,10),(2,6),(10),(2))),(67,"Ho",164.93,((2),(2,6),(2,6,10),(2,6,10),(2,6),(11),(2))),(68,"Er",167.26,((2),(2,6),(2,6,10),(2,6,10),(2,6),(12),(2))),(69,"Tm",168.93,((2),(2,6),(2,6,10),(2,6,10),(2,6),(13),(2))),(70,"Yb",173.05,((2),(2,6),(2,6,10),(2,6,10),(2,6),(14),(2))),(71,"Lu",174.97,((2),(2,6),(2,6,10),(2,6,10),(2,6),(14),(1),(2))),(72,"Hf",178.49,((2),(2,6),(2,6,10),(2,6,10),(2,6),(14),(2),(2))),(73,"Ta",180.95,((2),(2,6),(2,6,10),(2,6,10),(2,6),(14),(3),(2))),(74,"W",183.84,((2),(2,6),(2,6,10),(2,6,10),(2,6),(14),(4),(2))),(75,"Re",186.21,((2),(2,6),(2,6,10),(2,6,10),(2,6),(14),(5),(2))),(76,"Os",190.23,((2),(2,6),(2,6,10),(2,6,10),(2,6),(14),(6),(2))),(77,"Ir",192.22,((2),(2,6),(2,6,10),(2,6,10),(2,6),(14),(7),(2))),(78,"Pt",195.08,((2),(2,6),(2,6,10),(2,6,10),(2,6),(14),(9),(1))),(79,"Au",196.97,((2),(2,6),(2,6,10),(2,6,10),(2,6),(14),(10),(1))),(80,"Hg",200.59,((2),(2,6),(2,6,10),(2,6,10),(2,6),(14),(10),(2))),(81,"Tl",204.38,((2),(2,6),(2,6,10),(2,6,10),(2,6),(14),(10),(2,1))),(82,"Pb",207.2,((2),(2,6),(2,6,10),(2,6,10),(2,6),(14),(10),(2,2))),(83,"Bi",208.98,((2),(2,6),(2,6,10),(2,6,10),(2,6),(14),(10),(2,3))),(84,"Po",209,((2),(2,6),(2,6,10),(2,6,10),(2,6),(14),(10),(2,4))),(85,"At",210,((2),(2,6),(2,6,10),(2,6,10),(2,6),(14),(10),(2,5))),(86,"Rn",222,((2),(2,6),(2,6,10),(2,6,10),(2,6),(14),(10),(2,2))))

    def __init__(self,element,charge=0,multiplicity=1,pos=(0,0,0),label=''):
        super().__init__(*pos)
        ele = Atom.search(element)
        if ele == None:
            return
            #raise Exception('%s is not a chemical element' % element)
        self.order = ele[0]
        self.symbol = ele[1]
        self.weight = ele[2]
        self.orbit = ele[3]
        self.charge = charge
        if label == '':
            self.label = uuid.uuid1()
        else:
            self.label = label
        self.electrons = self.order-charge
        self.multiplicity = 2 if self.electrons % 2 else 1
        self.frozen = False

    def string(self,temp=" {symbol:>2s} {x:>-22.15f} {y:>-22.15f} {z:>-22.15f}"):
        return temp.format(**self.__dict__)

    def clone(self):
        return Atom(self.symbol,self.charge,self.multiplicity,self.pos())

    def freeze(self):
        self.frozen = True
        return self

    def setCharge(self,charge):
        self.charge = int(charge)
        self.electrons = self.order-charge
        self.multiplicity = 2 if self.electrons % 2 else 1
        return self

    def setMultiplicity(self,multiplicity):
        multiplicity = int(multiplicity)
        if (multiplicity - self.multiplicity) % 2 == 0:
            self.multiplicity = multiplicity
        else:
            return
            #raise Exception('multiplicity of %s is impossible!' % multiplicity)

    def setLabel(self,label):
        self.label = label
        return self
    
    def repWith(self,element):
        ele = Atom.search(element)
        self.order = ele[0]
        self.symbol = ele[1]
        self.weight = ele[2]
        self.orbit = ele[3]
        self.electrons = self.order-self.charge
        self.multiplicity = 2 if self.electrons % 2 else 1
    
    @staticmethod
    def search(symbol):
        symbol = symbol.strip().capitalize()
        for ele in Atom.periodic_table:
            if ele[1] == symbol:
                return ele