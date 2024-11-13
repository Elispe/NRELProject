import numpy as np
np.random.seed(8346)
from dataclasses import dataclass
from constant import Kg, Tsm

omega_syn = 2 * np.pi * 60
nb = 39
ng = 10
nb_single_area = 13


@dataclass
class Generator:
    bus: int
    Sgen: float
    H: float
    D: float
    Tsm: float
    Kg: float


@dataclass
class Line:
    bus_from: int
    bus_to: int
    R: float
    X: float
    B: float
    gen_from: int
    gen_to: int


@dataclass
class Transformer:
    bus_from: int
    bus_to: int
    R: float
    X: float
    gen_from: int
    gen_to: int


def create_dic_generator():
    dic_gen = {'Gen1': Generator(bus=30, Sgen=1000, H=4.200, D=1.5, Tsm=Tsm, Kg=Kg),
               'Gen2': Generator(bus=31, Sgen=1000, H=3.03, D=1.5, Tsm=Tsm, Kg=Kg),
               'Gen3': Generator(bus=32, Sgen=1000, H=3.580, D=1.5, Tsm=Tsm, Kg=Kg),
               'Gen4': Generator(bus=33, Sgen=1000, H=2.860, D=1.5, Tsm=Tsm, Kg=Kg),
               'Gen5': Generator(bus=34, Sgen=600, H=4.333, D=1.5, Tsm=Tsm, Kg=Kg),
               'Gen6': Generator(bus=35, Sgen=1000, H=3.480, D=1.5, Tsm=Tsm, Kg=Kg),
               'Gen7': Generator(bus=36, Sgen=1000, H=2.640, D=1.5, Tsm=Tsm, Kg=Kg),
               'Gen8': Generator(bus=37, Sgen=1000, H=2.430, D=1.5, Tsm=Tsm, Kg=Kg),
               'Gen9': Generator(bus=38, Sgen=1000, H=3.450, D=1.5, Tsm=Tsm, Kg=Kg),
               'Gen10': Generator(bus=39, Sgen=10000, H=5, D=1.5, Tsm=Tsm, Kg=Kg)
               }

    return dic_gen


def create_areaA():
    areaA = {'Gen1': Generator(bus=30, Sgen=1000+np.random.choice([100,-100,200,-200]), H=4.200, D=1.5, Tsm=Tsm*np.random.choice([1,0.8,0.7]), Kg=Kg),
             'Gen8': Generator(bus=37, Sgen=1000+np.random.choice([100,-100,200,-200]), H=2.430, D=1.5, Tsm=Tsm*np.random.choice([1,0.8,0.7]), Kg=Kg),
             'Gen9': Generator(bus=38, Sgen=1000+np.random.choice([100,-100,200,-200]), H=3.450, D=1.5, Tsm=Tsm*np.random.choice([1,0.8,0.7]), Kg=Kg),
             }

    return areaA


def create_areaB():
    areaB = {'Gen2': Generator(bus=31, Sgen=1000+np.random.choice([100,-100,200,-200]), H=3.03, D=1.5, Tsm=Tsm*np.random.choice([1,0.8,0.7]), Kg=Kg),
             'Gen3': Generator(bus=32, Sgen=1000+np.random.choice([100,-100,200,-200]), H=3.580, D=1.5, Tsm=Tsm*np.random.choice([1,0.8,0.7]), Kg=Kg),
             'Gen10': Generator(bus=39, Sgen=10000+np.random.choice([100,-100,200,-200]), H=5, D=1.5, Tsm=Tsm*np.random.choice([1,0.8,0.7]), Kg=Kg)
             }

    return areaB


def create_areaC():
    areaC = {'Gen4': Generator(bus=33, Sgen=1000+np.random.choice([100,-100,200,-200]), H=2.860, D=1.5, Tsm=Tsm*np.random.choice([1,0.8,0.7]), Kg=Kg),
             'Gen5': Generator(bus=34, Sgen=600+np.random.choice([100,-100,200,-200]), H=4.333, D=1.5, Tsm=Tsm*np.random.choice([1,0.8,0.7]), Kg=Kg),
             'Gen6': Generator(bus=35, Sgen=1000+np.random.choice([100,-100,200,-200]), H=3.480, D=1.5, Tsm=Tsm*np.random.choice([1,0.8,0.7]), Kg=Kg),
             'Gen7': Generator(bus=36, Sgen=1000+np.random.choice([100,-100,200,-200]), H=2.640, D=1.5, Tsm=Tsm*np.random.choice([1,0.8,0.7]), Kg=Kg),
             }

    return areaC


def create_dic_line():
    dic_line = {'Line1': Line(bus_from=1, bus_to=2, R=0.0035, X=0.0411, B=0.6987, gen_from=11, gen_to=12),
                'Line2': Line(bus_from=1, bus_to=39, R=0.001, X=0.025, B=0.75, gen_from=11, gen_to=1),
                'Line3': Line(bus_from=2, bus_to=3, R=0.0013, X=0.0151, B=0.2572, gen_from=12, gen_to=13),
                'Line4': Line(bus_from=2, bus_to=25, R=0.007, X=0.0086, B=0.146, gen_from=12, gen_to=25),
                'Line5': Line(bus_from=3, bus_to=4, R=0.0013, X=0.0213, B=0.2214, gen_from=13, gen_to=14),
                'Line6': Line(bus_from=3, bus_to=18, R=0.0011, X=0.0133, B=0.2138, gen_from=13, gen_to=28),
                'Line7': Line(bus_from=4, bus_to=5, R=0.0008, X=0.0128, B=0.1342, gen_from=14, gen_to=15),
                'Line8': Line(bus_from=4, bus_to=14, R=0.0008, X=0.0129, B=0.1382, gen_from=14, gen_to=24),
                'Line9': Line(bus_from=5, bus_to=6, R=0.0002, X=0.0026, B=0.0434, gen_from=15, gen_to=16),
                'Line10': Line(bus_from=5, bus_to=8, R=0.0008, X=0.0112, B=0.1476, gen_from=15, gen_to=18),
                'Line11': Line(bus_from=6, bus_to=7, R=0.0006, X=0.0092, B=0.113, gen_from=16, gen_to=17),
                'Line12': Line(bus_from=6, bus_to=11, R=0.0007, X=0.0082, B=0.1389, gen_from=16, gen_to=21),
                'Line13': Line(bus_from=7, bus_to=8, R=0.0004, X=0.0046, B=0.078, gen_from=17, gen_to=18),
                'Line14': Line(bus_from=8, bus_to=9, R=0.0023, X=0.0363, B=0.3804, gen_from=18, gen_to=19),
                'Line15': Line(bus_from=9, bus_to=39, R=0.001, X=0.025, B=1.2, gen_from=19, gen_to=1),
                'Line16': Line(bus_from=10, bus_to=11, R=0.0004, X=0.0043, B=0.0729, gen_from=20, gen_to=21),
                'Line17': Line(bus_from=10, bus_to=13, R=0.0004, X=0.0043, B=0.0729, gen_from=20, gen_to=23),
                'Line18': Line(bus_from=13, bus_to=14, R=0.0009, X=0.0101, B=0.1723, gen_from=23, gen_to=24),
                'Line19': Line(bus_from=14, bus_to=15, R=0.0018, X=0.0217, B=0.366, gen_from=24, gen_to=25),
                'Line20': Line(bus_from=15, bus_to=16, R=0.0009, X=0.0094, B=0.171, gen_from=25, gen_to=26),
                'Line21': Line(bus_from=16, bus_to=17, R=0.0007, X=0.0089, B=0.1342, gen_from=26, gen_to=27),
                'Line22': Line(bus_from=16, bus_to=19, R=0.0016, X=0.0195, B=0.304, gen_from=26, gen_to=29),
                'Line23': Line(bus_from=16, bus_to=21, R=0.0008, X=0.0135, B=0.2548, gen_from=26, gen_to=31),
                'Line24': Line(bus_from=16, bus_to=24, R=0.0003, X=0.0059, B=0.068, gen_from=26, gen_to=34),
                'Line25': Line(bus_from=17, bus_to=18, R=0.0007, X=0.0082, B=0.1319, gen_from=27, gen_to=28),
                'Line26': Line(bus_from=17, bus_to=27, R=0.0013, X=0.0173, B=0.3216, gen_from=27, gen_to=37),
                'Line27': Line(bus_from=21, bus_to=22, R=0.0008, X=0.014, B=0.2526, gen_from=31, gen_to=32),
                'Line28': Line(bus_from=22, bus_to=23, R=0.0006, X=0.0096, B=0.1846, gen_from=32, gen_to=33),
                'Line29': Line(bus_from=23, bus_to=24, R=0.0022, X=0.035, B=0.361, gen_from=33, gen_to=34),
                'Line30': Line(bus_from=25, bus_to=26, R=0.0032, X=0.0323, B=0.513, gen_from=35, gen_to=36),
                'Line31': Line(bus_from=26, bus_to=27, R=0.0014, X=0.0147, B=0.2396, gen_from=36, gen_to=37),
                'Line32': Line(bus_from=26, bus_to=28, R=0.0043, X=0.0474, B=0.7802, gen_from=36, gen_to=38),
                'Line33': Line(bus_from=26, bus_to=29, R=0.0057, X=0.0625, B=1.029, gen_from=36, gen_to=39),
                'Line34': Line(bus_from=28, bus_to=29, R=0.0014, X=0.0151, B=0.249, gen_from=38, gen_to=39)
                }
    return dic_line


def create_dic_transformer():
    dic_transformer = {'Trafo1': Transformer(bus_from=11, bus_to=12, R=0.0016, X=0.0435, gen_from=21, gen_to=22),
                       'Trafo2': Transformer(bus_from=13, bus_to=12, R=0.0016, X=0.0435, gen_from=23, gen_to=22),
                       'Trafo3': Transformer(bus_from=6, bus_to=31, R=0.0, X=0.025, gen_from=16, gen_to=2),
                       'Trafo4': Transformer(bus_from=10, bus_to=32, R=0.0, X=0.02, gen_from=20, gen_to=3),
                       'Trafo5': Transformer(bus_from=19, bus_to=33, R=0.0007, X=0.0142, gen_from=19, gen_to=4),
                       'Trafo6': Transformer(bus_from=20, bus_to=34, R=0.0009, X=0.018, gen_from=30, gen_to=5),
                       'Trafo7': Transformer(bus_from=22, bus_to=35, R=0.0, X=0.0143, gen_from=32, gen_to=6),
                       'Trafo8': Transformer(bus_from=23, bus_to=36, R=0.0005, X=0.0272, gen_from=33, gen_to=7),
                       'Trafo9': Transformer(bus_from=25, bus_to=37, R=0.0006, X=0.0232, gen_from=35, gen_to=8),
                       'Trafo10': Transformer(bus_from=2, bus_to=30, R=0.0, X=0.0181, gen_from=12, gen_to=10),
                       'Trafo11': Transformer(bus_from=29, bus_to=38, R=0.0008, X=0.0156, gen_from=39, gen_to=9),
                       'Trafo12': Transformer(bus_from=19, bus_to=20, R=0.0007, X=0.0138, gen_from=29, gen_to=30)
                       }
    return dic_transformer


def loadIEEE39bus():
    gen = create_dic_generator()
    line = create_dic_line()
    trafos = create_dic_transformer()
    gen_areaA = create_areaA()
    gen_areaB = create_areaB()
    gen_areaC = create_areaC()
    return gen, line, trafos, gen_areaA, gen_areaB, gen_areaC, ng, nb, nb_single_area
