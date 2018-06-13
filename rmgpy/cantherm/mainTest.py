#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2018 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

import unittest
import os
import shutil

from rmgpy.cantherm.statmech import StatMechJob
from rmgpy.cantherm.thermo import ThermoJob
from rmgpy.cantherm.input import load_species_from_database_file

from rmgpy.statmech.conformer import Conformer
from rmgpy.statmech.vibration import HarmonicOscillator
from rmgpy.transport import TransportData
from rmgpy.pdep.collision import SingleExponentialDown
from rmgpy.thermo import ThermoData
from rmgpy.molecule.translator import toInChIKey
from rmgpy.species import Species
from rmgpy.molecule import Molecule

################################################################################


class MainTest(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        """
        A method that is run before all unit tests in this class.
        """
        spc1 = Species()
        spc1.label = 'C2H4'
        self.path1 = os.path.join(os.path.dirname(os.path.abspath(__file__)),'data','C2H4.py')
        self.job1 = StatMechJob(spc1, self.path1)

        spc2 = Species()
        spc2.molecule = [Molecule(SMILES='S')]
        spc2.label = 'H2S'
        spc2.molecularWeight = (34.0809,'amu')
        spc2.transportData = TransportData(shapeIndex=2, epsilon=(1837.5,'J/mol'), sigma=(3.73,'angstroms'),
                                            dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'),
                                            rotrelaxcollnum=0.0, comment="""PrimaryTransportLibrary"""),
        spc2.energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
        spc2.thermo = ThermoData(
                        Tdata = ([300,400,500,600,800,1000,1500],'K'),
                        Cpdata = ([8.18,8.49,8.89,9.31,10.16,10.93,12.3],'cal/(mol*K)','+|-',[1,1,1,1,1,1,1]),
                        H298 = (-4.9,'kcal/mol','+|-',0.5),
                        S298 = (49.18,'cal/(mol*K)','+|-',1),
                        ),
        E0 = (-30.4684, 'kJ/mol')
        modes = [HarmonicOscillator(frequencies=([1215.16,2520.68,2520.68],'cm^-1'))]
        spc2.conformer = Conformer(E0=E0, modes=modes)
        self.job2 = ThermoJob(species=spc2, thermoClass='NASA')
        self.job2.author = 'C.U. There'
        self.job2.levelOfTheory = 'CCSD(T)-F12a/cc-pVTZ-F12//B3LYP/6-311G(2d,d,p), rotors at B3LYP/6-311G(2d,pd)'
        self.job2.modelChemistry = 'CCSD(T)-F12/cc-pVTZ-F12'
        self.job2.frequencyScaleFactor = 0.975
        self.job2.useHinderedRotors = False
        self.job2.useAtomCorrections = False
        self.job2.useBondCorrections = False
        self.job2.atomEnergies = None
        self.job2.path2 = os.path.join(os.path.dirname(os.path.abspath(__file__)),'data','')
        self.job2.path = os.path.join(self.job2.path2,'SpeciesDatabase','H2S.py')
        self.chem_string = """H2S                     H   2S   1          G    10.000  3000.000  764.98      1
        4.59283911E-01 1.81349959E-02-9.68093615E-06 2.54047040E-09-2.62097901E-13    2
        5.88617953E+03 1.86770188E+01 4.10264306E+00-6.91215990E-03 5.11902701E-05    3
        -6.07543686E-08 2.37716597E-11 5.50421547E+03 3.22278367E+00                   4"""

    def test_save_species_to_database_file(self):
        self.job2.save_species_to_database_file(path=self.job2.path2, chem_string=self.chem_string)
        load_species_from_database_file(self.job2)
        self.assertAlmostEquals(self.job2.species.conformer.E0.value_si, -30468.400)
        self.assertAlmostEquals(self.job2.species.molecularWeight.value_si, 0.034080899)
        self.assertEquals(toInChIKey(self.job2.species.molecule[0]), 'RWSOTUBLDIXVET-UHFFFAOYSA-N')
        self.assertEquals(self.job2.species.SMILES, 'S')
        self.assertAlmostEquals(self.job2.species.conformer.modes[0].frequencies.value_si[2], 2520.67999, 4)

    def test_load_species_from_database_file(self):
        load_species_from_database_file(self.job1)
        self.assertEquals(self.job1.species.conformer.E0.value_si, 45760.0)
        self.assertAlmostEquals(self.job1.species.molecularWeight.value_si, 0.0280532)
        self.assertEquals(self.job1.species.InChI, 'InChI=1S/C2H4/c1-2/h1-2H2')
        self.assertEquals(self.job1.species.SMILES, 'C=C')
        self.assertAlmostEquals(self.job1.species.conformer.modes[2].frequencies.value_si[0], 828.397)

    @classmethod
    def tearDownClass(self):
        """
        A method that is run after all unit tests in this class.
        """
        shutil.rmtree(os.path.join(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                'data','SpeciesDatabase','')))

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
