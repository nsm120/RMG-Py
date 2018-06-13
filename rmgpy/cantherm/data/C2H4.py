cantherm_species_dictionary = {'label': 'C2H4',
                               'date': '2018-06-13 14:25',
                               'timestamp': '1528914320.17',
                               'author': '',
                               'SMILES': 'C=C',
                               'adjacencyList': '1 C u0 p0 c0 {2,D} {3,S} {4,S}\n2 C u0 p0 c0 {1,D} {5,S} {6,S}\n3 H u0 p0 c0 {1,S}\n4 H u0 p0 c0 {1,S}\n5 H u0 p0 c0 {2,S}\n6 H u0 p0 c0 {2,S}\n',
                               'InChI': 'InChI=1S/C2H4/c1-2/h1-2H2',
                               'InChIKey': 'VGGSQFUCUMXWEO-UHFFFAOYSA-N',
                               'levelOfTheory': '',
                               'modelChemistry': 'CBS-QB3',
                               'frequencyScaleFactor': 0.99,
                               'HinderedRotors': True,
                               'AtomCorrections': True,
                               'BondCorrections': False,
}

CanthermSpecies(
    label = 'C2H4',
    SMILES = 'C=C',
    adjacencyList = '1 C u0 p0 c0 {2,D} {3,S} {4,S}\n2 C u0 p0 c0 {1,D} {5,S} {6,S}\n3 H u0 p0 c0 {1,S}\n4 H u0 p0 c0 {1,S}\n5 H u0 p0 c0 {2,S}\n6 H u0 p0 c0 {2,S}\n',
    InChI = 'InChI=1S/C2H4/c1-2/h1-2H2',
    conformer = Conformer(
        E0 = (45.76, 'kJ/mol'),
        modes = [
            IdealGasTranslation(mass=(28.0313, 'amu')),
            NonlinearRotor(
                inertia = ([3.41526, 16.6498, 20.065], 'amu*angstrom^2'),
                symmetry = 4,
            ),
            HarmonicOscillator(
                frequencies = ([828.397, 970.652, 977.223, 1052.93, 1233.55, 1367.56, 1465.09, 1672.25, 3098.46, 3111.7, 3165.79, 3193.54], 'cm^-1'),
            ),
        ],
        spinMultiplicity = 1,
        opticalIsomers = 1,
        mass = ([12, 1.00783, 1.00783, 12, 1.00783, 1.00783], 'amu'),
        number = [6, 1, 1, 6, 1, 1],
        coordinates = (
            [
                [0.001723, 0, 0.00107],
                [-0.000181, 0, 1.08623],
                [0.975039, 0, -0.478762],
                [-1.1246, 0, -0.700778],
                [-1.12099, 0, -1.78578],
                [-2.09702, 0, -0.219491],
            ],
            'angstroms',
        ),
    ),
    thermo = NASA(
        polynomials = [
            NASAPolynomial(
                coeffs = [4.10264, -0.00691216, 5.11903e-05, -6.07544e-08, 2.37717e-11, 5504.22, 3.22278],
                Tmin = (10, 'K'),
                Tmax = (764.984, 'K'),
            ),
            NASAPolynomial(
                coeffs = [0.459284, 0.018135, -9.68094e-06, 2.54047e-09, -2.62098e-13, 5886.18, 18.677],
                Tmin = (764.984, 'K'),
                Tmax = (3000, 'K'),
            ),
        ],
        Tmin = (10, 'K'),
        Tmax = (3000, 'K'),
        E0 = (45.7869, 'kJ/mol'),
        Cp0 = (33.2579, 'J/(mol*K)'),
        CpInf = (133.032, 'J/(mol*K)'),
    ),
    symmetryNumber = 1,
    collisionModel = TransportData(epsilon=(280.800, 'K'), sigma=(3.971, 'angstroms')),
    # molecularWeight = (28.05, 'amu'),  # Not specified to test calculation
    energyTransferModel = SingleExponentialDown(alpha0=(0.5718, 'kcal/mol'), T0=(300, 'K'), n=0.85),
    chemkin_thermo_string = 'C2H4                    H   4C   2          G    10.000  3000.000  764.98      1\n 4.59283911E-01 1.81349959E-02-9.68093615E-06 2.54047040E-09-2.62097901E-13    2\n 5.88617953E+03 1.86770188E+01 4.10264306E+00-6.91215990E-03 5.11902701E-05    3\n-6.07543686E-08 2.37716597E-11 5.50421547E+03 3.22278367E+00                   4\n',
)
