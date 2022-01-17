import numpy as np
"This file contains there species equation fit data and equation methods "


# oxygen:
cp_param_oxygen = (-8.66817221e-19, 4.73115002e-15, -1.15709215e-11,
                   1.66145697e-08, -1.50422620e-05, 8.23507193e-03,
                   - 2.12067742,   1.10964386e+03)

viscosity_param_oxygen = (1.18758866e-26,  -6.48183635e-23, 1.55753837e-19,
                          -2.18596122e-16,  2.02541399e-13, -1.38130567e-10,
                          1.02085148e-07, -1.50063345e-06)

lambda_param_oxygen = ((-1.58986421e-22, 8.03802084e-19, -1.67882604e-15,
                        1.84325862e-12, -1.11449899e-09, 3.57769046e-07,
                        2.11463976e-05, 8.31514294e-03),
                       (-1.79821722e-22, 9.07145474e-19, -1.89752639e-15,
                        2.10059733e-12, -1.29817069e-09, 4.38946254e-07,
                        -4.88815163e-07, 1.15869249e-02))

# hydrogen
cp_param_hydrogen = (2.75575856e-17, -1.58350769e-13, 3.93319791e-10,
                     -5.48239691e-07, 4.61978745e-04, -2.32058478e-01,
                     6.35361636e+01, 7.25459677e+03)

viscosity_param_hydrogen = (6.62607149e-27, -3.42854972e-23, 7.69171320e-20,
                            -9.86354577e-17, 8.10498717e-14, -4.74743879e-11,
                            3.50618247e-08, 1.14582528e-06)

lambda_param_hydrogen = ((-1.24645824e-21, 6.61764024e-18, -1.52054216e-14,
                         1.99690675e-11, -1.49509678e-08, 5.68819226e-06,
                         -4.82527146e-04, 7.12531055e-02),
                         (-1.27194170e-21, 6.74517329e-18, -1.54782877e-14,
                          2.02945953e-11, -1.51875168e-08, 5.79536974e-06,
                          -5.12222368e-04, 7.61461664e-02))

# nitrogen:
cp_param_nitrogen = (-5.97839654e-18, 2.58035023e-14, -4.51777701e-11,
                     4.09644416e-08, -2.06285776e-05, 6.06476999e-03,
                     -9.88549011e-01, 1.10768881e+03)

viscosity_param_nitrogen = (1.53556927e-26, -8.08312960e-23, 1.85403378e-19,
                            -2.44766231e-16, 2.08881853e-13, -1.27545734e-10,
                            8.53886431e-08, -3.89241700e-07)

lambda_param_nitrogen = ((-2.70705853e-22, 1.16336874e-18, -1.95546587e-15,
                          1.51768486e-12, -4.03713326e-10, -1.15746366e-07,
                          1.46652557e-04, -3.27592873e-03),
                         (-2.89265541e-22, 1.25583526e-18, -2.15215826e-15,
                          1.75051772e-12, -5.71060853e-10, -4.11730776e-08,
                          1.26579523e-04, -1.99639546e-04))


# water:
cp_param_water = (4.91358141e-18, -1.85174687e-14,  2.53707252e-11,
                  -1.22872163e-08, -4.19918931e-06, 7.45080766e-03,
                  -2.56601734e+00, 2.12709233e+03)

viscosity_param_water = (1.45183739e-25, -7.27081451e-22, 1.55686360e-18,
                         -1.85885301e-15, 1.35158418e-12, -6.12002995e-10,
                         2.08858479e-07, -1.90872183e-05)

lambda_param_water = ((-7.69988150e-22, 3.81045861e-18, -7.88736102e-15,
                       8.77358057e-12, -5.61460795e-09, 2.11880777e-06,
                       -3.22515696e-04, 2.21293426e-02),
                      (-9.94179604e-21, 4.66529326e-17, -9.19773736e-14,
                       9.85634165e-11, -6.18830008e-08, 2.27982875e-05,
                       -4.45126170e-03, 3.69981235e-01))


class Gas:

    def __init__(self, name, cp_param, viscosity_param, lambda_param):
        self.name = name
        self.cp_param = cp_param
        self.viscosity_param = viscosity_param
        self.lambda_param = lambda_param

    def calc_cp(self, temp):
        return np.polyval(self.cp_param, temp)

    def calc_visc(self, temp):
        return np.polyval(self.viscosity_param, temp)

    def calc_lambda(self, temp, p):
        lambda_1_bar = np.polyval(self.lambda_param[0], temp)
        lambda_10_bar = np.polyval(self.lambda_param[1], temp)
        return lambda_1_bar + (p - 1.e5) / 9.e5 * (lambda_10_bar-lambda_1_bar)


hydrogen = Gas('H2', cp_param_hydrogen, viscosity_param_hydrogen,
               lambda_param_hydrogen)
oxygen = Gas('O2', cp_param_oxygen, viscosity_param_oxygen, lambda_param_oxygen)
nitrogen = Gas('N2', cp_param_nitrogen, viscosity_param_nitrogen,
               lambda_param_nitrogen)
water = Gas('H2O', cp_param_water, viscosity_param_water, lambda_param_water)

species = [hydrogen, oxygen, nitrogen, water]
