import unittest
import numpy as np
from composite_plate import Ply, InputError, Laminae, Laminate

class TestCompositePlateClasses(unittest.TestCase):
    """ Defines a series of tests for the CompositePlate Module
    
    """
    def test_ply_E1_value(self):
        """ Checks to see if the input value for E1 is taken by the class
        """
        test_ply = Ply(E1=1.0,E2=1.0,G12=1.0,nu12=1.0,h=1.0);
        self.assertEqual(test_ply.E1, 1.0)

    def test_ply_h_gt_0(self):
        """ Checks to see that an InputError is raised for a h=0 input
        """
        with self.assertRaises(InputError):
            test_ply = Ply(E1=1.0,E2=1.0,G12=1.0,nu12=1.0,h=0.0);
            
    def test_laminae_matrix_orotropic(self):
        """ Test the laminae matrix against a known example
        
            'Fiber-Reinforced Composites' by Mallick (3rd edition), 
                Example 3.6
            
            Test will check to see that the laminae stiffness matrix 
            matches the expected stiffness matrix from Example 3.6 with 
            a maximum error of 1 decimal place (values in book given to 
            1 decimal place).
        """
        
        ply = Ply(E1=133.4, E2=8.78, nu12=0.26, G12=3.254, h=1.0)   # h is not used, just adding 1 as a placeholder
        laminae = Laminae(ply=ply, theta_rad=(45.0*np.pi/180.0))
        Q_bar = laminae.Q_bar
        Q_bar_expected = np.matrix([[40.11, 33.61, 31.3],
            [33.61, 40.11, 31.3],
            [31.3, 31.3, 34.57]])
        Q_max_diff = np.max(np.abs(Q_bar -Q_bar_expected))
        
        self.assertAlmostEqual(Q_max_diff,0,places=1)
        
        
    def test_laminate_matrix_angleply(self):
        """ Test an angle ply laminate matrix against a known example.
        
            'Fiber-Reinforced Composites' by Mallick (3rd edition), 
                Example 3.7a
                
            Test will check that the laminate stiffness matrix matches 
            the expected stiffness matrices from example 3.7a with a 
            maximum normalized error of 3 decimal places (< 0.1% error)
            
            This test will throw a 'RuntimeWarning: invalid value encountered
            in divide' because of the zero elements in the A, B and D 
            matrices, this is okay, ignore the warning.
        """
        
        ply = Ply(E1=133.4e9, E2=8.78e9, nu12=0.26, G12=3.254e9, h=0.006)   # h is in [mm]
        
        laminae_pos45 = Laminae(ply=ply, theta_rad=(45.0*np.pi/180.0))
        laminae_neg45 = Laminae(ply=ply, theta_rad=(-45.0*np.pi/180.0))
        
        laminae_list = [laminae_pos45, laminae_neg45]
        laminate = Laminate(laminae_list)
        
        A_expected = np.power(10.0,6.0)*np.matrix([[481.32, 403.32, 0.0],            
            [403.32, 481.32, 0.0],
            [0.0, 0.0, 414.84]]);
            
        A_diff_norm_max = np.nanmax(np.abs(A_expected -laminate.A)/laminate.A)

        B_expected = np.power(10.0,3.0)*np.matrix([[0.0, 0.0, -1126.8],
            [0.0, 0.0, -1126.8],
            [-1126.8, -1126.8, 0.0]])
            
        B_diff_norm_max = np.nanmax(np.abs(B_expected -laminate.B)/laminate.B)  
        
        D_expected = np.matrix([[5775.84, 4839.84, 0.0],
            [4839.84, 5775.84, 0.0],
            [0.0, 0.0, 4978.08]])
        D_diff_norm_max = np.nanmax(np.abs(D_expected -laminate.D)/laminate.D)
        
        max_norm_diff = np.max([A_diff_norm_max,B_diff_norm_max,D_diff_norm_max])
        
        self.assertAlmostEqual(max_norm_diff,0.0,places=3)  
        
    def test_laminate_matrix_symmetricbalanced(self):
        """ Test a symmetric balanced ply laminate matrix against a 
            known example.
        
            'Fiber-Reinforced Composites' by Mallick (3rd edition), 
                Example 3.7b
                
            Test will check that the laminate stiffness matrix matches 
            the expected stiffness matrices from example 3.7b with a 
            maximum normalized error of 3 decimal places (< 0.1% error)
            
            Symmetric Laminate should have A16, A26 = 0, B = 0
        """
        
        ply = Ply(E1=133.4e9, E2=8.78e9, nu12=0.26, G12=3.254e9, h=0.006)   # h is in [mm]
        
        laminae_pos45 = Laminae(ply=ply, theta_rad=(45.0*np.pi/180.0))
        laminae_neg45 = Laminae(ply=ply, theta_rad=(-45.0*np.pi/180.0))
        
        laminae_list = [laminae_pos45, laminae_neg45, 
            laminae_neg45, laminae_pos45]
            
        laminate = Laminate(laminae_list)
        
        A_expected = np.power(10.0,6.0)*np.matrix([[962.64, 806.64, 0.0],            
            [806.64, 962.64, 0.0],
            [0.0, 0.0, 829.68]]);
            
        A_diff_norm_max = np.nanmax(np.abs(A_expected -laminate.A)/laminate.A)

        B_expected = np.matrix([[0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0]])
            
        B_diff_norm_max = np.max(np.abs(B_expected -laminate.B))  
        
        D_expected = np.power(10.,3.)*np.matrix([[46.21, 38.72, 27.04],
            [38.72, 46.21, 27.04],
            [27.04, 27.04, 39.82]])
        D_diff_norm_max = np.nanmax(np.abs(D_expected -laminate.D)/laminate.D)
        
        max_norm_diff = np.max([A_diff_norm_max,B_diff_norm_max,D_diff_norm_max])
        
        self.assertAlmostEqual(max_norm_diff,0.0,places=3)  
        
        
    def test_laminate_matrix_inverse(self):
        """Checks the inverse (A1,B1,C1,D1) matrices to see they form properly.
        
            'Fiber-Reinforced Composites' by Mallick (3rd edition), 
                Example 3.13
                
            Test will use data from expample 3.13 to check the "inverse"
            relationship stiffness matrices with a maximum normalized 
            error of 2 decimal places (< 1% error).
        """
        
        ply = Ply(E1=133.4e9, E2=8.78e9, nu12=0.26, G12=3.254e9, h=0.006)   # h is in [mm]
        
        laminae_pos45 = Laminae(ply=ply, theta_rad=(45.0*np.pi/180.0))
        laminae_neg45 = Laminae(ply=ply, theta_rad=(-45.0*np.pi/180.0))
        
        laminae_list = [laminae_pos45, laminae_neg45]
        laminate = Laminate(laminae_list)
        
        A1_expected = np.power(10.0,-9.0)*np.matrix([[7.7385, -5.0715, 0.0],            
            [-5.0715, 7.7385, 0.0],
            [0.0, 0.0, 5.683]]);
            
        A1_diff_norm_max = np.nanmax(np.abs(A1_expected -laminate.A_1)/laminate.A_1)

        B1_expected = np.power(10.0,-9.0)*np.matrix([[0.0, 0.0, 603.54],
            [0.0, 0.0, 603.54],
            [602.74, 602.74, 0.0]])
            
        B1_diff_norm_max = np.nanmax(np.abs(B1_expected -laminate.B_1)/laminate.B_1)  
        
        C1_expected = np.power(10.0,-9.0)*np.matrix([[0.0, 0.0, 602.74],
            [0.0, 0.0, 602.74],
            [603.54, 603.54, 0.0]])
        C1_diff_norm_max = np.nanmax(np.abs(C1_expected -laminate.C_1)/laminate.C_1)
        
        D1_expected = np.power(10.,-4.0)*np.matrix([[6.45, -4.23, 0.0],
            [-4.23, 6.45, 0.0],
            [0.0, 0.0, 4.74]])
        D1_diff_norm_max = np.nanmax(np.abs(D1_expected -laminate.D_1)/laminate.D_1)
        
        max_norm_diff = np.nanmax([A1_diff_norm_max,B1_diff_norm_max,C1_diff_norm_max,D1_diff_norm_max])

        self.assertAlmostEqual(max_norm_diff,0.0,places=2)
        
        
    def test_laminate_applied_force_for_strains(self):
        """Apply a force and calculate the resultant laminate midplane strains.
        
            'Fiber-Reinforced Composites' by Mallick (3rd edition), 
                Example 3.13
                
            Test will check that the resultant strains match the expected 
            strains from example 3.13 with a maximum normalized error of
            2 decimal places (< 1% error)
        """
        
        ply = Ply(E1=133.4e9, E2=8.78e9, nu12=0.26, G12=3.254e9, h=0.006)   # h is in [mm]
        
        laminae_pos45 = Laminae(ply=ply, theta_rad=(45.0*np.pi/180.0))
        laminae_neg45 = Laminae(ply=ply, theta_rad=(-45.0*np.pi/180.0))
        
        laminae_list = [laminae_pos45, laminae_neg45]
        laminate = Laminate(laminae_list)
        
        N = np.matrix.transpose(np.matrix([100.0e3, 0.0, 0.0]));    # N[0] is in [N/m]
        M = np.matrix.transpose(np.matrix([0.0, 0.0, 0.0]));
        strain_dictionary = laminate.applied_stress(N,M)
        
        Epsilon = strain_dictionary['Epsilon']
        Kappa = strain_dictionary['Kappa']
        
        Epsilon_expected = np.matrix.transpose(np.matrix([77.385e-5, -50.715e-5, 0.0]))
        Kappa_expected = np.matrix.transpose(np.matrix([0.0, 0.0, 0.060354]))
        
        Epsilon_diff_norm_max = np.nanmax(np.abs(Epsilon_expected -Epsilon)/Epsilon)
        Kappa_diff_norm_max = np.nanmax(np.abs(Kappa_expected -Kappa)/Kappa)
        
        max_norm_diff = np.nanmax([Epsilon_diff_norm_max, Kappa_diff_norm_max])
        
        self.assertAlmostEqual(max_norm_diff,0.0,places=2)
        
        
    def test_laminate_applied_force_for_laminae_midplane_stress(self):
        """Apply a force to a laminate and determine the midplane stresses 
              of the laminae in the laminate.
              
            'Fiber-Reinforced Composites' by Mallick (3rd edition), 
                Example 3.13
                
            Test will look at the midplane stresses of the two laminae in
            the laminate.
        """
        
        ply = Ply(E1=133.4e9, E2=8.78e9, nu12=0.26, G12=3.254e9, h=0.006)   # h is in [mm]
        
        laminae_pos45 = Laminae(ply=ply, theta_rad=(45.0*np.pi/180.0))
        laminae_neg45 = Laminae(ply=ply, theta_rad=(-45.0*np.pi/180.0))
        
        laminae_list = [laminae_pos45, laminae_neg45]
        laminate = Laminate(laminae_list)
        
        N = np.matrix.transpose(np.matrix([100.0e3, 0.0, 0.0]));    # N[0] is in [N/m]
        M = np.matrix.transpose(np.matrix([0.0, 0.0, 0.0]));
        strain_dictionary = laminate.applied_stress(N,M)
        
        Epsilon = strain_dictionary['Epsilon']
        Kappa = strain_dictionary['Kappa']
        
        laminae_midplane_strains = laminate.laminae_midplane_strain(Epsilon, Kappa)  
        
        laminae_midplane_stresses = laminate.laminae_stress(Epsilon, Kappa)
        
        laminae_midplane_stress_expected = np.power(10.0,6.0)*np.matrix.transpose(np.matrix([8.33, 0.0, 2.09]))
        self.assertMatrixAlmostEqualPercent(laminae_midplane_stresses[0],laminae_midplane_stress_expected)
        
        laminae_midplane_stress_expected = np.power(10.0,6.0)*np.matrix.transpose(np.matrix([8.33, 0.0, -2.09]))
        self.assertMatrixAlmostEqualPercent(laminae_midplane_strains[1],laminae_midplane_stress_expected)
        
        
    def assertMatrixAlmostEqualPercent(self, matrix, matrix_expected, msg=None, percent=0.01, places=7):
        """Will compare a numpy matrix to see if all values are equal within 
        the defined percentage (default = 1%)
        """
        self.assertEqual(matrix.size,matrix_expected.size)
        
        if (matrix_expected == 0.0).any():      # Check any expected = 0 values against a fixed number of decimal places, not a percent (percent can be high when expecting = 0)
            matrix_expected_neq0 = []
            matrix_neq0 = []
            matrix_expected_eq0 = []
            matrix_eq0 = []
            
            matrix_expected_neq0 = matrix_expected.where(matrix_expected != 0.0)
            import pdb; pdb.set_trace()
            matrix_neq0 = [matrix[ind] for ind, x in enumerate(matrix_expected) if matrix_expected != 0.0]
            # FIX THIS TEST!!!
            
            if (matrix_expected_neq0):
                self.assertMatrixAlmostEqualPercent(np.matrix(matrix_neq0), np.matrix(matrix_expected_neq0), msg=msg, percent=percent);
            
            if (matrix_expected_eq0):
                assertMatrixAlmostEqual(np.matrix(matrix_eq0), np.matrix(matrix_expected_eq0), msg=msg, places=places);
            
            
        else:
            diff_norm_max = np.nanmax(np.abs(matrix_expected -matrix)/matrix_expected)
            self.assertAlmostEqual(diff_norm_max,0.0, msg=msg, delta=percent)
        
    
    def assertMatrixAlmostEqual(self, matrix, matrix_expected, places=7, msg=None):
        """Will compare a numpy matrix to see if all values are equal within 
        the defined number of decimal places (default = 7)
        """
        self.assertEqual(matrix.size,matrix_expected.size)
        
        diff_norm_max = np.nanmax(np.abs(matrix_expected -matrix))
        import pdb; pdb.set_trace()
        self.assertAlmostEqual(diff_norm_max,0.0,places=places, msg=msg, delta=delta)
        
        
if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestCompositePlateClasses)
    unittest.TextTestRunner(verbosity=2).run(suite)
