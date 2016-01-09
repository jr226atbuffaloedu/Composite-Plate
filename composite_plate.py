########################################################################
#   John Robinson                                                      #
#   2015-10-14                                                         #
#                                                                      #
#   Objects needed for composite plate theory.                         #
########################################################################

########################################################################
#   Module Level Imports
import numpy as np
import string
import unittest

class PlyFromFiber(object):
    """ Estimate ply properties from Fiber and Matrix properties
    
    Will provide estimates of various elastic and thermal ply properties
    based on the Fiber and matrix properties
    """
    
    def __init__(self):
        """ Initialize the PlyFromFiber object, all member variables to 
        NaN
        
        This initializes all possible member variables to NaN, since 
        different estimated properties require different variables to be
        defined.  When using, define the required fiber and matrix values 
        for your property of interest, and then calculate the property.
        
        Member Variables:
            E_f:    Elastic modulus of the fiber
            v_f:    Volume fraction of the ply that is fiber
            nu_f:   Poissons ratio of the fiber
            G_f:    Shear modulus of the fiber
            
            E_m:    Elastic modulus of the matrix
            v_m:    Volume fraction of the ply that is matrix
            nu_m:   Poissons ratio of the matrix
            G_m:    Shear modulus of the matrix
        """
        # Fiber Properties
        self.E_f = float('nan')
        self.v_f = float('nan')
        self.nu_f = float('nan')
        self.G_f = float('nan')
        
        # Matrix Properties
        self.E_m = float('nan')
        self.v_m = float('nan')
        self.nu_m = float('nan')
        self.G_m = float('nan')
        
    def elastic_props(self):
        """Calculates the elastic properties of the ply from the fiber 
            and matrix properties.
            
        """
        if (np.isnan([self.E_f,self.v_f,self.nu_f,self.G_f,
            self.E_m,self.v_m,self.nu_m,self.G_m]).any()):
            raise InputError('__init__()','A matrix of fiber elastic constant was not initialized!')
        else:        
            self.E11 = self.E_f*self.v_f +self.E_m*self.v_m
            self.E22 = E_f*E_m/(E_f*v_m +E_m*v_f)
            self.nu12 = self.nu_f*self.v_f +self.nu_m*self.v_m
            self.nu21 = self.E22/self.E11*self.nu12
            self.G12 = self.G_f*self.G_m/(self.G_f*self.v_m +self.G_m*self.v_f)

class Ply(object):    
    """ Defines a ply for a composite.  
    
    A ply is the material which goes into a composite.  Typically a ply 
    will have a different modulus in the primary (E1) direction and the 
    traverse (E2) direction.  The ply will have a shear modulus (G12) 
    and poissons ratio (nu12), as well as a thickness (h)
    
    Example:
        To instantiate a ply, you must specify the plys parameter values::
            
            >>> E1      = 133.44;   # Modulus in fiber direction        [GPa]
            >>> E2      = 8.78;     # Modulus in transverse direction   [GPa]
            >>> G12     = 3.254;    # Shear Modulus                     [GPa]
            >>> nu12    = 0.26;     # Poissons Ratio (principal)        [   ]
            >>> h       = 6;        # Ply thickness                     [mm]
            >>> ply = Ply(E1=E1,E2=E2,G12=G12,nu12=nu12,h=h)
    """
    
    def __repr__(self):
        return("E1:{0},E2:{1},G12:{2},nu12:{3},h:{4}"
            .format(self.E1,self.E2,self.G12,self.nu12,self.h))

    def __init__(self,E1,E2,G12,nu12,h):
        """Initializes the properties of the anisotropic ply.
        
        Anisotropic ply.  Need to define 4 of [E1,E2,nu12,nu21,G12] (4 
        of 5 independent variables).  
        Dependent constraint: nu21/E2 = nu12/E1 
     
        Units can be any set of self consistent units desired.
        
        Args:
            E1 (float):     Modulus in the primary direction
            E2 (float):     Modulus perpendicular to the primary direction
            G12 (float):    Shear Modulus
            nu12 (float):   Poissons Ratio
            h (float):      Ply thickness
            
        Exceptions:
            InputError:     If any input value is not > 0
        """
        if ((E1>0.0) and (E2>0.0) and (G12>0.0) and (nu12>0.0) and (h>0.0)):  # Check to be sure all inputs > 0
            self.E1     = E1
            self.E2     = E2
            self.G12    = G12
            self.nu12   = nu12
            self.h      = h
        else:
            raise InputError('__init__()','An input value is not greater than 0!')
            
    def degree_of_isotropy(self):
        """Returns the degree of Isotropy of the Ply (E1/E2)
        
        E1/E2 is a metric used to compare how isotropic plys are
        """
        
        return self.E1/self.E2

def rotation_matrix(_theta_rad):
    """Rotation Matrix
    """
    c_theta  = np.cos(_theta_rad);
    s_theta  = np.sin(_theta_rad);
    c2_theta = np.power(c_theta,2);
    s2_theta = np.power(s_theta,2);
    s_c_theta = s_theta*c_theta;
    c2minuss2_theta = c2_theta -s2_theta;

    return np.matrix([[c2_theta, s2_theta, 2*s_c_theta],
        [s2_theta, c2_theta, -2*s_c_theta],
        [-s_c_theta, s_c_theta, c2_theta-s2_theta]]);

class Laminae(object):
    """ Composite Laminae defined for each layer in a composite.
    
    A laminae is an individual layer in a composite laminate.  It 
    consists of a ply rotated some angle in the composite.  The angle is
    defined as the angle between the plys primary (1) axis and the 
    composites x axis.
    
    Example:
        To instantiate a laminae, you must start with a ply, and then use
        the ply instance and a specified angle to get a single laminae.::
            
            >>> E1      = 133.44;   # Modulus in fiber direction        [GPa]
            >>> E2      = 8.78;     # Modulus in transverse direction   [GPa]
            >>> G12     = 3.254;    # Shear Modulus                     [GPa]
            >>> nu12    = 0.26;     # Poissons Ratio (principal)        [   ]
            >>> h       = 6;        # Ply thickness                     [mm]
            >>> ply = Ply(E1=E1,E2=E2,G12=G12,nu12=nu12,h=h)
            >>> theta_deg = 45;     # Orientation (+ from X)            [deg]
            >>> theta_rad = theta_deg*np.pi/180.0;    #                 [rad]
            >>> laminae1 = Laminae(ply, theta_rad);
    """
    
    def __repr__(self):
        return ("Ply: {0}, Angle {1} [deg]".format(self.ply,(self.theta_rad*180.0/np.pi)))
    
    def __init__(self, ply,theta_rad):
        """ Initialize the composite Laminae
        
        Defines an individual laminae based on the ply and orientation 
        of the laminae x-axis to the ply 1 axis (theta_rad)
        
        Args:
            ply (composite_plate.Ply):        A ply object
            theta_rad (float):  The orientation of the plys primary (1) axis to the composites x axis [radians]        
        """
        self.ply = ply;
        self.theta_rad = theta_rad;

        ################################################################
        #   Compliance, Rotation and Reuters Matrices

        #   Compliance Matrix
        S = np.matrix([[1.0/self.ply.E1, -self.ply.nu12/self.ply.E1, 0.0],
            [-self.ply.nu12/self.ply.E1, 1.0/self.ply.E2, 0.0],
            [0.0, 0.0, 1.0/self.ply.G12]]);
            
        #   Rotation Matrix
        A = rotation_matrix(theta_rad);
            
        #   Reuters Matrix (Comes from classical definition of shear 
        #       strain)    
        R = np.matrix([[1, 0, 0],
            [0, 1, 0],
            [0, 0, 2]]);
        
        ################################################################
        #   Source: "Laminated Composite Plates" by David Roylance (MIT)
        #       Feb. 10, 2000    
        
        #   Transformed Compliance Matrix
        self.S_bar = R*np.linalg.inv(A)*np.linalg.inv(R)*S*A;

        self.Q_bar = np.linalg.inv(self.S_bar);

class Laminate(object):
    r""" A class defining a composite Laminate.  Input is as the composites 
    laminae, the laminae thickness and stacking order.
    
    Args:
        laminae_list (list of composite_plate.Laminae): A list of laminae which compose the laminate, in order from bottom to top
    
    Attributes:
        A,B,D (numpy.matrix): Compliance matrices will become class 
            variables after __init__.
            
        laminae_midplane_z (list[float]):   A list of heights from the 
            midplane of the laminate to the midplane of each laminae.
            
        Epsilon (numpy.matrix): If defined contains a vector [3x1] of the 
            midplane laminate strains.
            
        Kappa (numpy.matrix): If defined contains a vector [3x1] of the 
            laminate curvatures.
        
    Example:
        To instantiate a laminate, you must start with one or more laminae,
        and then use the laminae instance(s) in a list to create the 
        laminate.::
            
            >>> E1      = 133.44;   # Modulus in fiber direction        [GPa]
            >>> E2      = 8.78;     # Modulus in transverse direction   [GPa]
            >>> G12     = 3.254;    # Shear Modulus                     [GPa]
            >>> nu12    = 0.26;     # Poissons Ratio (principal)        [   ]
            >>> h       = 6;        # Ply thickness                     [mm]
            >>> ply = Ply(E1=E1,E2=E2,G12=G12,nu12=nu12,h=h)
            >>> theta_deg = 45;     # Orientation (+ from X)            [deg]
            >>> theta_rad = theta_deg*np.pi/180.0;    #                 [rad]
            >>> laminae1 = Laminae(ply, theta_rad);            
            >>> laminae2 = Laminae(ply,-theta_rad);
            >>> laminate = Laminate([laminae1, laminae2])
    
    Stress Strain Relationships
    ---------------------------
    
    .. math::
        [N](3x1) &= [A](3x3)*[\epsilon](3x1) +[B](3x3)*[\kappa](3x1) \\
        [M](3x1) &= [B](3x3)*[\epsilon](3x1) +[D](3x3)*[\kappa](3x1)
    
    Normal Force/ unit lenght [F/l] in x,y and xy directions:
    
    .. math:: 
        [N] &= [N_{xx} N_{yy} N_{xy}]^T \\
        [M] &= [M_{xx} M_{yy} M_{xy}]^T 
        
    -> Moment / unit length [F*l/l]
    
    .. math::
        [\epsilon] &= [\epsilon^0_{xx} \epsilon^0_{yy} \gamma^0_{xy}]^T \\
        [\kappa] &= [\kappa_{xx} \kappa_{yy} \kappa_{xy}]^T
    
    Midplane strains indicated by :math:`\epsilon^0`
    
    .. math::
    
        [A]_{mn} = \sum_{j=1}^N(\bar{Q}_{mn})_j(h_j -h_{j-1})\\
        [B]_{mn} = \frac{1}{2}\sum_{j=1}^N(\bar{Q}_{mn})_j(h_j^2 -h_{j-1}^2)\\
        [D]_{mn} = \frac{1}{3}\sum_{j=1}^N(\bar{Q}_{mn})_j(h_j^3 -h_{j-1}^3)\\
        
    Laminae stiffness matrix :math:`\bar{Q}`
    
    .. math::
        [\epsilon] &= [A_1]*[N] + [B_1]*[M] \\
        [\kappa] &= [C_1]*[N] + [D_1]*[M] \\ 
        \\
        [A_1] &= [A^{-1}]+[A^{-1}][B][D^{*-1}][B][A^{-1}] \\
        [B_1] &= -[A^{-1}][B][D^{*-1}] \\
        [C_1] &= -[D^{*-1}][B][A^{-1}] = [B_1]^T \\
        [D^*] &= [D] -[B][A^{-1}][B] \\
        [D_1] &= [D^{*-1}] 
    """

    def __repr__(self):
        repr_str = "";
        for laminae in self.laminae_list:
            repr_str += "{0}\n".format(laminae)

        return repr_str;
    
    def __init__(self,laminae_list):
        """Initialize class Laminate
        
        Args:
            self: This object instance.
            laminae_list: A list containing objects of type Laminae, 
                defined in Laminae.py
                
        Example:
            ply1 = CompositePlate.Ply(...)
            laminae1 = CompositePlate.Laminae
            laminate = CompositePlate.Laminate([laminae1, laminae2, laminae1]
        """
        self.laminae_list = laminae_list;
        
        # Determine total thickness
        h = [0];
        h_sum = 0;
        for laminae in self.laminae_list:
            h_i = laminae.ply.h;
            h_sum += h_i;
            h.append(h_sum);
        
        _laminae_z = [x-h_sum/2 for x in h]   # Make heights relative to midplane
        self.laminae_midplane_z = [0.5*(_laminae_z[ii+1]+_laminae_z[ii]) for ii in range(len(_laminae_z)-1)]  # Get height to midplane of each laminae
        
        # Calculate Composite matrix
        self.A = np.matrix([[0,0,0],[0,0,0],[0,0,0]]);
        self.B = self.A;
        self.D = self.A;

        for ii in range(0,len(self.laminae_list)):
            laminae = self.laminae_list[ii];
            D_laminae = np.linalg.inv(laminae.S_bar)
            self.A = self.A +D_laminae*(_laminae_z[ii+1] -_laminae_z[ii]);
            self.B = self.B +0.5*D_laminae*(np.power(_laminae_z[ii+1],2) -np.power(_laminae_z[ii],2));
            self.D = self.D +1.0/3.0*D_laminae*(np.power(_laminae_z[ii+1],3) -np.power(_laminae_z[ii],3));
            
        # Calculate Inverse Stiffness matrices
        _A_inv = np.linalg.inv(self.A)
        _D_inv = np.linalg.inv(self.D)
        _Dstar = self.D -self.B*_A_inv*self.B
        
        self.D_1 = np.linalg.inv(_Dstar)
        self.A_1 = _A_inv +_A_inv*self.B*self.D_1*self.B*_A_inv
        self.B_1 = -_A_inv*self.B*self.D_1
        self.C_1 = np.matrix.transpose(self.B_1)
        
        Epsilon = None
        Kappa = None
        
            
    def applied_strain(self, Epsilon, Kappa):
        r"""Method for applying midplane strains to the laminate instance and determine the resulting midplane stresses
        
        Args:
            Epsilon (numpy.matrix): A [3x1] vector of midplane strains.  
                :math:`[\epsilon]_1 = \epsilon_{xx}^0`, :math:`[\epsilon]_2 = \epsilon_{yy}^0`, :math:`[\epsilon]_3 = \gamma_{xy}^0`
            Kappa (numpy.matrix): A [3x1] vector of curvatures.
                :math:`[\kappa]_1 = \kappa_{xx}`, :math:`[\kappa]_2 = \kappa_{yy}`, :math:`[\kappa]_3 = \kappa_{xy}`
                
        Returns:
            numpy.matrix: A dictionary containing the Normal stresses (N) and Moments (M)
            
            N (numpy.matrix): A [3x1] vector containing the normal stresses
            
            M (numpy.matrix): A [3x1] vector containing the resulting moments
            
        """
        N = self.A*Epsilon +self.B*Kappa;
        M = self.B*Epsilon +self.D*Kappa;
        
        return {'N': N,'M': M}
            
    def applied_stress(self,N,M):
        r"""Method for applying forces normalized by the section width to the laminate at the midplane.
        
        Args:
            N (numpy.matrix): A [3x1] vector of normal stresses normalized by the section width
                :math:`[N]_1=N_{xx}`, :math:`[N]_2=N_{yy}`, :math:`[N]_3=N_{xy}`
            M (numpy.matrix): A [3x1] vector of moments normalized by the section width
                :math:`[M]_1=M_{xx}`, :math:`[M]_2=M_{yy}`, :math:`[M]_3=M_{xy}`        
                
        Returns:
            numpy.matrix: A dictionary containing the normal strains (Epsilon) and curvatures (Kappa)
            
            Epsilon (numpy.matrix): A [3x1] vector of laminate midplane strains.  
            
            Kappa (numpy.matrix): A [3x1] vector of laminate curvatures.
        """
        self.Epsilon = self.A_1*N +self.B_1*M;
        self.Kappa = self.C_1*N +self.D_1*M;
        
        return {'Epsilon': self.Epsilon, 'Kappa': self.Kappa}
        
    def laminae_midplane_strain(self, Epsilon=None, Kappa=None):
        r"""Method for determining the laminae midplane strains in 
        response to laminate nidplane strains and curvatures.
            
        If no laminate strain (Epsilon) or curvature (Kappa) is specified,
        the method will attempt to use the laminate's current strain
        and curvaure.  If the laminate's strain and curvature haven't 
        been set and no strain or curvature are specified, the method 
        will throw an :py:class:`InputError` exception.
        
        Note: If the strain and curvature are specified for the method, 
        it does not alter the laminates current state of strain and 
        curvature.  If they are not specified, the laminates strain and 
        curvature state will be used, and the laminae midplane strains 
        will be saved to the laminate state.
        
        Args:
            Epsilon (numpy.matrix): A [3x1] vector of laminate midplane 
                strains.  If a value is specified, this will be assigned
                to the laminate as the current state of strain.  
                
            Kappa (numpy.matrix): A [3x1] vector of laminate curvatures.
                If a value is specified, this will be assigned to the 
                laminate as the current state of curvature.
                
        Returns:
            numpy.matrix: A list of strain values corresponding to the 
                midplane strain in each laminae of the laminate.  
                
                Epsilon (numpy.matrix): A [3x1] vector of laminate midplane strains.  
        
        """
        
        if ((Epsilon is None) or (Kappa is None)) and ((self.Epsilon is None) or (self.Kappa is None)):
            raise InputError('midplane_strain()','Laminate midplane strains and curvatures must either be passed as parameters or defined for the laminate!')
        
        # If Epsilon and Kappa not specified, use laminate values
        if ((Epsilon is None) or (Kappa is None)):
            Epsilon = self.Epsilon
            Kappa = self.Kappa
            
        _laminae_midplane_strains = [];
        for ii, laminae in enumerate(self.laminae_list):
            _laminae_midplane_strains.append(Epsilon +self.laminae_midplane_z[ii]*Kappa)
        
        # If Epsilon and Kappa were not specified, save midplane strains to laminate state.
        if ((Epsilon is None) or (Kappa is None)):
            self.laminae_midplane_strains = _laminae_midplane_strains
        
        return _laminae_midplane_strains;
        
    
    def laminae_stress(self, Epsilon=None, Kappa=None):
        r"""Determine the stresses at the midplane of each laminae in the
        X & Y directions of the laminate.  Also translates these to the 
        principal 1 & 2 direction stresses of the ply.
        
        If no laminate strain (Epsilon) or curvature (Kappa) is specified,
        the method will attempt to use the laminate's current strain
        and curvaure.  If the laminate's strain and curvature haven't 
        been set and no strain or curvature are specified, the method 
        will throw an :py:class:`InputError` exception.
        
        Note: If the strain and curvature are specified for the method, 
        it does not alter the laminates current state of strain and 
        curvature.  If they are not specified, the laminates strain and 
        curvature state will be used, and the laminae stresses will be 
        saved to the laminate state.
        
        Args:
            Epsilon (numpy.matrix): A [3x1] vector of laminate midplane 
                strains.  If a value is specified, this will be assigned
                to the laminate as the current state of strain.  
                
            Kappa (numpy.matrix): A [3x1] vector of laminate curvatures.
                If a value is specified, this will be assigned to the 
                laminate as the current state of curvature.
                
        Returns:
            numpy.matrix: A list of stress values corresponding to the 
                midplane stress in each laminae of the laminate.  
                
                Sigma (numpy.matrix): A [3x1] vector of laminate midplane stresses.  
        """
        if ((Epsilon is None) or (Kappa is None)) and ((self.Epsilon is None) or (self.Kappa is None)):
            raise InputError('midplane_strain()','Laminate midplane strains and curvatures must either be passed as parameters or defined for the laminate!')
        
        # If Epsilon and Kappa not specified, use laminate values
        if ((Epsilon is None) or (Kappa is None)):
            self.Epsilon = Epsilon
            self.Kappa = Kappa
            
        _laminae_midplane_strains = self.laminae_midplane_strain(Epsilon, Kappa)
        
        _laminae_midplane_stresses = [];
        _laminae_principal_stresses = [];
        
        for ii, laminae in enumerate(self.laminae_list):
            _laminae_midplane_stresses.append(laminae.Q_bar*_laminae_midplane_strains[ii])
            _laminae_principal_stresses.append(np.multiply(_laminae_midplane_stresses[ii],rotation_matrix(laminae.theta_rad)))
            # TODO: Check Principal Stresses, Implement Rotation matrix method
        
        return _laminae_midplane_stresses

if __name__ == "__main__":
    import doctest
    doctest.testmod()

########################################################################
#   Module Exceptions        
########################################################################
class CompPlateError(Exception):
    """Base class for exceptions in this module."""
    pass

class InputError(CompPlateError):
    """Exception raised for errors in the input.

    Attributes:
        expr -- input expression in which the error occurred
        msg  -- explanation of the error
    """

    def __init__(self, expr, msg):
        self.expr = expr
        self.msg = msg
