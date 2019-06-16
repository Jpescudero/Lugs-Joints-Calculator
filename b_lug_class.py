import numpy as np
import math
from scipy import interpolate


class lug():
	
	def __init__ (self):
		""" Allocating attributes of the class"""

		# Material Properties
		# ---------------------------------------------------- 
		self.Ftu_L = 0
		self.Ftu_T = 0

		self.Fty_L = 0
		self.Fty_T = 0

		# Geometry
		# ----------------------------------------------------

		# Female thickness lug
		self.t01 = 0 
		self.t02 = 0

		# Male thickness lug
		self.ti = 0

		# Gaps between female and male lugs
		self.g1 = 0
		self.g2 = 0

		# Thickness
		self.t = 0
		
		# Lug Diameter 
		self.R0 = 0

		# Hole diameter
		self.D = 0

		# Equivalent lug width
		self.W = 0

		# Distance from the hole center to the end if the lug 
		self.a = 0

		# Tapered angle
		self.beta = 0

		# Input forces
		# ----------------------------------------------------
		self.P_ul = 0
		self.P_ll = 0
		self.theta = 0
		self.FF = 0

		# Solutions
		# ----------------------------------------------------
		#self.Pa = 0
		#self.Ptr = 0
		#self.P01 = 0
		#self.P02 = 0
		#self.Outplane_load = 0
		self.Ptu = 0
		self.Pbru = 0
		self.Py = 0
		self.P_axial_allowable = 0
		self.Ptru = 0
		self.Ptry = 0
		self.P_tr_allowable = 0
		self.Rf_tr = 0
		self.Pu_theta = 0
		self.Rf_theta = 0
		self.Py_theta = 0
		self.Rfy_theta = 0

	def calculate(self):
		""" Run function """

		# Calculate the rest of the geometry
		self.calculate_geometry()

		# Calculate allowables and reserve factors
		self.calculate_allowables()

	def calculate_geometry(self):
		""" Calculation of the rest of the geometry """

		# Convert to radians
		self.beta  = self.beta* np.pi/180
		self.theta = self.theta* np.pi/180

		# Minimum distance between the hole center to the tapered contour of the lug
		#self.R0 = self.W/2 * np.cos(self.beta/2)
		#self.W  = 2*self.R0/np.cos(self.beta/2)
		
		print("W: ",self.W, "D: ",self.D,"R0: ",self.R0)
		print("")

		# Tranverse sectional Areas measured around the hole 
		self.A1 = ((self.W/2)-(2**0.5*self.D/4)*(1-np.tan(self.beta/2)))*self.t
		self.A4 = self.A1

		# Transverse net section Area 
		self.A2 = ((self.W/2)-(self.D/2))*self.t

		# Least area of any radial section around the hole
		self.A3 = (self.R0 -(self.D/2)) *self.t

		# Bearing Area
		self.Abr = self.D * self.t

		# Net section Area
		self.At = (self.W-self.D)* self.t

		print("A1: ",self.A1,"A2: ", self.A2,"A3: ",self.A3,"A4: ",self.A4)
		print("")

		# Average transverse area
		self.Aav = 6/((3/self.A1)+(1/self.A2)+(1/self.A3)+(1/self.A4))
		

	def calculate_allowables (self):

		# In Plane Loads Calculation
		#self.inplane_loads()
		#self.outplane_loads()

		# Load Projections
		self.load_projections()

		# Axial Allowables
		self.tensile_failure()
		self.shear_bearing_failure()
		self.yield_failure_axial()
		self.axial_allowable()

		# Transverse Allowables
		self.transverse_failure()
		self.yield_failure_transverse()
		self.transverse_allowable()

		# Oblique Allowables
		self.oblique_ultimate_allowable()
		self.oblique_yield_allowable()

	
	## In Plane Loads
	## ---------------------------------------------------------------------------------------------

	def inplane_loads (self):

		pass
		#self.P01 = self.P * ((self.ti/2)+(self.t02/2)+self.g2)/((self.t01/2)+(self.t02/2)+self.ti+self.g1+self.g2)
		#self.P02 = self.P - self.P01


	## Out Plane Loads
	## ---------------------------------------------------------------------------------------------

	def outplane_loads (self):

		pass
		#self.outplane_loads = 0.10*self.P


	## Load Projections
	## ---------------------------------------------------------------------------------------------
	def load_projections(self):
		pass
		#self.Pa  = self.P_ul * np.cos(self.theta)
		#self.Ptr = self.P_ul * np.sin(self.theta)


	## Axial Loads Allowable
	## ---------------------------------------------------------------------------------------------

	def tensile_failure(self):

		print("Tensile Failure")
		print("---------------------------------------")
		print("W/D: ", round(self.W/self.D,2))

		W_D   = np.array((1.5 , 2.0, 2.5 , 3.0 , 3.5 , 4.0 , 4.5 , 5))
		Kt_4  = np.array((0.95, 0.9, 0.85, 0.81, 0.78, 0.73, 0.71, 0.69))

		f = interpolate.interp1d(W_D, Kt_4)

		Kt = f(self.W/self.D)

		self.Ptu = Kt * self.Ftu_L * self.At

		print("Kt: ", Kt)
		print("Ptu: ", self.Ptu)
		print("")

	def shear_bearing_failure(self):

		print ("Shear Bearing failure")
		print ("---------------------------------------")
		print ("D/t: ", self.D/self.t)
		print ("a/D: ", self.a/self.D)

		# For values of a/d greater than 2
		a_D  = np.array((1  , 1.2, 1.4, 1.6 , 1.8 , 2))
		Kbr  = np.array((0.82, 1.1, 1.3, 1.55, 1.7 , 1.85))

		f   = interpolate.interp1d(a_D, Kbr)

		Kbr = f(self.a/self.D)

		self.Pbru = Kbr * self.Ftu_L * self.Abr

		
		print("Kbr: ", Kbr)
		print("Pbru: ", self.Pbru)
		print("")

	def yield_failure_axial(self):

		print("Axial Yield Failure")
		print("---------------------------------------")
		print("Pumin_AbrFtumin: ",min([self.Ptu,self.Pbru])/(self.Abr*self.Ftu_L))

		Pumin_AbrFtu = np.array((0.8, 0.9, 1  , 1.5  , 2   , 2.5 , 3 ))
		C  			 = np.array((1.1, 1.1, 1.1, 1.03 , 0.93, 0.82, 0.71 ))

		f   = interpolate.interp1d(Pumin_AbrFtu, C)


		C   = f(min([self.Ptu,self.Pbru])/(self.Abr*self.Ftu_L))

		self.Py = C * (self.Fty_L/self.Ftu_L) * (np.minimum(self.Ptu,self.Pbru))

		
		print("C: ", min([self.Ptu,self.Pbru])/(self.Abr*self.Ftu_L))
		print("Py: ", self.Py)
		print("")


	def axial_allowable (self):

		self.P_axial_allowable = min ([self.Ptu,self.Pbru,1.5*self.Py])
		#self.P_axial_allowable = min ([self.Ptu,self.Pbru])
		#self.RF_a = self.P_axial_allowable/self.Pa


	## Transverse Loads Allowable
	## ---------------------------------------------------------------------------------------------

	def transverse_failure(self):

		print("Transverse Failure")
		print("---------------------------------------")
		print ("Aav/Abr: ", self.Aav/self.Abr)

		Aav_Abr  = np.array((0.1 , 0.2, 0.3  , 0.4  , 0.5 , 0.6 , 0.7 , 0.8 ,0.9 , 1   , 1.4))
		Ktru     = np.array((0.12, 0.23, 0.27, 0.285, 0.3 , 0.31, 0.32, 0.33,0.34, 0.35, 0.35))

		f   = interpolate.interp1d(Aav_Abr, Ktru)

		Ktru   = f(self.Aav/self.Abr)

		

		self.Ptru = Ktru * self.Ftu_T * self.Abr

		print ("Ktru: ", Ktru)
		print ("Ptru: ", self.Ptru)
		print ("")


	def yield_failure_transverse(self):

		print("Transverse Yiel Failure")
		print("---------------------------------------")
		print ("Aav/Abr: ", self.Aav/self.Abr)

		Aav_Abr  = np.array((0.1 , 0.2, 0.3  , 0.4 , 0.5 , 0.6 , 0.7 , 0.8 ,0.9 , 1   , 1.1 , 1.2 , 1.3 , 1.4))
		Ktry     = np.array((0.12, 0.28, 0.38, 0.5 , 0.62, 0.72, 0.82, 0.91,0.99, 1.05, 1.13, 1.19, 1.24, 1.31))

		f   = interpolate.interp1d(Aav_Abr, Ktry)

		Ktry   = f(self.Aav/self.Abr)

		self.Ptry = Ktry * self.Fty_T * self.Abr

		print ("Ktry: ", Ktry)
		print ("Ptry: ", self.Ptry)
		print ("")


	def transverse_allowable (self):

		self.P_tr_allowable = min([self.Ptru,1.5*self.Ptry])

		#self.RF_tr = self.P_tr_allowable/self.Ptr

	## Oblique Loads
	## ---------------------------------------------------------------------------------------------
	
	def oblique_ultimate_allowable (self):

		self.Pu_theta = 1/(np.cos(self.theta)/self.P_axial_allowable**1.6+((np.sin(self.theta)/self.Ptru))**1.6)**0.625

		self.Rf_theta = self.Pu_theta  / (self.FF * self.P_ul)

		print("oblique_ultimate_allowable (Putheta) : ", self.Pu_theta, "RF: ", self.Rf_theta)
		print ("")

	def oblique_yield_allowable (self):

		self.Py_theta = 1/((np.cos(self.theta)/self.Py)**1.6+(np.sin(self.theta)/self.Ptry)**1.6)**0.625

		self.Rfy_theta = self.Py_theta  / (self.P_ll)

		print("oblique_yield_allowable (Pytheta) : ", self.Py_theta, "RF: ", self.Rfy_theta)
		print ("")



class bearing():
	
	def __init__ (self):

		self.d = 0
		self.D = 0
		self.B = 0
		self.C = 0
		self.FF = 0
		self.P_radial = 0
		self.P_Axial  = 0
		self.P_radial_applied = 0
		self.FF = 0

	def calculate(self):

		self.RF = self.P_radial/(self.FF*self.P_radial_applied)
		print ("RF: ", self.RF)
		print ("")




class plain_bushing():
	
	def __init__ (self):

		# Material Properties
		self.Fcy = 0

		# Bush thickness
		self.tb = 0 

		# Lug Females Thickness
		self.tf = 0

		# Inner diameter 
		self.d = 0

		# Applied Load
		self.P_limit_applied = 0

		# Fitting Factor 
		self.FF = 0

	def calculate(self):

		# Bearing Area 
		self.Abr = self.d * self.tf

		# Bush Load allowable 
		self.P_bry = 1.85 * self.Fcy * self.Abr 

		# Reserve Factor
		self.RF = self.P_bry/(self.FF*self.P_limit_applied)

		print ("RF: ", self.RF)
		print ("")



class flange_bushing():
	
	def __init__ (self):

		# Material Properties
		self.Fcy = 0

		# Bush thickness
		self.tb = 0 

		# Lug Females Thickness
		self.tf = 0

		# Inner diameter 
		self.d = 0

		# Applied Load
		self.P_limit_applied = 0

		# Fitting Factor 
		self.FF = 0

	def calculate(self):

		# Bearing Area 
		self.Abr = self.d * self.tf

		# Bush Load allowable 
		self.P_bry = 1.85 * self.Fcy * self.Abr 

		# Reserve Factor
		self.RF = self.P_bry/(self.FF*self.P_limit_applied)

		print ("RF: ", self.RF)
		print ("")

class pin_bolt():

	def __init__ (self):
		

		# Material Properties
		# ---------------------------------------------------- 
		self.Ftu = 0	
		self.Fty = 0
		self.Fsu = 0
		self.Fsy = 0
		self.Fbu = 0 
		self.k   = 0		# plastic bending

		# Geometry
		# ---------------------------------------------------- 
		self.d   = 0		# Pin Diameter
		self.gap = 0		# Gap
		self.ti  = 0 		# Male Lug thickness
		self.tf  = 0		# Female Lug thickness

		# Applicated Force Inputs
		# ---------------------------------------------------- 
		self.P_limit = 0
		self.P_ultimate = 0
		self.FF = 0

	def calculate(self):

		# Pin Area
		self.A  = np.pi * self.d**2 / 4

		# Pin Inertia
		self.I  = np.pi * (self.d/2)**4 / 4

		# Reduction factor 
		gamma = 1

		# Arm
		self.b  = self.tf/2 + self.gap + gamma*self.ti/4

		# Calculate Reserve factors
		self.allowable("ultimate")
		self.allowable("limit")


	def allowable(self,case):


		# Bending Moment 
		M = getattr(self,"P_"+case) * self.b/2

		print(case,M)
		# Bending Stress
		sigma = M * self.d / 2 / self.I


		# Reserve Factor Bending
		if case == "ultimate":
			setattr(self, "RFbn_" + case, self.Fbu/sigma)
		else: 
			setattr(self, "RFbn_" + case, self.Fty/sigma)

		#shear stress 
		tau = getattr(self,"P_" + case)/self.A/2

		# Reserve Factor shear 
		if case == "ultimate":
			setattr(self, "RFsh_" + case, self.Fsu/tau)
		else:
			setattr(self, "RFsh_" + case, self.Fsy/tau)

		# Bending shear combination
		setattr(self, "RFpin_" + case, 1/(self.FF* ((1/getattr(self, "RFbn_" + case))**2 + (1/getattr(self, "RFsh_" + case))**2)**0.5))

