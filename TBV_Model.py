'''
Turbine Bypass Valve Model v5 - January 6, 2018
Justin Gomes

UPDATES: 
    - Updated for Python 3
    - Switched from freesteam to iapws Python module

Model simulates bypass to condenser for HRH and LP to determine required Cv, maximum flow rate, or 
pressure required.  Also calculates values such as attemperator spray flow rate, total flow rate, and 
bypass valve back pressure.  Can be used to model HP bypass to CRH if user manually enters back pressure.

Seems to differ from vendor Cv values because of different calculations for:
	- ratio of specific heats
	- steam density at TBV inlet
	- x = (P1 - P2)/P1

ASSUMPTIONS: 
	- Perfect vacuum in condenser
	- Joule-Thompson expansion through TBV
	- All processes are adiabatic

'''

import matplotlib.pyplot as plt
from pylab import *
from iapws import IAPWS97
from math import pi, log10, sqrt

def float_range(start, end, interval):
    '''
    Creates list of floats in "interval" increments
	
    input: 
        start: FLOAT or INT, first number of the list
        end: FLOAT or INT, last number to be in the list
        interval: FLOAT, fraction of difference between "start" and "end" over which to iterate
		
    output: 
        list: LIST of floats from "start" to "end"
    '''
    l = range(int((end - start)*1.0/interval) + 1)
    l2 = [start]
    for n in l[1:]:
        l2.append(round(interval + l2[n-1], 7))
    return l2

class back_P(object):
    ''' 
    Create instance of a simulated condenser sparger or other downstream fixture.
    Determine back pressure felt by TBV assuming:
        - default: condenser sparger reduces pressure to 0
        - optional inputs for P12 and P22, downstream pressures for non-zero inputs (i.e. HP bypass to CRH)
        - any downstream fixture can be treated as a pipe
        - pressure drop through fixture is adiabatic
	
    All inputs in imperial units - functions convert to SI as required
    '''
	
    def __init__(self, back_input):
        self.m_dot1 = back_input[0]*0.000125998 # mass flow rate for reference conditions in kg/s
        self.P11 = back_input[1]*6894.76/1000000 # upstream pressure for reference conditions in MPa
        self.T11 = (back_input[2]+491.67)*0.5555556 # upstream temperature for reference conditions in K
        self.m_dot2 = back_input[3]*0.000125998 # upstream flow rate for unknown conditions in kg/s
        self.D = back_input[4]*.0254 # pipe diameter in m
        self.P12 = back_input[5]*6894.76/1000000 # downstream pressure for reference case in MPa
        self.P22 = back_input[6]*6894.76/1000000 # downstream pressure for unknown case in MPa
	
    def calc_fric_coeff(self,rho,V,mu):
        '''
        calc_fric_coeff() uses bisection search to iterate through potential solutions to the 
        Colebrook equation to determine the friction coefficient.
		
        Due to the fact that the friction coefficient cannot be determined for transient flow
        conditions, the function treats transient flow as turbulent, adding in some error.
		
        INPUTS:
            - rho: density (kg/m^3)
            - V: flow velocity (m/s)
            - mu: viscosity
		
        OUTPUT: 
            - fd: friction coefficient
        '''
        Re = self.D*V*rho/mu # pipe Reynolds number in m/s
        if Re <= 2300: # laminar flow
            fd = 64/Re
        elif Re > 2300 and Re <4000: # transient flow conditions
#            raise ValueError, "transient flow - friction coefficient cannot be calculated"
            '''
            For the sake of keeping the simulation running, pretend transient flow is turbulent
            '''
            k = 0.05*10**-3 # absolute roughness in m - generic value for commercial steel pipe
			
            # Use bisection search to iterate through potential solutions for the Colebrook 
            # equation to find the Darcy Friction Factor, fd.
            fd_high = 0.1
            fd_low = 0
            fd_guess = (fd_high + fd_low)/2.
            fd = (1/(-2*log10(k/(3.7*self.D)+2.51/(Re*sqrt(fd_guess)))))**2

            while abs(fd-fd_guess) > fd_guess/100.:
                if fd > fd_guess:
                    fd_low = fd_guess
                    fd_guess = (fd_high + fd_low)/2.
                else:
                    fd_high = fd_guess
                    fd_guess = (fd_high + fd_low)/2.
				
            fd = (1/(-2*log10(k/(3.7*self.D)+2.51/(Re*sqrt(fd_guess)))))**2
			
        elif Re >= 4000: # turbulent flow
            k = 0.05*10**-3 # absolute roughness in m - generic value for commercial steel pipe
			
            # Use bisection search to iterate through potential solutions for the Colebrook 
            # equation to find the Darcy Friction Factor, fd.
            fd_high = 0.1
            fd_low = 0
            fd_guess = (fd_high + fd_low)/2.
            fd = (1/(-2.*log10(k/(3.7*self.D)+2.51/(Re*sqrt(fd_guess)))))**2

            while abs(fd-fd_guess) > fd_guess/100.:
                if fd > fd_guess:
                    fd_low = fd_guess
                    fd_guess = (fd_high + fd_low)/2.
                else:
                    fd_high = fd_guess
                    fd_guess = (fd_high + fd_low)/2.
				
            fd = (1/(-2.*log10(k/(3.7*self.D)+2.51/(Re*sqrt(fd_guess)))))**2
		
        return fd
	
    def calc_vel(self,m_dot,D,v):
        V = m_dot*v/((D/2.)**2*pi) # pipe flow velocity in m/s
        return V
	
    def calc_length(self, fd, rho, V, delta_P, D):
        length = delta_P/(fd*rho/2.*V**2/D)
        return length
	
    def calc_P_drop(self, length, fd, rho, V, D):
        P_drop = length*fd*rho/2.*V**2/D
        return P_drop
	
    def solve_back_P(self, P_high, P_low, T21, length):
        P_guess = (P_high + P_low) / 2. # initial pressure guess
        S2 = IAPWS97(P=P_guess, T=T21) # create instance for guessed steam conditions
        V2 = self.calc_vel(self.m_dot2, self.D, S2.v)
        fd2 = self.calc_fric_coeff(S2.rho, V2, S2.mu)
        P_drop = self.calc_P_drop(length, fd2, S2.rho, V2, self.D)
		
        # iterate using bisection search until P_guess matches P_drop
        while abs(P_guess - P_drop) >= P_guess/100.:
            if P_guess > P_drop:
                P_high = P_guess
                P_guess = (P_high + P_low) / 2.
            elif P_guess < P_drop:
                P_low = P_guess
                P_guess = (P_high + P_low) / 2.
            try:
                print('P = %sMpa and T = %sK' %(P_guess, T21))
                S2 = IAPWS97(P=P_guess, T=T21) # create instance for guessed steam conditions
            except ValueError: print('Pressure or temperature out of bounds')
            
            V2 = self.calc_vel(self.m_dot2, self.D, S2.v)
            fd2 = self.calc_fric_coeff(S2.rho, V2, S2.mu)
            P_drop = self.calc_P_drop(length, fd2, S2.rho, V2, self.D)
			
        return P_guess
	
    def calc_back_P(self, P11 = None):
        '''
        calc_P_drop guesses a pressure upstream of the condenser and finds the corresponding density.
        It uses the density and an assumed backpressure of 0 to calculate a pressure drop through the 
        sparger.  The function iterates until the guess and the pressure drop match.
		
        All calculations in SI units.
		
        returns P21 - pressure upstream of condenser/back pressure for TBV in psia
        '''
        if P11 == None: P11 = self.P11
        T21 = self.T11 # upstream temperature for unknown conditions in K
        S = IAPWS97(P=P11, T=self.T11) # create instance for reference steam conditions
        V = self.calc_vel(self.m_dot1, self.D, S.v) # call function to return reference flow upstream velocity in m/s
        fd = self.calc_fric_coeff(S.rho, V, S.mu) # call function to return reference flow upstream Reynolds number
		
        # equivalent length of pipe to represent pressure drop through condenser sparger
        length = self.calc_length(fd, S.rho, V, self.P11-self.P12, self.D)
		
        if self.m_dot2 < self.m_dot1:
            # set initial guess bounds
            P_high = P11
            P_low = 0
            P_guess = self.solve_back_P(P_high, P_low, T21, length)
        elif self.m_dot2 == self.m_dot1: pass
        elif self.m_dot2 > self.m_dot1:
            # set initial guess bounds - P_high is an arbitrarily high number to avoid exceeding the upper bound
            P_high = 6894.76/100
            P_low = P11
            P_guess = self.solve_back_P(P_high, P_low, T21, length)

        return (P_guess+self.P22)/(6894.76/1000000) # return guess that gives matching P_drop in psia

class TBV_model(back_P):
	
    def __init__(self, TBV_input, back_input):
        self.TBV_input = TBV_input
        self.back_input = back_input
        self.P1 = TBV_input[0]
        self.P1_SI = self.P1*6894.76/1000000 # in MPa
        self.T1 = TBV_input[1]
        self.T1_SI = (self.T1+491.67)*0.5555556 # in K
        self.Cv_max = TBV_input[2]
        self.Cv_desired = TBV_input[3]
        self.valve_travel = TBV_input[4]
        self.x_t = TBV_input[5]
        self.d_inlet = TBV_input[6]
        self.d_outlet = TBV_input[7]
        self.d_valve = TBV_input[8]
        self.m_dot = TBV_input[9]
        self.h2 = TBV_input[10]
        self.T_s = TBV_input[11]
        self.back_input = back_input
        self.back = back_P(self.back_input)
	
    def calc_sigmaK(self): # calculate velocity head loss coefficient
        K1 = 0.5*(1-self.d_valve**2./self.d_inlet**2)**2
        K2 = (1-self.d_valve**2./self.d_outlet**2)**2
        KB1 = 1 - (self.d_valve/self.d_inlet)**4
        KB2 = 1 - (self.d_valve/self.d_outlet)**4
        sigmaK = K1+K2+KB1-KB2
        return sigmaK
	
    def calc_Fp(self, sigmaK): # calculate piping geometry factor
        Fp = 1./(sqrt(1+sigmaK/890.*(self.Cv_max/self.d_valve**2)**2))
        return Fp
	
    def calc_x_tp(self, Fp): # calculate pressure drop ratio factor
        K1 = 0.5*(1-self.d_valve**2./self.d_inlet**2)**2
        KB1 = 1 - (self.d_valve/self.d_inlet)**4
        Ki = K1 + KB1
        x_tp = self.x_t/Fp**2./(1+(self.x_t*Ki/1000.)*(self.Cv_max/self.d_valve**2)**2)
        return x_tp
	
    def calc_Y(self, x, x_tp, Fk): # calculate Y expansion factor
        Y = 1 - x/(3.*Fk*x_tp)
        return Y
	
    def calc_Cv(self, P1, P2, m_dot): # calculate Cv for given pressure drop and flow rate.  raises exception if backpressure is entered incorrectly
        '''
        P1, P2, and m_dot in english units
        '''
        sigmaK = self.calc_sigmaK()
        Fp = self.calc_Fp(sigmaK)
        x = (P1 - P2)/P1 # vendor may provide different value that may result in results more consistent with the vendor's
        x_tp = self.calc_x_tp(Fp)
        stm = IAPWS97(P=P1*6894.76/1000000, T=self.T1_SI)
#        Fk = 1.33/1.4 # vendor may provide different value that may result in results more consistent with the vendor's - 1.33 is a mid range value of cp/cv for superheated steam
        Fk = stm.cp/stm.cv/1.4 # more accurate calculation of Fk
        Y = self.calc_Y(x, x_tp, Fk)
        rho = stm.rho*.062428 # vendor may provide different value that may result in results more consistent with the vendor's
		
        try:
            Cv = m_dot/(63.3*Fp*Y*sqrt(P1*rho*x))
        except ValueError: print("Backpressure exceeds upstream pressure.")
		
        return Cv
	
    def Cv_required(self, P1 = None, P2 = None): # Calculate Cv required to pass the given flow rate at the given steam conditions
        '''
		
        '''
        if P1 == None: P1 = self.P1
        if P2 == None: P2 = self.back.calc_back_P()
        Cv_req = self.calc_Cv(P1, P2, self.m_dot)
		
        return Cv_req
	
    def calc_P2(self, P1, m_dot): # Determine outlet pressure from TBV for given flow rate and steam conditions
        '''
        Given the valve described by the given instance of TBV_capacity, a chosen inlet 
        pressure, and a mass flow rate, calc_P2 uses bisection search to iterate through
        potential outlet pressures until the calculated Cv for the pressure drop converges
        on the actual Cv defined by the instance.
		
        Uses english units
		
        If valve capacity is exceeded, function will not converge and will return a ValueError
        '''
        P_high = P1
        P_low = 0
        P_guess = (P_high + P_low)/2.
		
        Cv_guess = self.calc_Cv(P1, P_guess, m_dot)
		
        n = 0
        while abs(Cv_guess - self.Cv_desired) > self.Cv_desired/100.:
            n += 1 # count iterations
            if Cv_guess > self.Cv_desired:
                P_high = P_guess
                P_guess = (P_high + P_low)/2.
            elif Cv_guess < self.Cv_desired:
                P_low = P_guess
                P_guess = (P_high + P_low)/2.
            Cv_guess = self.calc_Cv(P1, P_guess, m_dot)
            if n > 1000: # escape - halts runaway function if valve capacity is exceeded
                raise ValueError('Function will not converge.  Valve capacity exceeded.')
        return P_guess
	
    def calc_m_dot_total(self,P1 = None, m_dot = None):
        '''
        Calculate spray water mass flow rate and total mass flow rate after attemperation.
		
        Inputs and output are in english units, but function converts to SI for all calculations.
        '''
        if P1 == None: P1 = self.P1
        if m_dot == None: m_dot = self.m_dot
        P1 = P1*6894.76/1000000 # convert P1 to metric
        stm1 = IAPWS97(P=P1, T=self.T1_SI)
        stm_s = IAPWS97(T=(self.T_s+491.67)*0.5555556,x=0)
        h2 = self.h2*1000/0.429923 # convert desired enthalpy to J/kg
        # Assume Joule-Thompson expansion through TBV - no change in enthalpy
        m_dot_s = m_dot*0.000125998*(h2-stm1.h)/(stm_s.h-h2) # calculate spray water flow rate (kg/s)
        m_dot_total = m_dot_s + m_dot*0.000125998 # Total flow rate after attemperation in kg/s
		
        return m_dot_total/0.000125998 # return total flow rate in lbm/hr
	
    def calc_P1_required(self):
        P_high = 5000. # uses arbitrarily high upper limit for pressure to ensure solution can be found
        P_low = 0
        P_guess = (P_high+P_low)/2.
		
        Cv_guess = self.Cv_required(P_guess)
		
        while abs(Cv_guess - self.Cv_desired) > self.Cv_desired/1000.:
            if Cv_guess > self.Cv_desired:
                P_low = P_guess
                P_guess = (P_high+P_low)/2.
            elif Cv_guess < self.Cv_desired:
                P_high = P_guess
                P_guess = (P_high+P_low)/2.
            Cv_guess = self.Cv_required(P_guess)
		
        return P_guess
	
    def calc_flow_rate_1(self, p, m_dot_guess, m_dot_low = 0):
        '''
        Picks a flow rate guess by taking an exception if flow rate exceeds valve capacity
        and recursively calling itself to pick a new flow rate until it finds one that does
        not exceed capacity.
		
        Inputs are in english units
        '''

        try: P2 = self.calc_P2(p, m_dot_guess)
        except ValueError:
            m_dot_high = m_dot_guess
            m_dot_guess = (m_dot_high + m_dot_low)/2.
            m_dot, P2 = self.calc_flow_rate_1(p, m_dot_guess, m_dot_low)
            return m_dot, P2
        return m_dot_guess, P2
	
    def calc_flow_rate(self, p, m_dot_high = 1000000.):
        '''
        Calculate flow rate for given pressure starting at a high flow rate guess.  It starts with
        an arbitrarily high upper limit for the flow rate unless another input is provided.
		
        Inputs are in english units
        '''
		
        m_dot_low = 0
        m_dot_guess, P2 = self.calc_flow_rate_1(p,m_dot_high)
		
        m_dot_total = self.calc_m_dot_total(p, m_dot_guess)
        back = self.back_input[:]
        back[3] = m_dot_total
        backTest = back_P(back)
        P_test = backTest.calc_back_P()
        n = 0
        while abs(P_test-P2) > p/100.: # converge back pressure calcs to within 1% of upstream pressure
            n += 1
            if n > 100 and abs(P_test-P2) < p/10.: break # escape - breaks loop if iterations exceed 100, ensures accuracy within 10%
            elif n > 100 and abs(P_test-P2) >= p/10.: 
                raise ValueError("TBV outlet pressure will not converge")
            if P_test < P2: 
                m_dot_low = m_dot_guess
                m_dot_guess = (m_dot_high + m_dot_low)/2.
                m_dot_guess, P2 = self.calc_flow_rate_1(p,m_dot_guess, m_dot_low)
            elif P_test > P2:
                m_dot_high = m_dot_guess
                m_dot_guess = (m_dot_high + m_dot_low)/2.
                m_dot_guess, P2 = self.calc_flow_rate_1(p,m_dot_guess, m_dot_low)
			
            m_dot_total = self.calc_m_dot_total(p, m_dot_guess)
            back = self.back_input[:]
            back[3] = m_dot_total
            backTest = back_P(back)
            P_test = backTest.calc_back_P()
		
        return m_dot_guess
	
    def plotPvsFlow(self,maxP):
        m_dot_list = []
        P_list = []
		
        for p in float_range(100, maxP, maxP/1000.):
            m_dot = self.calc_flow_rate(p) # uses default max flow rate upper bound
            print(p) # only here to indicate to the user that the function is still working
            m_dot_list.append(m_dot)
            P_list.append(p)
		
        plt.figure(1)
        P_vs_mdot = plt.plot(P_list, m_dot_list)
		
        plt.setp(P_vs_mdot, color='k', lw=1.0)
        plt.xlabel('Upstream Pressure (psia)')
        plt.ylabel('Mass Flow Rate (lbm/hr)')
        plt.grid(True)

