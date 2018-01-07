from TBV_Model import *

HP_input = [
	170, # P1, bypass valve inlet pressure (psia)
	943, # T1, bypass valve inlet temperature (F)
	2338.7, # Cv_max, max Cv at 100% travel
	2165, # Cv_desired, Cv at desired percent travel
	90, # valve_travel, desired percent of max valve travel (%)
	.65, # x_t, pressure drop ratio factor without fittings
	20, # d_inlet, nominal diameter of inlet piping (in)
	28, # d_outlet, nominal diameter of outlet piping (in)
	20, # d_valve, nominal valve diameter (in)
	600000, # m_dot, mass flow rate before attemperation (lbm/hr)
	1225, # h2, target downstream enthalpy after attemperation (Btu/hr)
	125 # T_s, temperature of attemperator spray water (F)
	]

HRH_input = [
	170, # P1, bypass valve inlet pressure (psia)
	943, # T1, bypass valve inlet temperature (F)
	2338.7, # Cv_max, max Cv at 100% travel
	2165, # Cv_desired, Cv at desired percent travel
	90, # valve_travel, desired percent of max valve travel (%)
	.65, # x_t, pressure drop ratio factor without fittings
	20, # d_inlet, nominal diameter of inlet piping (in)
	28, # d_outlet, nominal diameter of outlet piping (in)
	20, # d_valve, nominal valve diameter (in)
	600000, # m_dot, mass flow rate before attemperation (lbm/hr)
	1225, # h2, target downstream enthalpy after attemperation (Btu/hr)
	125, # T_s, temperature of attemperator spray water (F)
	]

LP_input = [
	170, # P1, bypass valve inlet pressure (psia)
	943, # T1, bypass valve inlet temperature (F)
	2338.7, # Cv_max, max Cv at 100% travel
	2165, # Cv_desired, Cv at desired percent travel
	90, # valve_travel, desired percent of max valve travel (%)
	.65, # x_t, pressure drop ratio factor without fittings
	20, # d_inlet, nominal diameter of inlet piping (in)
	28, # d_outlet, nominal diameter of outlet piping (in)
	20, # d_valve, nominal valve diameter (in)
	600000, # m_dot, mass flow rate before attemperation (lbm/hr)
	1225, # h2, target downstream enthalpy after attemperation (Btu/hr)
	125, # T_s, temperature of attemperator spray water (F)
	]

HP_back_input = [
	762975, # m_dot1, mass flow rate for reference case (lbm/hr)
	150, # P11, inlet pressure for reference case (psia)
	409, # T11, inlet temperature for reference case (F)
	650000, # m_dot2, mass flow rate for unknown case (lbm/hr)
	28, # D, pipe diameter (in) - can basically be anything
	0, # P12, condenser sparger back pressure for reference case (psia)
	0, # P22, condenser sparger back pressure for unknown case (psia)
	]

HRH_back_input = [
	762975, # m_dot1, mass flow rate for reference case (lbm/hr)
	150, # P11, inlet pressure for reference case (psia)
	409, # T11, inlet temperature for reference case (F)
	650000, # m_dot2, mass flow rate for unknown case (lbm/hr)
	28, # D, pipe diameter (in) - can basically be anything
	0, # P12, condenser sparger back pressure for reference case (psia)
	0, # P22, condenser sparger back pressure for unknown case (psia)
	]

LP_back_input = [
	762975, # m_dot1, mass flow rate for reference case (lbm/hr)
	150, # P11, inlet pressure for reference case (psia)
	409, # T11, inlet temperature for reference case (F)
	650000, # m_dot2, mass flow rate for unknown case (lbm/hr)
	28, # D, pipe diameter (in) - can basically be anything
	0, # P12, condenser sparger back pressure for reference case (psia)
	0, # P22, condenser sparger back pressure for unknown case (psia)
	]

HP = TBV_model(HP_input, HP_back_input)
HRH = TBV_model(HRH_input, HRH_back_input)
LP = TBV_model(LP_input, LP_back_input)

HRH.plotPvsFlow(200)