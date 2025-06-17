"""
@author: Jack Charles   jack@jackcharlesconsulting.com
"""

#*************************************************************************
#pi = 3.14159265358979
GRAVITY_CONSTANT = 32.174    #gravitational constant, ft/s**2
#to check all _CONST used, input only the units into MathCAD and solve for desired output units
#example, for calcVelocity input bbl/min * 1/(in**2) = _CONST ft/s
#*************************************************************************
#VELOCITY, RHEOLOGY, FRICTION EQUATIONS
#DELTA-PRESSURE EQUATIONS
#SETTLING AND TRANSPORT EQUATIONS
#GRAVEL PACK SLURRY EQUATIONS
#SANDOUT EQUATIONS
#FRACTURE MECHANICS EQUATIONS
#ALPHA BETA PACKING
#SURVEY EQUATIONS
#TUBULAR STRESS ANALYSIS EQUATIONS
#PRODUCTION EQUATIONS
#PVT CALCULATIONS
#OTHER EQUATIONS
#UNIT CONVERSION FUNCTION
#OTHER MATH FUNCTIONS
#*************************************************************************

import math
import numpy as np
import matplotlib.pyplot as plt
import thermopy as tp

#*************************************************************************
#VELOCITY, RHEOLOGY, FRICTION EQUATIONS
#*************************************************************************
def calc_fluid_velocity(fluid_rate, casing_diameter, tubing_diameter=0):
    #flow_rate:bbl/min   tubing_diameter, casing_diameter:in   output:ft/s
    #this def only needs one diameter input, other can be 0. No preference on order
    _CONST = 13.475
    if abs(casing_diameter ** 2 - tubing_diameter ** 2) == 0: 
        fluid_velocity = 0
    else: 
        fluid_velocity = _CONST * 4 * fluid_rate / math.pi / abs(casing_diameter ** 2 - tubing_diameter ** 2)
    return fluid_velocity

def calc_shear_rate(nPrime, fluid_velocity, hydraulic_diameter, geometry: str):
    #fluid_velocity:ft/s    hydraulic_diameter:in   geometry:Pipe or Slot   output:1/s
    _CONST = 12
    if geometry == "Pipe": 
        shear_rate = ((3 * nPrime + 1) / (4 * nPrime)) * _CONST * 8 * fluid_velocity / hydraulic_diameter
    elif geometry == "Slot": 
        shear_rate == ((2 * nPrime + 1) / (3 * nPrime)) * _CONST * 6 * fluid_velocity / hydraulic_diameter
    else: 
        shear_rate = 0
    return shear_rate

def calc_power_viscosity_pipe(nPrime, kPrime, fluid_velocity, hydraulic_diameter):
    #kPrime:lbf-sn/ft2, generalized K, not geometry specific!    fluid_velocity:ft/s   hydraulic_diameter:in     output:cP
    _CONST = 12
    shear_rate = ((3 * nPrime + 1) / (4 * nPrime)) * _CONST * 8 * fluid_velocity / hydraulic_diameter
    fluid_viscosity = (kPrime * ((3 * nPrime + 1) / (4 * nPrime)) ** nPrime) * shear_rate ** (nPrime - 1)
    power_viscosity_pipe = unit(fluid_viscosity, "lbf-s/ft2", "cP")
    return power_viscosity_pipe

def calc_power_viscosity_slot(nPrime, kPrime, fluid_velocity, hydraulic_diameter):
    #kPrime:lbf-sn/ft2, generalized K, not geometry specific!    fluid_velocity:ft/s   hydraulic_diameter(slot width):in     output:cP
    _CONST = 12
    shear_rate = ((2 * nPrime + 1) / (3 * nPrime)) * _CONST * 6 * fluid_velocity / hydraulic_diameter
    fluid_viscosity = (kPrime * ((2 * nPrime + 1) / (3 * nPrime)) ** nPrime) * shear_rate ** (nPrime - 1)
    power_viscosity_slot = unit(fluid_viscosity, "lbf-s/ft2", "cP")
    return power_viscosity_slot

def calc_power_viscosity_apparent(nPrime, kPrime, shear_rate):      
    #kPrime:lbf-sn/ft2, generalized K, not geometry specific!    shearRate:1/s     output:cP
    #use with bob data
    _CONST = 1
    #fluid_viscosity = 47880 * kPrime / shearRate ** (1 - nPrime)   #no conversion needed, kPrime is geometry-dependent
    fluid_viscosity = _CONST * kPrime * ((3 * nPrime + 1) / (4 * nPrime)) ** nPrime * shear_rate ** (nPrime - 1)
    power_viscosity_apparent = unit(fluid_viscosity, "lbf-s/ft2", "cP")
    return power_viscosity_apparent

def calc_kPrime(nPrime, fluid_viscosity, shear_rate):
    #fluid_viscosity:cP    shear_rate:1/s     output:lbf-sn/ft2
    _CONST = 1
    kPrime = _CONST * unit(fluid_viscosity, "cP", "lbf-s/ft2") / (((3 * nPrime + 1) / (4 * nPrime)) ** nPrime * shear_rate ** (nPrime - 1))
    return kPrime

def calc_slurry_viscosity_keck(nPrime, fluid_velocity, hydraulic_diameter, fluid_density, solid_density, solid_loading):
    #fluid_velocity:ft/s   hydraulic_diameter:in    solid_density,fluid_density,solid_loading:ppg(a)  output:cP
    #SPE 19771 & SPE 104253 - multiply result by base fluid viscosity
    c = solid_loading / (solid_density + solid_loading)
    shear_rate = 12 * 8 * fluid_velocity / hydraulic_diameter
    relative_viscosity = (1 + (0.75 * (math.exp(1.5 * nPrime) - 1) * math.exp(-(1 - nPrime) * shear_rate / 1000)) * (1.25 * c / (1 - 1.5 * c))) ** 2
    relative_density = (1 + solid_loading / fluid_density) / (1 + solid_loading / solid_density)
    slurry_viscosity_keck = relative_viscosity ** 0.55 * relative_density ** 0.45
    
    #slurry_viscosity_keck = (1 + solid_loading / fluid_density) / (1 + solid_loading / solid_density)
    
    #slurry_viscosity_mult = relative_viscosity ** 0.6 * relative_density ** 0.4    #laminar
    #slurry_viscosity_mult = relative_viscosity ** 0.5 * relative_density ** 0.5    #transition
    #slurry_viscosity_mult = relative_viscosity ** 0.2 * relative_density ** 0.8    #turbulent
    return slurry_viscosity_keck

def calc_slurry_viscosity(solid_loading, solid_density, fluid_density):  
    #solid_density,fluid_density,solid_loading:ppg   hydraulic_diameter:in   output:cP
    #generalized collection. multiply result by base fluid viscosity
    c = solid_loading / (solid_density + solid_loading)
    #relative_viscosity = 1 + 2.5 * c                                    #Einstein
    #relative_viscosity = math.exp(2.7 * c / (1 - S * c))                #AT, S = 1.5-1.7
    relative_viscosity = (1 + 2.61 * (1.25 * c) / (1 - 1.5 * c)) ** 2    #AT
    
    #cmax = 0.58     #42% porosity
    #relative_viscosity = 1 / (1 - c / cmax) ** (2.5 * nPrime)           #Nolte 1988
    #relative_viscosity = 1 / (1 - c / cmax) ** (2.5 * cmax)             #Krieger-Dogherty (Hackley and Ferraris, 2001)
    return relative_viscosity

def calc_NRe_newton(fluid_velocity, hydraulic_diameter, fluid_density, fluid_viscosity):
    #fluid_velocity:ft/s  hydraulic_diameter:in     fluid_density:ppg     fluid_viscosity:cP
    _CONST = 927.6866
    NRe_newton = _CONST * fluid_velocity * hydraulic_diameter * fluid_density / fluid_viscosity
    return NRe_newton

def calc_NRe_power(fluid_velocity, hydraulic_diameter, fluid_density, nPrime, kPrime):
    #fluid_velocity:ft/s  hydraulic_diameter:in     fluid_density:ppg     kPrime:lbf-s**n/ft2
    #KPrime here is from concentric cylinder test. it is not the same as the generalized K
    _CONST = 7.48052 / 32.2 / 12 ** nPrime
    NRe_power = _CONST * hydraulic_diameter ** nPrime * fluid_velocity ** (2 - nPrime) * fluid_density / ((kPrime * 8 ** (nPrime - 1)) * ((3 * nPrime + 1) / (4 * nPrime)) ** nPrime)
    return NRe_power

def calc_friction_chen(hydraulic_diameter, NRe, roughness):
    #hydraulic_diameter, roughness:in
    #Fanning friction factor, Chen explicit solution
    #FrictionChen = 1 / 4 / (-2 * log10(e / 3.7065 - 5.0452 / NRe * log10(e ** 1.1098 / 2.8257 + 5.8506 / (NRe ** 0.8981)))) ** 2
    if NRe < 2100:
        friction_chen = 16 / NRe
    else:
        _e = roughness / hydraulic_diameter
        friction_chen = 1 / (-4 * math.log10(_e / 3.7065 - 5.0452 / NRe * math.log10(_e ** 1.1098 / 2.8257 + (7.149 / NRe) ** 0.8981))) ** 2
    return friction_chen

def calc_friction_colebrook(hydraulic_diameter, NRe, roughness):    
    #hydraulic_diameter, roughness:in
    #Fanning friction factor, Colebrook-White implicit solution
    #FrictionColebrook = 1 / (f ** 0.5) + 4 * log10(e / 3.7065 + 1.2613 / (NRe * (f ** 0.5)))
    if NRe < 2100:
        friction_colebrook = 16 / NRe
    else:
        _e = roughness / hydraulic_diameter
        #Newton-Raphson solution
        f1 = 0.00001
        f01 = 1
        dF = 0.000005
        while f01 > 0.000000001:
            f0 = f1
            #ff1 = 1 / (f0 ** 0.5) + 2 * log10(_e / 3.7065 + 2.51 / (NRe * (f0 ** 0.5)))               #Darcy friction factor
            #ff2 = 1 / ((f0 + dF) ** 0.5) + 2 * log10(_e / 3.7065 + 2.51 / (NRe * ((f0 + dF) ** 0.5)))
            ff1 = 1 / (f0 ** 0.5) + 4 * math.log10(_e / 3.7065 + 1.2613 / (NRe * (f0 ** 0.5)))              #Fanning friction factor
            ff2 = 1 / ((f0 + dF) ** 0.5) + 4 * math.log10(_e / 3.7065 + 1.2613 / (NRe * ((f0 + dF) ** 0.5)))
            f1 = f0 - ff1 / ((ff2 - ff1) / dF)
            f01 = abs(f1 - f0)
        friction_colebrook = f1
    return friction_colebrook

def calc_friction_power_explicit(NRe, nPrime):         
    #Fanning friction factor, power law from Economides
    if NRe < 2100:
        friction_power_explicit = 16 / NRe
    else:
        _B = (1.4 - math.log10(nPrime)) / 7
        c = (math.log10(nPrime) + 2.5) / 50
        friction_power_explicit = c / NRe ** _B
    return friction_power_explicit
    
def calc_friction_power_implicit(NRe, nPrime):        
    #Fanning friction factor, power law from SPE handbook
    #calc_friction_power_implicit = 1 / f ** 0.5 - ((4 / nPrime ** 0.75) * log10(NRe * f ** (1 - nPrime / 2)) - 0.4 / nPrime ** 1.2)
    if NRe < (3250 - 1150 * nPrime):
        friction_power_implicit = 16 / NRe
    else:
    #Newton-Raphson solution
        f1 = 0.00001
        f01 = 1
        dF = 0.000005
        while f01 > 0.000000001:
            f0 = f1
            ff1 = 1 / f0 ** 0.5 - ((4 / nPrime ** 0.75) * math.log10(NRe * f0 ** (1 - nPrime / 2)) - 0.4 / nPrime ** 1.2)
            ff2 = 1 / (f0 + dF) ** 0.5 - ((4 / nPrime ** 0.75) * math.log10(NRe * (f0 + dF) ** (1 - nPrime / 2)) - 0.4 / nPrime ** 1.2)
            f1 = f0 - ff1 / ((ff2 - ff1) / dF)
            f01 = abs(f1 - f0)
        
    friction_power_implicit = f1
    return friction_power_implicit
    
def calc_standoff(tubing_diameter, casing_diameter, centralizer_diameter, blade_count):
    angle = math.pi / blade_count
    standoff1 = (centralizer_diameter * math.cos(angle) - tubing_diameter) / 2
    #both standoff2 are correct, just different forms. http://mathworld.wolfram.com/CircularSegment.html
    #standoff2 = (2 * casing_diameter / 2 - ((2 * casing_diameter / 2) ** 2 - 4 * (centralizer_diameter / 2 * sin(angle)) ** 2) ** 0.5) / 2
    standoff2 = casing_diameter / 2 - ((casing_diameter / 2) ** 2 - (centralizer_diameter / 2 * math.sin(angle)) ** 2) ** 0.5
    standoff = standoff1 + standoff2
    return standoff

def calc_eccentricity(tubing_diameter, casing_diameter, standoff_radius):
    #range from 0 (fully concentric) to 1 (fully eccentric)
    #standoff_radius is used instead of centralizer_diameter because the standoff on a bladed centralizer will be less. See calc_standoff
    if standoff_radius < 0:
        eccentricity = 1
    else: 
        eccentricity = (casing_diameter / 2 - (standoff_radius + tubing_diameter / 2)) / (casing_diameter / 2 - tubing_diameter / 2)
        #eccentricity = (casing_diameter - (centralizer_diameter) / (casing_diameter - tubing_diameter)
    return eccentricity

def calc_eccentricity_factor(tubing_diameter, casing_diameter, eccentricity):  
    #SPE 31147
    #This is multiplied by Dh, per paper Deff = S*Dh. Two different solutions for eccentricity are provided in the versions of the papers
    if eccentricity == 0:
        eccentricity_factor = 0.66667
    else:
        #eccentricity_factor = 4.33 + (tubing_diameter / casing_diameter) * (10.54 * tubing_diameter / casing_diameter - 12.26)
        eccentricity_factor = 1.004 + 1.873 * (tubing_diameter / casing_diameter) - 1.826 * (tubing_diameter / casing_diameter) ** 2 + 0.617 * (tubing_diameter / casing_diameter) ** 3
    return eccentricity_factor
    
def calc_eccentricity_factor_powerlaw(nPrime, kPrime, NRe, tubing_diameter, casing_diameter, eccentricity):  
    #this is multiplied by f or dP
    diameter_ratio = tubing_diameter / casing_diameter
    if nPrime > 0:  #check for error
        if NRe < (3250 - 1150 * nPrime):
            #Haciislamoglu & Langlinais, 1990 - as described in SPE 111514
            eccentricity_factor_powerlaw = 1 - (0.072 * eccentricity / nPrime * diameter_ratio ** 0.8454) - (1.5 * eccentricity ** 2 * math.sqrt(nPrime) * diameter_ratio ** 0.1852) + (0.96 * eccentricity ** 3 * math.sqrt(nPrime) * diameter_ratio ** 0.2527)
        else:
            #Aadnoy, Cooper, Miska, Mitchel, & Payne, 2009 
            eccentricity_factor_powerlaw = 1 - (0.048 * eccentricity / nPrime * kPrime ** 0.8454) - (0.67 * eccentricity ** 2 * math.sqrt(nPrime) * kPrime ** 0.1852) + (0.28 * eccentricity ** 3 * math.sqrt(nPrime) * kPrime ** 0.2527) 
    else: 
        eccentricity_factor_powerlaw = 1
    return eccentricity_factor_powerlaw
    
    
    
    
    
#*************************************************************************
#DELTA-PRESSURE EQUATIONS
#*************************************************************************
def calc_DPpe(fluid_density, length, angle):
    #fluid_density:ppg    angle:deg   length:ft   output:psi
    #pressure drop due to potential energy change
    _CONST = 0.001615
    DPpe = GRAVITY_CONSTANT * _CONST * fluid_density * length * math.cos(angle * math.pi / 180)
    #_CONST = 0.052
    #DPpe = _CONST * fluid_density * length * math.cos(angle * math.pi / 180)
    return DPpe

def calc_DPke(fluid_density, fluid_velocity_inlet, fluid_velocity_outlet, Cd):    
    #fluid_density:ppg    fluid_velocity_inlet,fluid_velocity_utlet:ft/s    output:psi
    #pressure drop due to kinetic energy change
    _CONST = 0.0016146
    DPke = _CONST * fluid_density / (Cd * 2) * (fluid_velocity_outlet ** 2 - fluid_velocity_inlet ** 2)
    return DPke

def calc_DPf(friction_factor, fluid_density, fluid_velocity, hydraulic_diameter, length):
    #fluid_density:ppg    fluid_velocity:ft/s   hydraulic_diameter:in     length:ft    output:psi
    #pressure drop due to friction
    _CONST = 0.01938
    DPf = _CONST * 2 * friction_factor * fluid_density * length * fluid_velocity ** 2 / hydraulic_diameter
    return DPf

def calc_DPf_full(model: str, fluid_density, fluid_viscosity, fluid_velocity, hydraulic_diameter, length, eccentricity, roughness, nPrime = 1, kPrime = 0.0000208854342245726):
    #fluid_density:ppg    fluid_viscosity:cP   fluid_velocity:ft/s   hydraulic_diameter,roughness:in     length:ft    output:psi
    #pressure drop due to friction, selects model type
    if model == "Newtonian":
        kPrime = calc_kPrime(nPrime, fluid_viscosity, 170)
        NRe = calc_NRe_newton(fluid_velocity, hydraulic_diameter, fluid_density, fluid_viscosity)
        fluid_friction = calc_friction_colebrook(hydraulic_diameter, NRe, roughness)
        #fluidFriction = calcFrictionChen(hydraulic_diameter, NRe, roughness)
    elif model == "Power Law":
        NRe = calc_NRe_power(fluid_velocity, hydraulic_diameter, fluid_density, nPrime, kPrime)
        fluid_friction = calc_friction_power_explicit(NRe, nPrime)
        #fluid_friction = calc_friction_power_implicit(NRe, nPrime)
    
    _CONST = 0.01938         #needs to be here since we calculated NRe, friction and eccentricity
    DPf_full = _CONST * 2 * fluid_friction * fluid_density * length * fluid_velocity ** 2 / hydraulic_diameter * calc_eccentricity_factor_powerlaw(nPrime, kPrime, NRe, eccentricity)
    return DPf_full







#*************************************************************************
#SETTLING AND TRANSPORT EQUATIONS
#*************************************************************************
def calc_settling_velocity(solid_density, fluid_density, hydraulic_diameter, NRe_Newton):    #general solution based on NRe
    _CONST = (1 / 12) ** 0.5
    NRe = NRe_Newton
    if NRe <= 2:
        Cd = 24 / NRe
    elif NRe > 2 and NRe <= 1000:
        Cd = 18 * NRe ** -0.6
    elif NRe > 1000:
        Cd = 0.44
    
    calc_settling_velocity = (4 / 3 * GRAVITY_CONSTANT * hydraulic_diameter / Cd * (solid_density - fluid_density) / fluid_density) ** 0.5
    return calc_settling_velocity

def calc_settling_Stokes(solid_density, fluid_density, hydraulic_diameter, fluid_viscosity):   #solution given by Stoke#s equation, Cd = 24/Re
    #solid_density,fluid_density:ppg      hydraulic_diameter:in     fluid_viscosity:cP    output:ft/s
    _CONST = 77.307
    calc_settling_Stokes = _CONST * ((solid_density - fluid_density) * GRAVITY_CONSTANT * hydraulic_diameter ** 2) / (18 * fluid_viscosity)
    return calc_settling_Stokes

def calc_settling_intermediate(solid_density, fluid_density, hydraulic_diameter, fluid_viscosity): #solution given by Cd = 18*Re**-0.6
    #solid_density,fluid_density:ppg      hydraulic_diameter:in     fluid_viscosity:cP    output:ft/s
    _CONST = 3.169
    calc_settling_intermediate = _CONST * (2 * GRAVITY_CONSTANT / 27 * (solid_density - fluid_density) / fluid_density) ** (5 / 7) * (hydraulic_diameter ** (8 / 7)) / ((fluid_viscosity / fluid_density) ** (3 / 7))
    return calc_settling_intermediate

def calc_settling_Newtonian(solid_density, fluid_density, hydraulic_diameter):       #solution given by Cd = 0.44
    #solid_density,fluid_density:ppg      hydraulic_diameter:in     output:ft/s
    _CONST = 0.289
    f = 0.44            #for Re > 1000
    calc_settling_Newtonian = _CONST * (3 * GRAVITY_CONSTANT * hydraulic_diameter * (solid_density - fluid_density) / fluid_density) ** 0.5
    return calc_settling_Newtonian

def calc_settling_power_law(solid_density, fluid_density, hydraulic_diameter, nPrime, kPrime):    #power law fluid settling
    #solid_density,fluid_density:ppg      hydraulic_diameter:in   output:ft/s
    #calc_settling_power_law = (0.8667 * (solid_density - fluid_density) / kPrime) * (pipeDiameter ** (1 / nPrime) * (2 * nPrime + 1) * solidDiameter) / (108 * nPrime)       #Halliburton CT manual
    #calc_settling_power_law = (((solid_density - fluid_density) * GRAVITY_CONSTANT * hydraulic_diameter ** (nPrime + 1)) / (2 * 3 ** (nPrime + 1) * kPrime)) ** (1 / nPrime)    #Economides Modern Hydraulic Fracturing
    calc_settling_power_law = (((solid_density - fluid_density) * GRAVITY_CONSTANT * hydraulic_diameter ** (nPrime + 1)) / (18 * 3 ** (nPrime - 1) * kPrime)) ** (1 / nPrime)    #Economides Reservoir Stimulation
    return calc_settling_power_law

def calc_settling_Budryck(solid_density, fluid_density, hydraulic_diameter): #intermediate solution with Budryck equation
    #solid_density,fluid_density:ppg      hydraulic_diameter:in     output:ft/s
    _CONST = 1 / 304.8
    calc_settling_Budryck = _CONST * 8.925 / (hydraulic_diameter * 25.4) * (math.sqrt(1 + 95 * (hydraulic_diameter * 25.4) ** 3 * (solid_density - fluid_density) / fluid_density) - 1)
    #calc_settling_Budryck = _CONST * 8.925 / (fydraulicDiameter) * (sqrt(1 + 95 * (hydraulic_diameter) ** 3 * (solid_density - fluid_density) / fluid_density) - 1)
    return calc_settling_Budryck

def calc_settling_Rittinger(solid_density, fluid_density, hydraulic_diameter): #turbulent solution with Rittinger equation
    #solid_density,fluid_density:ppg      hydraulic_diameter:in     output:ft/s
    _CONST = 1 / 304.8
    calc_settling_Rittinger = _CONST * 87 * ((solid_density - fluid_density) / fluid_density * hydraulic_diameter * 25.4) ** 0.5
    return calc_settling_Rittinger

def calc_settling_McCabeSmith(solid_density, fluid_density, hydraulic_diameter, fluid_viscosity): #intermediate solution with McCabe-Smith equation, 2<Re<500. SPE 187498
    #solid_density,fluid_density:ppg      hydraulic_diameter:in     output:ft/min
    _CONST = 1 / 60
    calc_settling_McCabeSmith = _CONST * (0.072 * GRAVITY_CONSTANT * unit(hydraulic_diameter, "m", "in") ** 1.6 * (solid_density - fluid_density) / (fluid_density ** 0.4 * fluid_viscosity ** 0.6)) ** 0.71
    return calc_settling_McCabeSmith

def calc_settling_turbulent(solid_density, fluid_density, hydraulic_diameter): 
    #turbulent solution with turbulent equation, 500<Re. SPE 187498
    #solid_density,fluid_density:ppg      hydraulic_diameter:in     output:ft/min
    _CONST = 1 / 60
    calc_settling_turbulent = _CONST * 1.74 * (GRAVITY_CONSTANT * (solid_density - fluid_density) / fluid_density * unit(hydraulic_diameter, "m", "in")) ** 0.5
    return calc_settling_turbulent

def calc_void_fraction(bulk_density, solid_density):
    #bulk_density,solid_density in same units
    calc_void_fraction = 1 - bulk_density / solid_density
    return calc_void_fraction

def calc_solids_fraction(solid_loading, solid_density):
    #also known as solids concentration, c
    calc_solids_fraction = solid_loading / (solid_loading + solid_density)
    return calc_solids_fraction

def calc_fluidization(hydraulic_diameter, solid_density, fluid_density, fluid_viscosity, porosity, sphericity):
    #hydraulic_diameter:in   solid_density,fluid_density:ppg   fluid_viscosity:cP    porosity,sphericity:fraction    output:ft/s
    #only good for small Re<10
    #according to McCabe, 30x fluidization velocity is ok to produce at
    _CONST = 77.307
    calc_fluidization = _CONST * ((solid_density - fluid_density) * GRAVITY_CONSTANT * hydraulic_diameter ** 2) / (150 * fluid_viscosity) * (porosity ** 3 * sphericity ** 2 / (1 - porosity))
    return calc_fluidization

def calc_angle_correction(inclination):
    #calcAngleCorrection = 1 + 2 * Inc / 45                             #Rubiandini model, for transport UP
    #calcAngleCorrection = 0.0342 * Inc - 0.000233 * Inc ** 2 - 0.213   #SPE 25872 Larsen, for transport UP
    calc_angle_correction = 0.022 * inclination - 0.9793                #Hang model
    return calc_angle_correction

def calc_horizontal_transport_Oroskar(casing_ID, solid_diameter, solid_density, fluid_density, fluid_viscosity, c, 
                                      x = 0.95, _CONSTY = 1.85,  _CONSTn1 = 0.1536,  _CONSTn2 = 0.3564,  _CONSTn3 = -0.378 ,  _CONSTn4 = 0.09,  _CONSTn5 = 0.3):
    #solid_density,fluid_density:ppg      casing_ID,solid_diameter:in     fluid_viscosity:cP    output:ft/s
    #c: solids concentration, loading/(loading+density)
    #for non-concentric flow
    _CONST1 = math.sqrt(1/12)   #constant, sqrt(1/12)
    _CONST2 = 927.6866          #constant for modified NRe to multiply by ft/s. See calc_NRe_newton function

    horizontal_transport_Oroskar = _CONST1 * math.sqrt(GRAVITY_CONSTANT * solid_diameter * (solid_density / fluid_density - 1)) * \
        _CONSTY * c ** _CONSTn1 * \
        (1 - c) ** _CONSTn2 * \
        (solid_diameter / casing_ID) ** _CONSTn3 * \
        (_CONST1 * _CONST2 * casing_ID * fluid_density / fluid_viscosity * math.sqrt(GRAVITY_CONSTANT * solid_diameter * (solid_density / fluid_density - 1))) ** _CONSTn4 * \
        x ** _CONSTn5
    return horizontal_transport_Oroskar

def calc_horizontal_transport_OroskarMod(hydraulic_diameter, solid_diameter, solid_density, fluid_density, fluid_viscosity, c, 
                                      x = 0.95, _CONSTY = 1.85,  _CONSTn1 = 0.1536,  _CONSTn2 = 0.3564,  _CONSTn3 = -0.378 ,  _CONSTn4 = 0.09,  _CONSTn5 = 0.3):
    #solid_density,fluid_density:ppg      hydraulic_diameter,solid_diameter:in     fluid_viscosity:cP    output:ft/s
    #c: solids concentration, loading/(loading+density)
    #modified Oroskar for gravel packing
    _CONST1 = math.sqrt(1/12)   #constant, sqrt(1/12)
    _CONST2 = 927.6866          #constant for modified NRe to multiply by ft/s. See calc_NRe_newton function
     
    #_CONSTY = 3.52, _CONSTn1 = -0.111, _CONSTn2 = -2.97, _CONSTn3 = -0.357, _CONSTn4 = -0.0595, _CONSTn5 = 0.3         #updated constants from regression fit, uncomment to compare
            
    horizontal_transport_OroskarMod = _CONST1 * math.sqrt(GRAVITY_CONSTANT * solid_diameter * (solid_density / fluid_density - 1)) * \
        _CONSTY * c ** _CONSTn1 * \
        (1 - c) ** _CONSTn2 * \
        (solid_diameter / hydraulic_diameter) ** _CONSTn3 * \
        (_CONST1 * _CONST2 * hydraulic_diameter * fluid_density / fluid_viscosity * math.sqrt(GRAVITY_CONSTANT * solid_diameter * (solid_density / fluid_density - 1))) ** _CONSTn4 * \
        x ** _CONSTn5
    return horizontal_transport_OroskarMod

def calc_horizontal_transport_SGS(equivalent_diameter, solid_diameter, solid_density, fluid_density, fluid_viscosity, c, dune_height, hole_diameter, 
                                    _CONSTY = 3.445, _CONSTn1 = -3.743, _CONSTn2 = 0.3075, _CONSTn3 = 0.065):
    #solid_density,fluid_density:ppg      equivalent_diameter,solid_diameter,hole_diameter:in     viscosity:cP  height:height of dune    output:ft/s
    #c: solids concentration, loading/(loading+density)
    _CONST1 = 1 #math.sqrt(1/12)   #constant, sqrt(1/12)
    _CONST2 = 927.6866          #constant for modified NRe to multiply by ft/s. See calc_NRe_newton function
              
    horizontal_transport_SGS = _CONST1 * math.sqrt(GRAVITY_CONSTANT * solid_diameter * (solid_density / fluid_density - 1)) * \
        _CONSTY * (1 - c) ** _CONSTn1 * \
        (solid_diameter / equivalent_diameter * (1 - dune_height / hole_diameter)) ** _CONSTn2 *  \
        (_CONST1 * _CONST2 * equivalent_diameter * fluid_density / fluid_viscosity * math.sqrt(GRAVITY_CONSTANT * solid_diameter * (solid_density / fluid_density - 1))) ** _CONSTn3
    return horizontal_transport_SGS

def calc_horizontal_transport_SGSl_alt(equivalent_diameter, solid_diameter, solid_density, fluid_density, fluid_viscosity, c, bed_width, wetted_perimeter, 
                                    _CONSTY = 1.15915872803019, _CONSTn1 = -3.58822426562324, _CONSTn2 = 0.377485123473551, _CONSTn3 = -0.272080143824705, _CONSTn4 = 0.05439676041480):
    #solid_density,fluid_density:ppg      hydraulic_diameter,solid_diameter,hole_diameter:in     viscosity:cP  height:height of dune    output:ft/s
    #c: solids concentration, loading/(loading+density)
    _CONST1 = math.sqrt(1/12)   #constant, sqrt(1/12)
    _CONST2 = 927.6866          #constant for modified NRe to multiply by ft/s. See calc_NRe_newton function
   
    horizontal_transport_SGS_alt = _CONST1 * math.sqrt(GRAVITY_CONSTANT * solid_diameter * (solid_density / fluid_density - 1)) * \
        _CONSTY * (1 - c) ** _CONSTn1 * \
        (equivalent_diameter / solid_diameter) **_CONSTn2  * \
        max(bed_width / wetted_perimeter, 0.05) ** _CONSTn3 * \
        (_CONST1 * _CONST2 * equivalent_diameter * fluid_density / fluid_viscosity * math.sqrt(GRAVITY_CONSTANT * solid_diameter * (solid_density / fluid_density - 1))) ** _CONSTn4
    return horizontal_transport_SGS_alt

def calc_horizontal_transport_Hang(hydraulic_diameter, solid_diameter, solid_density, fluid_density, slurryDensity, fluid_viscosity):
    #hydraulic_diameter:in   solid_diameter:in   solid_density:ppg    fluid_density:ppg   slurryDensity:ppg   fluid_viscosity:cP    output:ft/s
    _CONST1 = 9.912
    _CONST2 = 0.289

    V1 = _CONST1 * (0.0251 * GRAVITY_CONSTANT * solid_diameter * (solid_density / fluid_density - 1) * ((hydraulic_diameter * slurryDensity) / fluid_viscosity) ** 0.775) ** 0.816
    V2 = _CONST2 * (1.35 * 2 * GRAVITY_CONSTANT * hydraulic_diameter * (solid_density / fluid_density - 1)) ** 0.5
    if V1 > V2:
        horizontal_transport_Hang = V1
    else: 
        horizontal_transport_Hang = V2
    return horizontal_transport_Hang

def calc_horizontal_transport_SPE(hydraulic_diameter, solid_diameter, solid_density, fluid_density, fluid_viscosity): 
    #solid_density,fluid_density:ppg      hydraulic_diameter,solid_diameter:in     fluid_viscosity:cP      output:ft/s
    #from SPE97142
    _CONST = 3.281
    d1 = math.pi / 6 * solid_diameter * (solid_density - fluid_density) * GRAVITY_CONSTANT * (2 * fluid_viscosity / fluid_density) ** 0.5
    d2 = 162.4 * fluid_viscosity / hydraulic_diameter * (2 - solid_diameter / hydraulic_diameter) * (solid_diameter / hydraulic_diameter) * (hydraulic_diameter - solid_diameter) ** 0.5
    horizontal_transport_SPE = _CONST * (d1 / d2) ** (2 / 3)
    return horizontal_transport_SPE






#*************************************************************************
#GRAVEL PACK SLURRY EQUATIONS
#*************************************************************************
def calc_clean_volume(slurry_volume, solid_absVol, solid_loading, solid_loading_end=0):
    #slurry_volume,sandVolume:ft3    bulkDensity:lbm/ft3  solid_absVol:gal/lbm  solid_loading:ppg     output:gal
    #can also be used for clean volume add rate (gal/min) when using proppant feed rate (ft3/min)
    #clean_volume = sand_volume * bulkDensity / solid_loading
    if solid_loading_end == 0: 
        solid_loading_end = solid_loading
    clean_volume = slurry_volume / (1 + solid_absVol * (solid_loading + solid_loading_end) / 2)
    return clean_volume

def calc_slurry_volume(clean_volume, solid_absVol, solid_loading, solid_loading_end=0):
    #cleanVolume:volume    solid_absVol:gal/lbm   solid_loading:ppg     output:volume
    if solid_loading_end == 0: 
        solid_loading_end = solid_loading
    slurry_volume = clean_volume * (1 + solid_absVol * (solid_loading + solid_loading_end) / 2)
    return slurry_volume

def calc_slurry_density(clean_density, solid_absVol, solid_loading):
#clean_density:ppg   solid_absVol:gal/lbm   solid_loading:ppga   output:ppg
    slurry_density = (clean_density + solid_loading) / (1 + solid_absVol * solid_loading)
    return slurry_density

def calc_proppant_mass(slurry_volume, solid_absVol, solid_loading, solid_loading_end=0):
    #slurry_volume:gal   solid_absVol:gal/lbm   solid_loading:ppga   #output:lbm
    #can also be used for proppant mass add rate (lbm/min) when using slurry rate (gal/min)
    _CONST = 1
    if solid_loading_end == 0: 
        solid_loading_end = solid_loading
    proppant_mass = _CONST * slurry_volume / (1 / ((solid_loading + solid_loading_end) / 2) + solid_absVol)
    return proppant_mass

def calc_solids_concentration(solid_density, solid_loading):
    #solid_loading,solid_density:ppg(a)
    #commonly used for solids fraction in transport equations as c
    solids_concentration = solid_density / (solid_density + solid_loading)
    return solids_concentration

def create_Nolte_schedule(fluid_efficiency, solid_loading, solid_absVol, solid_mass):
    #solid_loading:ppga  absVol:gal/lbm  solid_mass:lbm

    pad_fraction = (1 - fluid_efficiency) ** 2
    #padFraction = (1 - fluid_efficiency) / (1 + fluid_efficiency)
    spatial_exponent = (1 - pad_fraction - fluid_efficiency) / fluid_efficiency

    proppant_concentration = []
    relative_concentration = []
    relative_time = []
    sand_mass = []
    slurry_volume = []
    output = []

    for i in range(solid_loading + 2):  #extra for 0 and 0.5ppga stages
        if i == 0:
            proppant_concentration.append(0)
            relative_time.append(0)
        elif i == 1:
            proppant_concentration.append(0.5)
            relative_time.append(math.exp(math.log(proppant_concentration[i] / solid_loading) / spatial_exponent))
        else:
            proppant_concentration.append(i - 1)
            relative_time.append(math.exp(math.log(proppant_concentration[i] / solid_loading) / spatial_exponent))

    #Newton method to find total proppant-slurry volume that contains specified sand volume
    #Mass Proppant in a stage = dT * TotalSlurry / (1/PPGA + absVol) * (PPGA/PPGA). To avoid error at 0 PPGA multiply by PPGA/PPGA
    exitConverged = 1
    loopCounter = 1
    total_slurry_volume = solid_mass / 100
    d_total_slurry_volume = 1
    while abs(exitConverged) > 0.1 and loopCounter < 1000:
        total_sand_mass1 = 0
        total_sand_mass2 = 0     
        sand_mass = []
        for i in range(solid_loading + 2):
            sand_mass.append((relative_time[i] - relative_time[i - 1]) * total_slurry_volume * proppant_concentration[i] / (1 + proppant_concentration[i] * solid_absVol))
            total_sand_mass1 = total_sand_mass1 + sand_mass[i]
        sand_mass = []
        for i in range(solid_loading + 2):
            sand_mass.append((relative_time[i] - relative_time[i - 1]) * (total_slurry_volume + d_total_slurry_volume) * proppant_concentration[i] / (1 + proppant_concentration[i] * solid_absVol))
            total_sand_mass2 = total_sand_mass2 + sand_mass[i]

        total_sand_mass1 = total_sand_mass1 - solid_mass
        total_sand_mass2 = total_sand_mass2 - solid_mass

        totalSlurryVolume1 = total_slurry_volume - total_sand_mass1 / ((total_sand_mass2 - total_sand_mass1) / d_total_slurry_volume)
        exitConverged = totalSlurryVolume1 - total_slurry_volume
        loopCounter = loopCounter + 1
        total_slurry_volume = totalSlurryVolume1

    #create slurry volume schedule
    output.append([0,total_slurry_volume * pad_fraction]) #/ (1 - pad_fraction)
    for i in range(1, solid_loading + 2):
        output.append([proppant_concentration[i],sand_mass[i] / proppant_concentration[i] * (1 + proppant_concentration[i] * absVol)])
    Nolte_schedule = output[:]
    return Nolte_schedule






#*************************************************************************
#SANDOUT EQUATIONS
#*************************************************************************
def calc_Darcy_sandout(casing_diameter, tubing_diameter, fluid_rate, fluid_viscosity, permeability, height):
#casing_diameter, tubing_diameter:in     flowRate:bbl/min   fluid_viscosity:cP    permeability:D  height:ft   output:psi
    area = unit(calc_area(casing_diameter, tubing_diameter), "in2", "ft2")
    calc_Darcy_sandout = height * fluid_viscosity * fluid_rate / (0.000783 * permeability * area)
    return calc_Darcy_sandout

def calc_Forcheimer_sandout(casing_diameter, tubing_diameter, fluid_rate, fluid_density, fluid_viscosity, permeability, height):
#casing_diameter, tubing_diameter:in     fluid_rate:bbl/min   fluid_density:ppg  fluid_viscosity:cP    permeability:D  height:ft  output:psi
    area = unit(calc_area(casing_diameter, tubing_diameter), "in2", "ft2")
#    calc_Forcheimer_sandout = height * (1279 * fluid_viscosity * fluid_rate * area ** 2 + 2.82 * fluid_density * fluid_rate ** 2 * permeability ** 0.698) / (permeability * area ** 2)
    calc_Forcheimer_sandout = height * (1279 * fluid_viscosity * fluid_rate * area + 2.82 * fluid_density * fluid_rate ** 2 * permeability ** 0.698) / (permeability * area ** 2)
#    calc_Forcheimer_sandout = height * (1279 * fluid_viscosity * fluid_rate * area + 4.63 * fluid_density * fluid_rate ** 2 * permeability ** 0.45) / (permeability * area ** 2)
    return calc_Forcheimer_sandout

def calc_Ergun_sandout(casing_diameter, tubing_diameter, solid_diameter, fluid_rate, fluid_viscosity, fluid_density, height, void_fraction, sphericity):
#casing_diameter, tubing_diameter:in     fluid_rate:bbl/min   fluid_viscosity:cP    fluid_density:ppg  height:ft   output:psi
    _CONST1 = 1.954 * 10 ** -6
    _CONST2 = 1.697 * 10 ** -4
    area = unit(calc_area(casing_diameter, tubing_diameter), "in2", "ft2")
    calc_Ergun_sandout = height * (_CONST1 * 150 * (fluid_rate / area * fluid_viscosity * (1 - void_fraction) ** 2) / ((sphericity * solid_diameter) ** 2 * void_fraction ** 3) + \
                                   (_CONST2 * 1.75 * (fluid_density * (fluid_rate / area) ** 2) / (sphericity * solid_diameter) * (1 - void_fraction) / void_fraction))
    return calc_Ergun_sandout

def calc_sand_fill_rate(casing_diameter, screen_diameter, fluid_rate, solid_loading, absVol, bulk_density, proppant_perf_weight):
#casing_diameter,screen_diameter:in   fluid_rate:bpm    solid_loading:ppga    absVol:gal/lbm   bulk_density:lbm/ft3  proppant_perf_weight:lbm/ft in perfs
    area = calc_area(casing_diameter, screen_diameter) / 144
    sand_rate = calc_proppant_mass(fluid_rate * 42, absVol, solid_loading)
    sand_volume = proppant_perf_weight + bulk_density * area
    calc_sand_fill_rate = sand_rate / sand_volume
    return calc_sand_fill_rate

def calc_settling_factor(solid_loading, bulk_density, absVol):   #BHI tech facts
#solid_loading:ppg    bulk_density:lbm/ft3  absVol:gal/lbm   output:%
#this gives the bulk settling factor. this is distinct from the volume concentration c, in that this uses (1 / bulkDensity)
#refer to calcSolidsFraction
    _CONST = 7.4608
    calc_settling_factor = _CONST * solid_loading * (1 / bulk_density) / (1 + absVol * solid_loading) * 100
    return calc_settling_factor




#*************************************************************************
#FRACTURE MECHANICS EQUATIONS
#*************************************************************************
def calc_perm_estimate(solid_diameter, porosity):  #Blake-Kozeny, SPE 31141
    #solid_diameter:mm
    calc_perm_estimate = solid_diameter ** 2 * porosity ** 3 / (150 * (1 - porosity) ** 2)
    return calc_perm_estimate

def calc_fracture_rate(permeability, formation_height, frac_gradient, formation_TVD, fluid_density, fluid_viscosity, wellbore_ID, Bo, skin):
    #permeability:mD    formation_height,formation_TVD:ft  frac_gradient:psi/ft    fluid_density:ppg    fluid_viscosity:cP   wellbore_ID:in   Bo:stb/rb
    _CONST = 1 / 24 / 60
    calc_fracture_rate = _CONST * calc_Darcy_IPR(permeability, formation_height, frac_gradient * formation_TVD, 0, fluid_density * formation_TVD * 0.052, Bo, fluid_viscosity, calc_shape_radial(500 * 12, wellbore_ID / 2, skin, "Pseudoradial"))
    return calc_fracture_rate

def calc_injectivity_index(rate, pressure):
#rate:bpm   pressure:psi    output:bpd/psi
#pressure should be total pressure over BHP, e.g. surface annulus + overbalance
    _CONST = 1440
    calc_injectivity_index = rate * _CONST / pressure
    return calc_injectivity_index

def calcEatonClosure(overburdenPressure, porePressure, poissonsRatio, biotNumber):
#overburdenPressure:psi     porePressure:psi    output:psi      #psi/ft can be used also but must be consistent!
    calcEatonClosure = poissonsRatio / (1 - poissonsRatio) * (overburdenPressure - biotNumber * porePressure) + biotNumber * porePressure
    return calcEatonClosure
    
def calcEatonFracImpermeable(overburdenPressure, porePressure, poissonsRatio, biotNumber):
#overburdenPressure:psi     porePressure:psi    output:psi      #psi/ft can be used also but must be consistent!
    calcEatonFracImpermeable = 2 * poissonsRatio / (1 - poissonsRatio) * (overburdenPressure - biotNumber * porePressure) + biotNumber * porePressure
    return calcEatonFracImpermeable

def calcEatonFracPermeable(overburdenPressure, porePressure, poissonsRatio, biotNumber):
#overburdenPressure:psi     porePressure:psi    output:psi      #psi/ft can be used also but must be consistent!
    calcEatonFracPermeable = 2 * poissonsRatio * (overburdenPressure - biotNumber * porePressure) + biotNumber * porePressure
    return calcEatonFracPermeable

def calcSGSClosure(overburdenPressure, porePressure):
#overburdenPressure:psi     porePressure:psi
    calcSGSClosure = 0.488 * (overburdenPressure - porePressure) + porePressure
    return calcSGSClosure

def calcSGSFrac(overburdenPressure, porePressure):
#overburdenPressure:psi     porePressure:psi
    calcSGSFrac = 0.6934 * (overburdenPressure - porePressure) + porePressure
    return calcSGSFrac

def calcOverburden(seawaterDensity, rockDensity, liquidDensity, formationDepth, seabedDepth, phi_0):    #from Applied Drilling Engineering
#seawaterDensity, rockDensity, liquidDensity:ppg     formationDepth, seabedDepth:ft     output:psi
#typical constants: phi_0~0.40  seawaterDensity~1.029SG     rockDensity~2.6SG   liquidDensity~1.074SG
#equation will provide overburden pressure in offshore environments. the resultant gradient by dividing by TVD will give gradient to surface
    bedThickness = formationDepth - seabedDepth

    if bedThickness // 1000 == 0:     #calculate to nearest rounded thousand
        bulkDensity = 1.95
        phi = 0.43
    elif bedThickness // 1000 == 1:
        bulkDensity = 2.02
        phi = 0.38
    elif bedThickness // 1000 == 2:
        bulkDensity = 2.06
        phi = 0.35
    elif bedThickness // 1000 == 3:
        bulkDensity = 2.11
        phi = 0.32
    elif bedThickness // 1000 == 4:
        bulkDensity = 2.16
        phi = 0.29
    elif bedThickness // 1000 == 5:
        bulkDensity = 2.19
        phi = 0.27
    elif bedThickness // 1000 == 6:
        bulkDensity = 2.24
        phi = 0.24
    elif bedThickness // 1000 == 7:
        bulkDensity = 2.27
        phi = 0.22
    elif bedThickness // 1000 == 8:
        bulkDensity = 2.29
        phi = 0.2
    elif bedThickness // 1000 == 9:
        bulkDensity = 2.33
        phi = 0.18
    elif bedThickness // 1000 == 10:
        bulkDensity = 2.35
        phi = 0.16
    elif bedThickness // 1000 == 11:
        bulkDensity = 2.37
        phi = 0.15
    elif bedThickness // 1000 == 12:
        bulkDensity = 2.38
        phi = 0.14
    elif bedThickness // 1000 == 13:
        bulkDensity = 2.4
        phi = 0.13
    elif bedThickness // 1000 == 14:
        bulkDensity = 2.41
        phi = 0.12
    elif bedThickness // 1000 == 15:
        bulkDensity = 2.43
        phi = 0.11
    elif bedThickness // 1000 == 16:
        bulkDensity = 2.44
        phi = 0.1
    elif bedThickness // 1000 == 17:
        bulkDensity = 2.45
        phi = 0.098
    elif bedThickness // 1000 == 18:
        bulkDensity = 2.46
        phi = 0.092
    elif bedThickness // 1000 == 19:
        bulkDensity = 2.47
        phi = 0.085
    elif bedThickness // 1000 == 20:
        bulkDensity = 2.48
        phi = 0.079
    else:
        bulkDensity = -1 * 10 ** 9 * (bedThickness) ** 2 + 5 * 10 ** -5 * (bedThickness) + 1.9676  #from polynomial curve fit
        phi = 9 * 10 ** -10 * (bedThickness) ** 2 - 3 * 10 ** -5 * (bedThickness) + 0.4164         #from polynomial curve fit
    
    if phi == 0 or phi_0 == 0:
        k = 0
    else:
        k = math.log(phi_0 / phi) / bedThickness
    
    calcOverburden = 0.052 * (seawaterDensity * seabedDepth + rockDensity * bedThickness) - (rockDensity - liquidDensity) * 0.052 * phi_0 / k * (1 - math.exp(-k * bedThickness))
    return calcOverburden

def calc_formation_mean_stress(overburden_pressure, closure_pressure, pore_pressure):      #from GWong
    #overburden_pressure,closure_pressure,pore_pressure:psi/ft    output:psi
    calc_formation_mean_stress = (overburden_pressure + 2 * closure_pressure) / 3 - pore_pressure
    return calc_formation_mean_stress

def calc_formation_Youngs_modulus(mean_stress):      #from GWong
    #mean_stress:psi     output:psi
    calc_formation_Youngs_modulus = 85.053 * mean_stress + 515477
    return calc_formation_Youngs_modulus

def calc_formation_toughness(Youngs_modulus, poissons_ratio, gamma):
    #Youngs_modulus:psi      gamma:in*lbm/in2, estimates are 0.5-2.0      output:psi*sqrt(in)
    calc_formation_toughness = math.sqrt(gamma * Youngs_modulus / (1 - poissons_ratio ** 2))
    return calc_formation_toughness

def calcElasticShearModulus(Youngs_modulus, poissons_ratio):      #from Economides, Petroleum Production Systems
    calcElasticShearModulus = Youngs_modulus / (2 * (1 + poissons_ratio))
    return calcElasticShearModulus

def calcPlaneStrainModulus(Youngs_modulus, poissons_ratio):      #from Economides, Petroleum Production Systems
    calcPlaneStrainModulus = Youngs_modulus / (1 - poissons_ratio ** 2)
    return calcPlaneStrainModulus

def calc_fracture_width_PKN(fluid_rate, fluid_viscosity, poissons_ratio, halfLength, Youngs_modulus, widthType):     #from Economides, Petroleum Production Systems
#fluidRate:bpm   fluid_viscosity:cP    halfLength:ft   Youngs_modulus:psi   WidthType:max or average    #output:in
    _CONST = 0.13
    if widthType == "max": gamma = 1
    else: gamma = math.pi / 4 * 0.75
    
    calc_fracture_width_PKN = _CONST * 2.31 * (fluid_rate * fluid_viscosity * 2 * (1 - poissons_ratio ** 2) * halfLength / Youngs_modulus) ** 0.25 * gamma
    return calc_fracture_width_PKN

def calcPKNWidthPower(fluidRate, nPrime, kPrime, halfLength, fracHeight, Youngs_modulus, widthType):     #from Economides, Petroleum Production Systems
#fluidRate:bpm   halfLength, fracHeight:ft   Youngs_modulus:psi   WidthType:max or average    #output:in
    _CONST = 1
    if widthType == "max": gamma = 1
    else: gamma = math.pi / 4 * 0.75
    
    calcPKNWidthPower = _CONST * 12 * ((128 / 3 / math.pi) * (nPrime + 1) * ((2 * nPrime + 1) / nPrime) ** nPrime * (0.9775 / 144) * (5.61 / 60) ** nPrime) ** (1 / (2 * nPrime + 2)) * \
    (fluidRate ** nPrime * kPrime * halfLength * fracHeight ** (1 - nPrime) / Youngs_modulus) ** (1 / (2 * nPrime + 2)) * gamma
    return calcPKNWidthPower
    
def calcKGDWidth(fluidRate, fluid_viscosity, poissons_ratio, halfLength, Youngs_modulus, fracHeight, widthType):     #from Economides, Petroleum Production Systems
#fluidRate:bpm   fluid_viscosity:cP    halfLength, fracHeight:ft   Youngs_modulus:psi     widthType:max or average    #output:in
    _CONST = 0.128
    if widthType == "max": gamma = 1
    else: gamma = math.pi / 4
    
    calcKGDWidth = _CONST * 2.27 * (fluidRate * fluid_viscosity * 2 * (1 - poissons_ratio ** 2) * halfLength ** 2 / (Youngs_modulus * fracHeight)) ** 0.25 * gamma
    return calcKGDWidth

def calcPennyWidth(fluidRate, fluid_viscosity, poissons_ratio, fracRadius, Youngs_modulus, widthType):       #from Guo & Ghalambor
#fluidRate:bpm   fluid_viscosity:cP    fracRadius:ft   Youngs_modulus:psi   WidthType:max or average    #output:in
    if widthType == "max": _CONST = 0.85
    else: _CONST = 2.56 / 12
    
    calcPennyWidth = _CONST * (fluid_viscosity * fluidRate * (1 - poissons_ratio) * fracRadius / Youngs_modulus) ** 0.25
    return calcPennyWidth
    
def calcNetPressureWidth(Youngs_modulus, poissons_ratio, characteristicLength, netPressure, model, widthType): #from Economides, Reservoir Stimulation
#Youngs_modulus,netPressure:psi  characteristicLength:ft     output:in
    _CONST = 1 / 12
    
    if widthType == "max": gamma = 1
    else: gamma = math.pi / 4 * 0.75
    
    planeStrainModulus = calcPlaneStrainModulus(Youngs_modulus, poissons_ratio)
    
    if model == "PKN": fractureStiffness = _CONST * 2 * planeStrainModulus / math.pi / characteristicLength      #PKN, characteristicLength is height
    elif model == "KGD": fractureStiffness = _CONST * planeStrainModulus / math.pi / characteristicLength         #KGD, characteristicLength is length
    elif model == "Penny": fractureStiffness = _CONST * 3 * math.pi * planeStrainModulus / 16 / characteristicLength #Penny, characteristicLength is radius
    
    calcNetPressureWidth = netPressure / fractureStiffness * gamma
    return calcNetPressureWidth

def calcFCD(fractureWidth, halfLength, formationPermeability, fracturePermeability, formationHeight, fractureHeight):
#fractureWidth:in   halfLength,fractureHeight,formationHeight:ft   formationPermability,fracturePermeability:mD
#if using kf-w,: use width of 12 in
    _CONST = 1 / 12
    calcFCD = _CONST * (fractureWidth * fracturePermeability) / (halfLength * formationPermeability) * (fractureHeight / formationHeight)
    return calcFCD
    
def calcFractureSkin(FCD, halfLength, wellboreRadius):
#halfLength,wellboreRadius:ft
    u = math.log(FCD)
    calcFractureSkin = (1.65 - 0.328 * u + 0.116 * u ** 2) / (1 + 0.18 * u + 0.064 * u ** 2 + 0.005 * u ** 3) - math.log(halfLength / wellboreRadius)
#   calcFractureSkin = log(wellboreRadius / halfLength * (2 + (pi / 2) ** 2 / FCD))
    return calcFractureSkin

def calcEquivalentWellboreRadius(wellboreRadius, skin):
    calcEquivalentWellboreRadius = wellboreRadius * math.exp(-skin)
    return calcEquivalentWellboreRadius

def calcFOI(reservoirRadius, wellboreRadius, equivalentwellboreRadius):
    calcFOI = math.log(reservoirRadius / wellboreRadius) / math.log(reservoirRadius / equivalentwellboreRadius)
    return calcFOI

def calcFractureVolume(fractureWidth, halfLength, fractureHeight):
#fractureWidth:in   halfLength,fractureHeight:ft    output:ft3
    _CONST = 1 / 12
    calcFractureVolume = 2 * _CONST * fractureWidth * halfLength * fractureHeight
    return calcFractureVolume

def calcFractureMass(fractureWidth, halfLength, fractureHeight, bulkDensity):
#fractureWidth:in   halfLength,fractureHeight:ft    bulkDensity:lbm/ft3      output:lbm proppant
    _CONST = 1 / 12
    calcFractureMass = 2 * _CONST * fractureWidth * halfLength * fractureHeight * bulkDensity
    return calcFractureMass

def calcUFD(FCD, halfLength, reservoirLength, proppantVolume, fracturePermeability, formationPermeability, formationHeight):
#halfLength,reservoirLength,formationHeight:ft      proppantVolume:ft3      fracturePermeability,formationPermeability:mD       output:[ft,in]
    _CONST = 12
    Ix = 2 * halfLength / reservoirLength
    Nprop = FCD * Ix ** 2
    xfOpt = ((proppantVolume / 2 * fracturePermeability) / (FCD * formationHeight * formationPermeability)) ** 0.5
    wOpt = _CONST * ((FCD * proppantVolume / 2 * formationPermeability) / (formationHeight * fracturePermeability)) ** 0.5
    
    output = [Ix, Nprop, xfOpt, wOpt]
    calcUFD = output[:]
    return calcUFD

def calcFractureWidth(arealProppantConcentration, solid_density, porosity):
#fractureWidth:in   solid_density:ppg    porosity:fraction
    calcFractureWidth = arealProppantConcentration * 1.604 / (solid_density * (1 - porosity))
    return calcFractureWidth

def calcArealProppantConcentration(fractureWidth, solid_density, porosity):
#fractureWidth:in   solid_density:ppg    porosity:fraction
    calcArealProppantConcentration = fractureWidth * solid_density * (1 - porosity) / 1.604
    return calcArealProppantConcentration

def calcEquivalentTime(dt, timePumped):
    calcEquivalentTime = dt * timePumped / (dt + timePumped)
    return calcEquivalentTime

def calcGdefTime(timeTotal, timePumped, alpha):
    dt = (timeTotal - timePumped) / (timePumped)
    if alpha == 0.5:
        g0 = math.pi / 2
        gdTd = (1 + dt) * math.asin((1 + dt) ** -0.5) + dt ** 0.5
    elif alpha == 1:
        g0 = 4 / 3
        gdTd = 4 / 3 * ((1 + dt) ** 1.5 - dt ** 1.5)
    else:
        g0 = 4 / 3 * ((1 + (0)) ** 1.5 - (0) ** 1.5 - 1)
        gdTd = 4 / 3 * ((1 + dt) ** 1.5 - dt ** 1.5 - 1)
    
    calcGdefTime = 4 / math.pi * (gdTd - g0)
    return calcGdefTime
    




#*************************************************************************
#ALPHA BETA PACKING
#*************************************************************************
#Alpha-Beta calculations, Hang solution
def calc_Cfunction(openhole_ID, dune_height):
    #openhole_ID,dune_height:ft      output:ft2
    Cfunction = 2 * math.sqrt(dune_height * (openhole_ID - dune_height))
    return Cfunction

def calc_Sfunction(openhole_ID, dune_height):
    #openhole_ID,dune_height:ft      output:ft2
    Sfunction = openhole_ID * math.acos((2 * dune_height - openhole_ID) / openhole_ID)
    return Sfunction

def calc_Gfunction(openhole_ID, dune_height):
    #outer_diameter,dune_height:ft      output:ft2
    Gfunction = 0.25 * (openhole_ID * calc_Sfunction(openhole_ID, dune_height) - (2 * dune_height - openhole_ID) * calc_Cfunction(openhole_ID, dune_height))
    return Gfunction

def calc_Pfunction(openhole_ID, dune_height, H1, H2, screen_OD):
    #outer_diameter,dune_height,H1,H2,inner_diameter:in      output:in
    if dune_height < H1: 
        Pfunction = calc_Cfunction(openhole_ID, dune_height) + calc_Sfunction(openhole_ID, dune_height) + math.pi * screen_OD
    elif dune_height > H2: 
        Pfunction = calc_Cfunction(openhole_ID, dune_height) + calc_Sfunction(openhole_ID, dune_height)
    else: 
        Pfunction = calc_Cfunction(openhole_ID, dune_height) - calc_Cfunction(screen_OD, dune_height - H1) + calc_Sfunction(openhole_ID, dune_height) + calc_Sfunction(screen_OD, dune_height - H1)
    return Pfunction

def calc_Afunction(openhole_ID, dune_height, H1, H2, screen_OD):
    #outer_diameter,dune_height,H1,H2,inner_diameter:in    output:in2
    if dune_height < H1:
        Afunction = calc_Gfunction(openhole_ID, dune_height) - math.pi / 4 * screen_OD ** 2
    elif dune_height > H2:
        Afunction = calc_Gfunction(openhole_ID, dune_height)
    else:
        Afunction = calc_Gfunction(openhole_ID, dune_height) - calc_Gfunction(screen_OD, dune_height - H1)
    return Afunction 

def calc_alphawave_DH_Hang(openhole_ID, screen_OD, centralizer_OD, dune_height_ratio):
#openhole_ID, screen_OD, centralizer_OD:in    DHR:dimensionless       output:in
    H1 = (centralizer_OD - screen_OD) / 2
    H2 = H1 + screen_OD
    dune_height = dune_height_ratio * openhole_ID
    alphawave_DH_Hang = 4 * calc_Afunction(openhole_ID, dune_height, H1, H2, screen_OD) / calc_Pfunction(openhole_ID, dune_height, H1, H2, screen_OD)
    return alphawave_DH_Hang
    
#****************************************************************
def calc_tubing_dune_height(diameter, dune_height_ratio):     
    #height of dune in circular pipe. from Marks Standard Handbook, mechanics of fluids
    #diameter:in    dune_height_ratio:dimensionless     output:in
    theta = 2 * math.acos(2 * dune_height_ratio - 1)

    if dune_height_ratio >= 0.5:
        area = diameter ** 2 / 8 * (theta - math.sin(theta))
    else:
        area = diameter ** 2 / 8 * (theta + math.sin(theta))  

    perimeter = (math.sin(theta / 2) + theta / 2) * diameter
    hydraulic_diameter = 4 * area / perimeter
    equivalent_diameter = math.sqrt(4 * area / math.pi)

    output = []
    output = [theta * 180 / math.pi, area, perimeter, hydraulic_diameter, equivalent_diameter]
    return output

#Alpha-Beta Calculations
def calc_alphawave_dune_height(casing_ID, screen_OD, centralizer_OD, dune_height_ratio):
    #casing_ID, screen_OD, centralizer_OD,output:in
    h = dune_height_ratio * casing_ID
    g_i = min((casing_ID - screen_OD), (centralizer_OD - screen_OD)) / 2
    
    h_o = min(max(h, 0), casing_ID)
    h_i = min(max(h - g_i, 0), screen_OD)
    
    theta_o = 2 * math.asin((casing_ID - 2 * h_o) / casing_ID) + math.pi
    theta_i = 2 * math.asin((screen_OD - 2 * h_i) / screen_OD) + math.pi

    area_o = casing_ID ** 2 / 8 * (theta_o + math.sin(theta_o - math.pi))
    area_i = screen_OD ** 2 / 8 * (theta_i + math.sin(theta_i - math.pi))

    perimeter_o = casing_ID * theta_o / 2 + casing_ID * math.sin((theta_o - math.pi) / 2)
    perimeter_i = screen_OD * theta_i / 2 - screen_OD * math.sin((theta_i - math.pi) / 2)

    width_o = casing_ID * math.sin(theta_o - math.pi)/2
    width_i = screen_OD * math.sin(theta_i - math.pi)/2

    hydraulic_diameter = 4 * (area_o - area_i) / (perimeter_o + perimeter_i)
    equivalent_diameter = math.sqrt(4 * (area_o - area_i) / math.pi)

    return hydraulic_diameter, equivalent_diameter, area_o, area_i, perimeter_o, perimeter_i, width_o, width_i

#*************************************************************************
#SURVEY EQUATIONS
#*************************************************************************

def calc_dogleg_angle(a1, a2, e1, e2):
    _B = math.acos(math.cos(math.radians(a2 - a1)) - math.sin(math.radians(a1)) * math.sin(math.radians(a2)) * (1 - math.cos(math.radians(e2 - e1))))        
    return _B       #reminder that beta is used in literature, and references later functions

def calc_tvd(md1, md2, a1, a2, B, prior):   #angles in degree
    if B == 0: _f = 1    #satisfies the case of no inclination
    else:   _f = (2 / B) * math.tan(B / 2)
    tvd = (md2 - md1) / 2 * (math.cos(math.radians(a1)) + math.cos(math.radians(a2))) * _f + prior     #TVD
    return tvd

def calc_north_south(md1, md2, a1, a2, e1, e2, B, prior):    #angles in degree
    if B == 0: _f = 1    #satisfies the case of no inclination
    else:   _f = (2 / B) * math.tan(B / 2)
    northsouth = (md2 - md1) / 2 * (math.sin(math.radians(a1)) * math.cos(math.radians(e1)) + math.sin(math.radians(a2)) * math.cos(math.radians(e2))) * _f + prior   #northing
    return northsouth
        
def calc_east_west(md1, md2, a1, a2, e1, e2, B, prior):     #angles in degree
    if B == 0: _f = 1    #satisfies the case of no inclination
    else:   _f = (2 / B) * math.tan(B / 2)
    eastwest = (md2 - md1) / 2 * (math.sin(math.radians(a1)) * math.sin(math.radians(e1)) + math.sin(math.radians(a2)) * math.sin(math.radians(e2))) * _f + prior  #easting
    return eastwest

def calc_departure(northsouth, eastwest):
    departure = math.sqrt(northsouth ** 2 + eastwest ** 2)
    return departure

def calc_dls(md1, md2, a1, a2, e1, e2):      #angles in degree
    _x = math.sqrt(math.sin((math.radians(a2 - a1)) / 2) ** 2 + math.sin(math.radians(e2 - e1) / 2) ** 2 * math.sin(math.radians(a2 + a1) / 2) ** 2)
    _B = 2 * math.asin(_x)     #asin(x) = Atn(x / Sqr(-x ^ 2 + 1))
    dls = math.degrees(_B) * 100 / (md2 - md1)    #convert to deg/100ft
    """old code
    b = ((a2 - a1) - (sin(a1) * sin(a2) * (1 - cos(e2 - e1))))
    if b == 0: f = 1    #satisfies the case of no inclination
    else:   f = (2 / b) * tan(b / 2)"""
    return dls





#*************************************************************************
#TUBULAR STRESS ANALYSIS EQUATIONS
#*************************************************************************
#Burst and Collapse Equations
def calc_pipe_temperature_derating(temperature):
    #temperature:deg F      output:dimensionless
    pipe_temperature_derating = (1 - 0.0003 * (temperature - 75))
    return pipe_temperature_derating

def calc_API_tensile(Yp, tubing_OD, tubing_ID):
    API_tensile = Yp * math.pi * (tubing_OD**2 - tubing_ID**2) / 4
    return API_tensile

def calc_API_burst(Yp, tubing_OD, tubing_ID, eccentricity=0.875):
    #Yp,output:stress unit     tubing_OD,tubing_ID:length
    #Barlow's equation
    API_burst = Yp * eccentricity * (tubing_OD - tubing_ID) / tubing_OD
    return API_burst

def calc_Lame_stress(Ro, Ri, r, Po, Pi):
    #Ro,Ri,r:inch        Po,Pi,output:psi
    Lame_hoop_stress = -(Ri**2 * Ro**2 * (Po - Pi))/(Ro**2 - Ri**2)*(1/r**2) + (Pi*Ri**2 - Po*Ro**2)/(Ro**2 - Ri**2)
    Lame_radial_stress = (Ri**2 * Ro**2 * (Po - Pi))/(Ro**2 - Ri**2)*(1/r**2) + (Pi*Ri**2 - Po*Ro**2)/(Ro**2 - Ri**2)
    return Lame_hoop_stress, Lame_radial_stress

def calc_VM_stress_envelope(Yp, tubing_OD, tubing_ID):
    #Yp,output:psi      tubing_OD,tubing_ID:inch
    hoop_stress = []
    axial_stress = []
    Ro = tubing_OD / 2
    Ri = tubing_ID / 2
    _C = (1 - Ro**2 / Ri**2) / (Ro**2 / Ri**2 + 1)
    hoop_stress.append(2 * Yp / math.sqrt(3*(_C**2 - 2*_C + 1)))
    axial_stress.append(0.5*(_C + 1)*hoop_stress[-1] - 0.5*math.sqrt(4*Yp**2 - 3*(_C**2 - 2*_C + 1)*hoop_stress[-1]**2))
    
    for _x in range(0, 100):
        hoop_stress.append(hoop_stress[-1] - hoop_stress[0]/50)
        axial_stress.append(0.5*(_C + 1)*hoop_stress[-1] - 0.5*math.sqrt(4*Yp**2 - 3*(_C**2 - 2*_C + 1)*hoop_stress[-1]**2))
    for _x in range(0, 100):
        hoop_stress.append(hoop_stress[-1] + hoop_stress[0]/50)
        axial_stress.append(0.5*(_C + 1)*hoop_stress[-1] + 0.5*math.sqrt(4*Yp**2 - 3*(_C**2 - 2*_C + 1)*hoop_stress[-1]**2))
    return hoop_stress, axial_stress

def calc_VM_envelope(Yp, tubing_OD, tubing_ID, radius, eccentricity=0.875, temperature=75.0, DF_burst=1.0, DF_tension=1.0):
    #Yp:psi      tubing_OD,tubing_ID,radius:inch        temperature:deg F      output:psi & lbf/1000
    #use API eccentricity for default wall thickness
    hoop_stress = []
    axial_stress = []
    pressure = []
    tension = []
    Ro = tubing_OD / 2
    Ri = tubing_ID / 2
    _C = (1 - Ro**2 / radius**2) / (Ro**2 / radius**2 + 1)
    Yp = Yp * calc_pipe_temperature_derating(temperature)
    cross_section_area = math.pi * (Ro**2 - Ri**2)
    uniaxial_burst = calc_API_burst(Yp, tubing_OD, tubing_ID, eccentricity)
    hoop_stress.append(2 * Yp / math.sqrt(3*(_C**2 - 2*_C + 1)))
    axial_stress.append(0.5*(_C + 1)*hoop_stress[-1] - 0.5 * math.sqrt(4*Yp**2 - 3*(_C**2 - 2*_C + 1)*hoop_stress[-1]**2))
    
    for _x in range(0, 100):
        hoop_stress.append(hoop_stress[-1] - hoop_stress[0]/50)
        axial_stress.append(0.5*(_C + 1)*hoop_stress[-1] - 0.5 * math.sqrt(4*Yp**2 - 3*(_C**2 - 2*_C + 1)*hoop_stress[-1]**2))
    for _x in range(0,100):
        hoop_stress.append(hoop_stress[-1] + hoop_stress[0]/50)
        axial_stress.append(0.5*(_C + 1)*hoop_stress[-1] + 0.5 * math.sqrt(4*Yp**2 - 3*(_C**2 - 2*_C + 1)*hoop_stress[-1]**2))
    for _x in range(len(hoop_stress)):
        pressure.append(uniaxial_burst * hoop_stress[_x] / (Yp / math.sqrt(_C**2 - _C + 1)) / DF_burst)
        tension.append(axial_stress[_x] * cross_section_area / DF_tension / 1000)
    return pressure, tension

def calc_API_collapse(Yp, tubing_OD, tubing_ID, Youngs_modulus, Poissons_ratio, temperature=75.0, axial_stress=0.0, DF_collapse=1.0):
    #Yp,axial_stress,output:psi      outer_diameter,inner_diameter:inch     temperature:deg F
    Dt = tubing_OD / ((tubing_OD - tubing_ID) / 2) 
    Ypt = Yp * calc_pipe_temperature_derating(temperature)
    Ypa = Ypt * (math.sqrt(1-0.75 * (axial_stress/Ypt)**2) - 0.5*axial_stress/Yp)
    _A = 2.8762 + Ypa * 0.10679 * 10**-5 + Ypa**2 * 0.21302 * 10**-10 - Ypa**3 * 0.53132 * 10**-16
    _B = 0.026233 + Ypa * 0.50609 * 10**-6
    _C = -465.93 + Ypa * 0.030867 - Ypa**2 * 0.10483 * 10**-7 + Ypa**3 * 0.36989 * 10**-13
    _F = (46.95 * 10**6 * ((3*_B/_A)/(2+_B/_A))**3) / (Ypa * ((3*_B/_A)/(2+_B/_A)-_B/_A) * (1-(3*_B/_A)/(2+_B/_A))**2)
    _G = _F * _B / _A

    if Dt <= (math.sqrt((_A-2)**2 + 8*(_B+_C/Ypa))+(_A-2)) / (2*(_B+_C/Ypa)):
        API_collapse = 2 * Ypa * (Dt - 1)/(Dt**2)
    elif Dt <= (Ypa * (_A - _F)) / (_C + Ypa * (_B - _G)):
        API_collapse = Ypa * (_A/Dt - _B) - _C
    elif Dt <= (2+_B/_A) / (3*_B/_A):
        API_collapse = Ypa * (_F/Dt - _G)
    else:
        API_collapse = (2*Youngs_modulus)/(1-Poissons_ratio**2) * 1/(Dt * (Dt - 1)**2) / DF_collapse
    return API_collapse    

#Buckling and Tubemove Equations
def calc_dL_HookesLaw(length, force, youngsModulus, area):   #units of L
    dL_HookesLaw = -length * force / youngsModulus / area
    return dL_HookesLaw

def calc_force_HookesLaw(length, dL, youngsModulus, area):   #lbf
    force_HookesLaw = dL * youngsModulus * area / length
    return force_HookesLaw

def calc_buoyed_weight(pipeWt, tubingfluid_density, IDArea, casingfluid_density, ODArea): #lb/ft
    buoyed_weight = pipeWt + (tubingfluid_density * IDArea - casingfluid_density * ODArea) * 12 / 231
    return buoyed_weight

def calc_force_Pasley(buoyWt, angle, youngsModulus, momentInertia, rc):        #lbf
    criticalAngle = math.asin((1.94 / 2) ** 2 * rc * ((buoyWt / 12) / (youngsModulus * momentInertia)) ** (1 / 3))
    if (angle) >= criticalAngle:
        force_Pasley = math.sqrt(4 * (buoyWt / 12) * math.sin(angle) * youngsModulus * momentInertia / rc)
    else:
        force_Pasley = 1.94 * (youngsModulus * momentInertia * (buoyWt / 12) ** 2) ** (1 / 3)
    return force_Pasley
    
def calcNormalBucklingForce(rc, force, youngsModulus, momentInertia):        #lbf/ft
    calcNormalBucklingForce = rc * force ** 2 / (4 * youngsModulus * momentInertia) * 12
    return calcNormalBucklingForce

def bucklingDogleg(force, pasleyForce, youngsModulus, momentInertia, bucklingState):
    if bucklingState == "lateral":
        bucklingDogleg = 1.1227 / math.sqrt(2 * youngsModulus * momentInertia) * force ** 0.04 * (force - pasleyForce) ** 0.46
    elif bucklingState == "helical":
        bucklingDogleg = math.sqrt(force / (2 * youngsModulus * momentInertia))
    else:
        bucklingDogleg = 0
    return bucklingDogleg
    
def bucklingStrain(force, pasleyForce, youngsModulus, momentInertia, rc, bucklingState):  #unit strain
    if bucklingState == "lateral":
        bucklingStrain = -0.7285 * rc ** 2 / (4 * youngsModulus * momentInertia) * force ** 0.08 * (force - pasleyForce) ** 0.92
    elif bucklingState == "helical":
        bucklingStrain = -rc ** 2 / (4 * youngsModulus * momentInertia) * force
    else:
        bucklingStrain = 0
    return bucklingStrain

def calcDLBuckling(force, force1, axialBuoyWt, pasleyForce, youngsModulus, momentInertia, rc, bucklingState):  #from Mitchell, ft
    if bucklingState == "lateral":
        calcDLBuckling = -(rc ** 2 / (4 * youngsModulus * momentInertia * axialBuoyWt) * (force - pasleyForce) * (0.3771 * force - 0.3668 * pasleyForce))
    elif bucklingState == "helical":
        calcDLBuckling = -(rc ** 2 / (8 * youngsModulus * momentInertia * axialBuoyWt) * (force ** 2 - force1 ** 2))
    else:
        calcDLBuckling = 0
    return calcDLBuckling
    
def calcDLBalloon(length, poisson, youngsModulus, dPTbg, dDensityTbg, dPAnn, dDensityAnn, R, d):    #in
    #D is pressure drop of fluid through the pipe, r = OD/ID of tubing
    #calcDLBalloon = -1 * Length ** 2 * poisson / YoungsModulus * (dDensityTbg - r ** 2 * dDensityAnn - (1 + 2 * poisson) * D / (2 * poisson)) / (r ** 2 - 1) - 2 * Length * poisson / YoungsModulus * (dPTbg - r ** 2 * dPAnn) / (r ** 2 - 1)
    calcDLBalloon = -2 * length * poisson / youngsModulus * (dPTbg - R ** 2 * dPAnn) / (R ** 2 - 1)
    return calcDLBalloon

def calcDLTemperature(alpha, dTemp, length):   #ft
    calcDLTemperature = alpha * dTemp * length
    return calcDLTemperature

def calcDFPistonPacker(Ap, Ai, Ao, dPi, dPo):   #lbf, negative is pushing down against the packer (tensile force)
    calcDFPistonPacker = (Ap - Ai) * dPi - (Ap - Ao) * dPo
    return calcDFPistonPacker

def calcDFPistonUpset(Ai1, Ai2, Ao1, Ao2, dPi, dPo):    #lbf
    calcDFPistonUpset = (Ai1 - Ai2) * dPi - (Ao1 - Ao2) * dPo
    return calcDFPistonUpset

def calcDFBend(F2, ff, angle):  #angle in radians
#    calcDFBend = F2 * (exp(ff * angle) - 1) #calculates dF across a bend
    calcDFBend = 2 * F2 * math.sin(angle / 2)    #alternative for normal force calculation
#    calcDFBend = 0
    return calcDFBend

def calcSlackoffForce(force, axialBuoyWt, normalBuoyWt, normalBucklingWt, pulleyWt, friction, length):
    calcSlackoffForce = force - (axialBuoyWt - (normalBuoyWt + normalBucklingWt + pulleyWt) * friction) * length
    return calcSlackoffForce

def calcPickupForce(force, axialBuoyWt, normalBuoyWt, normalBucklingWt, pulleyWt, friction, length):
    calcPickupForce = force - (axialBuoyWt + (normalBuoyWt + normalBucklingWt + pulleyWt) * friction) * length
    return calcPickupForce

#def xFictForceSlackoff(F2, wBuoyedAxial, wBuoyedNormal, wBucklingNormal, wPulley, ffriction, L):
#    f1 = F2 - (wBuoyedAxial - (wBuoyedNormal + wBucklingNormal + wPulley) * ffriction) * L
#    FictForceSlackoff = f1
#
#def xFictForcePickup(F2, wBuoyedAxial, wBuoyedNormal, wBucklingNormal, wPulley, ffriction, L):
#    f1 = F2 - (wBuoyedAxial + (wBuoyedNormal + wBucklingNormal + wPulley) * ffriction) * L
#    FictForcePickup = f1
#





#*************************************************************************
#PRODUCTION EQUATIONS
#*************************************************************************
def calc_shape_radial(re, rw, S, state):
    #re and rw same units
    if state == "Steady State":
        calc_shape_radial = math.log(re / rw) + S          #steady state
    else:
        calc_shape_radial = math.log(0.472 * re / rw) + S  #pseudosteady state
    return calc_shape_radial

def calc_shape_horizontal(kh, kv, re, rw, h, L, S, state):
    #kh & kv same units, re & rw same units, h & L same units
    if state == "Steady State":
        St = 1          #steady state
    else:
        St = 0.472     #pseudosteady state
    Iani = math.sqrt(kh / kv)
    a = L / 2 * (0.5 + (0.25 + (re / (0.5 * L)) ** 4) ** 0.5) ** 0.5
    calc_shape_horizontal = math.log((a + math.sqrt(a ** 2 - (L / 2) ** 2)) / (L / 2)) + Iani * h / L * (math.log(St * Iani * h / (rw * (Iani + 1))) + S)
    return calc_shape_horizontal

def calc_shape_dietz(c, a, rw, S):
    #A and rw in same units
    gamma = 1.78    #constant
    calc_shape_dietz = 0.5 * math.log(4 * a / (gamma * c * rw ** 2)) + S
    return calc_shape_dietz

def calc_Darcy_IPR(k, h, pr, pb, pwf, B, m, shape):   #calculate oil well IPR using Darcy & Vogel correlation
    #k:mD    h:ft    pr,pb,pwf:psi   m:cP    B:bbl/stb   #output stb/day
    _CONST = 0.001127           #constant to convert units to stb/day

    if pwf < pb:  #Vogel solution
        qb = (2 * math.pi * _CONST * k * h * (pr - pb)) / (B * m * shape)
        qv = qb * (pb / 1.8) / (pr - pb)
        q = qb + qv * (1 - 0.2 * (pwf / pb) - 0.8 * (pwf / pb) ** 2)
    else:            #Darcy solution
        q = (2 * math.pi * _CONST * k * h * (pr - pwf)) / (B * m * shape)
    calc_Darcy_IPR = q
    return calc_Darcy_IPR

def calc_Wiggins_IPR(k, ko, kw, h, pr, pb, pwf, Bo, Bw, Mo, mw, shape, a):    #calculate 3 phase IPR using Wiggins# correlation
    #k:mD    h:ft    pr,pb,pwf:psi   m:cP    B:bbl/stb   #output stb/day
    _CONST = 0.001127           #constant to convert units stb/day

    if a == 0:       #oil solution
        if pwf < pb:  #saturated solution
            qb = (2 * math.pi * _CONST * ko * k * h * (pr - pb)) / (Bo * Mo * shape)
            qv = qb * (pb / 1.48) / (pr - pb)
            q = qb + qv * (1 - 0.52 * (pwf / pb) - 0.48 * (pwf / pb) ** 2)
        else:            #unsaturated solution
            q = (2 * math.pi * _CONST * ko * k * h * (pr - pwf)) / (Bo * Mo * shape)
        
    else:                #water solution
        if pwf < pb:  #saturated solution
            qb = (2 * math.pi * _CONST * kw * k * h * (pr - pb)) / (Bw * mw * shape)
            qv = qb * (pb / 1.28) / (pr - pb)
            q = qb + qv * (1 - 0.72 * (pwf / pb) - 0.28 * (pwf / pb) ** 2)
        else:            #unsaturated solution
            q = (2 * math.pi * _CONST * kw * k * h * (pr - pwf)) / (Bw * mw * shape)
        
    
    calc_Wiggins_IPR = q
    return calc_Wiggins_IPR

def calc_gas_IPR(k, h, pr, pwf, m, T, Z, shape):  #calculate gas well IPR
    #input D, ft, psi, cP, F    #output Mscf/day
    _CONST = 0.111807        #constant to convert units

    q = (2 * math.pi * _CONST * k * h * (pr ** 2 - pwf ** 2)) / (T * Z * m * shape)
    calc_gas_IPR = q
    return calc_gas_IPR

def calc_perforating_skin(wellRadius, PhaseAngle, perfLength, perfHeight, perfEHD, horizontalPerm, verticalPerm): 
    #calculate perforating skin
    if PhaseAngle == 0:
        a_theta = 0.25
        a1 = -2.091
        a2 = 0.0453
        b1 = 5.1313
        b2 = 1.8672
        c1 = 0.16
        c2 = 2.675
    elif PhaseAngle == 180:
        a_theta = 0.5
        a1 = -2.025
        a2 = 0.0943
        b1 = 3.0373
        b2 = 1.8115
        c1 = 0.026
        c2 = 4.532
    elif PhaseAngle == 120:
        a_theta = 0.648
        a1 = -2.018
        a2 = 0.0634
        b1 = 1.6136
        b2 = 1.777
        c1 = 0.0066
        c2 = 5.32
    elif PhaseAngle == 90:
        a_theta = 0.726
        a1 = -1.905
        a2 = 0.1038
        b1 = 1.5674
        b2 = 1.6935
        c1 = 0.0019
        c2 = 6.155
    elif PhaseAngle == 60:
        a_theta = 0.813
        a1 = -1.898
        a2 = 0.1023
        b1 = 1.3654
        b2 = 1.649
        c1 = 0.0003
        c2 = 7.509
    elif PhaseAngle == 45:
        a_theta = 0.86
        a1 = -1.788
        a2 = 0.2398
        b1 = 1.1915
        b2 = 1.6392
        c1 = 0.000046
        c2 = 8.791

    #SH, plane flow effect
    if PhaseAngle == 0:
        rwPrime = perfLength / 4
    else: 
        rwPrime = a_theta * (wellRadius + perfLength)
    
    SH = math.log(wellRadius / rwPrime)

    #SV, vertical converging effect
    hD = perfHeight / perfLength * (horizontalPerm / verticalPerm) ** 0.5
    rD = perfEHD / (4 * perfHeight) * (1 + (verticalPerm / horizontalPerm) ** 0.5)
    a = a1 * math.log10(rD) + a2
    B = b1 * rD + b2
    SV = 10 ** a * hD ** (B - 1) * rD ** B

    #SWB, wellbore effect
    rwD = wellRadius / (perfLength + wellRadius)
    SWB = c1 * math.exp(c2 * rwD)
    calc_perforating_skin = SH + SV + SWB
    return calc_perforating_skin

def calc_VcVs(rate, formationVolumeFactor, formationPerm, formationHeight, fluid_density, fluid_viscosity, skin, OD, ID, perfIntervalLength, perf_diameter, proppantPerm, betaFactor, minAnnGap):
    _CONST = 141.2
    DPSkinMechanical = rate * _CONST * formationVolumeFactor * fluid_viscosity * skin / (formationPerm * formationHeight)
   
    perfLength = (OD - ID) / 2
    perfArea = calc_area(perf_diameter, 0)
    B = 10646.7 * fluid_viscosity * perfLength / proppantPerm
    a = 1.57501 * 10 ** -9 * fluid_density * perfLength * betaFactor
    effectiveSPF = (formationVolumeFactor * rate) / ((-B + (B ** 2 + 4 * a * DPSkinMechanical) ** 0.5) / (2 * a) * perfArea) / perfIntervalLength
    
    casingVelocity = ((-B + (B ** 2 + 4 * a * DPSkinMechanical) ** 0.5) / (2 * a)) * 144 * 5.6146 / 24 / 60 / 60
    screenVelocity = casingVelocity / (1 + 5.1285 * minAnnGap + 9.0898 * minAnnGap ** 2 - 0.7236 * minAnnGap ** 3)
    
    output = [DPSkinMechanical, effectiveSPF, casingVelocity, screenVelocity]
    calc_VcVs = output[:]
    return calc_VcVs






#****************************************************************
#PVT CALCULATIONS
#****************************************************************
def calcBubblePoint(method, Rs, oilDensity, gasGravity, temperature, pressure=None, Tsp=None, Psp=None):
#Rs:scf/stb    OilDensity:API  GasGravity:SG  Temperature:F    Pressure:psi
#Recommended methods are Glaso, Standing, Lasater, Vasquez & Beggs, Petrosky, Al-Marhoun, Velarde

    if method == "Standing":
        pb = 18.2 * ((Rs / gasGravity) ** 0.83 * 10 ** (0.00091 * temperature - 0.0125 * oilDensity) - 1.4)
    
    elif method == "Elam":
        pb = Rs ** 0.702 / gasGravity ** 0.514 * math.exp(0.00348 * temperature - 0.0282 * oilDensity + 3.58)
    
    elif method == "Lasater":
        Mo = 6084 / (oilDensity - 5.9)
        xg = (1 + unit(oilDensity, "API", "liquid sg") / (0.000007521 * Rs * Mo)) ** -1
        pf = math.exp((xg - 0.15649) / 0.33705) - 0.59162
        pb = pf * (temperature + 459.67) / gasGravity
    
    elif method == "Vazquez & Beggs":
        if oilDensity <= 30:
            a = 27.64
            B = -11.172
            c = 0.9143
        else:
            a = 56.06
            B = -10.393
            c = 0.8425
        
        GasGravityc = gasGravity * (1 + 0.00005912 * oilDensity * Tsp * math.log10(Psp / 114.7))
        pb = (a * (Rs / GasGravityc) * 10 ** (B * oilDensity / (temperature + 459.67))) ** c
    
    elif method == "Glaso - non-volatile":
        x = (Rs / gasGravity) ** 0.816 * (temperature ** 0.172 / oilDensity ** 0.989)    #non-volatile oils
        pb = 10 ** (1.7669 + 1.7447 * math.log10(x) - 0.30218 * (math.log10(x)) ** 2)
    
    elif method == "Glaso - volatile":
        x = (Rs / gasGravity) ** 0.816 * (temperature ** 0.13 / oilDensity ** 0.989)    #volatile oil
        pb = 10 ** (1.7669 + 1.7447 * math.log10(x) - 0.30218 * (math.log10(x)) ** 2)
    
    elif method == "Labedi":
        pb = 6.0001 / gasGravity * (Rs ** 0.6714 * (temperature / oilDensity) ** 0.7097 * Tsp ** 0.08929) / (10 ** (Rs * 7.995 * 10 ** -5))
    
    elif method == "Owolabi":
        #pb = 55 + 0.8643 * ((Rs / GasGravity) ** 1.255 * (Temperature ** 0.172 / OilDensity ** 0.178))   #old solution
        pb = -987.56359 + 179.58816 * ((Rs / gasGravity) ** 0.48088266 * (temperature ** 0.09353815 / oilDensity ** 0.16648326))
    
    elif method == "Al-Marhoun 1985":
        oilDensity = unit(oilDensity, "API", "liquid sg")
        x = Rs ** 0.722569 * oilDensity ** 3.04659 / gasGravity ** 1.879109 * (temperature + 459.67) ** 1.302347
        pb = -64.13891 + 7.02362 * 10 ** -3 * x - 2.278475 * 10 ** -9 * x ** 2
    
    elif method == "Obomanu & Okpobiri":
        pb = ((Rs * temperature ** 0.497 * 10 ** 0.811) / (1.01136371 * gasGravity ** 2.15 * oilDensity ** 1.27)) ** 1.0787
    
    elif method == "Al-Marhoun 1988":
        oilDensity = unit(oilDensity, "API", "liquid sg")
        pb = 5.38088 * 10 ** -3 * Rs ** 0.715082 * oilDensity ** 3.1437 * (temperature + 459.67) ** 1.32657 / gasGravity ** 1.87784
    
    elif method == "Asgarpour-Viking":
        temperature = unit(temperature, "F", "C")
        a1 = 70.9815
        a2 = -0.0101
        a3 = -0.2514
        a4 = 0.4593
        a5 = 0.5093
        pb = a1 * gasGravity ** a2 * oilDensity ** a3 * temperature ** a4 * (Rs / 5.614583) ** a5
    
    elif method == "Asgarpour-Nisku":
        temperature = unit(temperature, "F", "C")
        a1 = 83.7883
        a2 = -0.4114
        a3 = -0.2697
        a4 = 0.281
        a5 = 0.6259
        pb = a1 * gasGravity ** a2 * oilDensity ** a3 * temperature ** a4 * (Rs / 5.614583) ** a5
    
    elif method == "Asgarpour-Leduc":
        temperature = unit(temperature, "F", "C")
        a1 = 193.77
        a2 = 0.0928
        a3 = -1.1369
        a4 = 0.7899
        a5 = 0.6519
        pb = a1 * gasGravity ** a2 * oilDensity ** a3 * temperature ** a4 * (Rs / 5.614583) ** a5
    
    elif method == "Al-Najjar, Al-Soof & Al-Khalisy":
        if oilDensity <= 30:
            a = 7.92
            B = 1.025
            c = -24.244
        else:
            a = 30.91
            B = 0.816
            c = -19.748
        
        pb = a * (Rs / gasGravity) ** B * math.exp(c * oilDensity / (temperature + 459.67))
    
    elif method == "Dokla & Osma":
        oilDensity = unit(oilDensity, "API", "liquid sg")
        pb = 8363.86 * Rs ** 0.724047 * oilDensity ** 0.107991 / (gasGravity ** 1.01049 * (temperature + 459.67) ** 0.952584)
    
    elif method == "Petrosky":
        pb = 112.727 * (Rs ** 0.577421 / (gasGravity ** 0.8439 * 10 ** (7.916 * 10 ** -4 * oilDensity ** 1.541 - 4.561 * 10 ** -5 * temperature ** 1.3911)) - 12.34)
    
    elif method == "Kartoatmodjo & Schmidt":
        if oilDensity <= 30:
            a = 0.05958
            B = 0.7972
            c = 13.1405
            d = 0.998602
        else:
            a = 0.0315
            B = 0.7587
            c = 11.289
            d = 0.914328
        
        GasGravityc = gasGravity * (1 + 0.1595 * oilDensity ** 0.4078 / Tsp ** 0.2466 * math.log10(Psp / 114.7))
        pb = (Rs / (a * GasGravityc ** B * 10 ** (c * oilDensity / (temperature + 459.67)))) ** d
    
    elif method == "Farshad":
        pb = 64.14 * (Rs ** 0.6343 / (gasGravity ** 1.15036 * 10 ** (7.97 * 10 ** -3 * oilDensity - 3.35 * 10 ** -4 * temperature)) - 7.2818)
    
    elif method == "Macary":
        pb = 204.257 * math.exp(7.7 * 10 ** -4 * temperature - 9.7 * 10 ** -3 * oilDensity - 0.4003 * gasGravity) * (Rs ** 0.51 - 4.7927)
    
    elif method == "Omar & Todd":
        pb = 0
    
    elif method == "Hasan":
        pb = 18.3 * ((Rs / gasGravity) ** 0.83 * 10 ** (0.00091 * temperature - 0.0125 * oilDensity) + 2.2)
    
    elif method == "De Ghetto":    #questionable if GasGravity or GasGravityc
        if oilDensity <= 10:
            pb = 10.7025 * (Rs / gasGravity) ** 0.8986 * 10 ** (0.00156 * temperature - 0.00169 * oilDensity)
        elif (oilDensity > 10) and (oilDensity <= 22.3):
            pb = (56.434 * Rs / (gasGravity * 10 ** (10.9267 * oilDensity / (temperature + 459.67)))) ** 0.8294
        elif (oilDensity > 22.3) and (oilDensity <= 31.1):
            pb = (Rs / (0.10084 * gasGravity ** 0.2556 * 10 ** (7.4576 * oilDensity / (temperature + 459.67)))) ** 1.0134
        else:
            pb = (Rs / (0.1347 * gasGravity ** 0.3873 * 10 ** (12.753 * oilDensity / (temperature + 459.67)))) ** 0.8536
        
    
    elif method == "Agip":    #questionable if GasGravity or GasGravityc
        pb = (37.966 * Rs / (gasGravity * 10 ** (9.441 * oilDensity / (temperature + 459.67)))) ** 0.8669
    
    elif method == "Almedhaideb":
        Bo = 1.122018 + (1.41 * 10 ** -6 * Rs * temperature / oilDensity ** 2)
        pb = -620.592 + 6.23087 * Rs * oilDensity / (gasGravity * Bo ** 1.38559) + 2.89868 * temperature
    
    elif method == "Elsharkawy":
        if oilDensity <= 30:
            pb = (Rs / (gasGravity ** 0.04439 * oilDensity ** 1.1394 * 10 ** (8.392 * 10 ** -4 * temperature - 2.188))) ** 1.0551194
        else:
            pb = (Rs / (gasGravity * 10 ** (0.4636 * oilDensity / temperature - 1.2179))) ** 0.847271
        
    
    elif method == "Khairy":
        pb = 49.3647 * Rs ** 0.5774 * temperature ** 0.6641 / (gasGravity ** 1.4676 * oilDensity ** 1.0305)
    
    elif method == "Al-Shammasi":
        oilDensity = unit(oilDensity, "API", "liquid sg")
        pb = oilDensity ** 5.527215 * (gasGravity * Rs * (temperature + 459.67)) ** 0.783716 / math.exp(1.841408 * oilDensity * gasGravity)
    
    elif method == "Levitan & Murtha":
        oilDensity = unit(oilDensity, "API", "liquid sg")
        pb = 14.7 * (Rs / gasGravity) ** 0.85 * oilDensity ** 5 * ((temperature + 459.67) / 519.67) ** 1.5
    
    elif method == "Velarde":
        #A = 9.73 * 10 ** -7 * GasGravity ** 1.1672608 * OilDensity ** 0.92987 * Temperature ** 0.247235 * pb ** 1.056052
        #B = 0.022339 * GasGravity ** -1.00475 * OilDensity ** 0.337711 * Temperature ** 0.132795 * pb ** 0.302065
        #C = 0.725167 * GasGravity ** -1.48548 * OilDensity ** -0.164741 * Temperature ** -0.09133 * pb ** 0.047094
        #pr = Pressure / pb
        #Rs = Rsb * (A * pr ** B + (1 - A) * pr ** C)
        #use above to solve for Rs in future
        x = 0.013098 * temperature ** 0.282372 - 8.2 * 10 ** -6 * oilDensity ** 2.176124
        pb = 1091.47 * (Rs ** 0.081465 * 10 ** x / gasGravity ** 0.161488 - 0.740152) ** 5.354891
    
    elif method == "Dindoruk & Christman":
        a = (1.42828 * 10 ** -10 * temperature ** 2.844591797 - 6.74896 * 10 ** -4 * oilDensity ** 1.225226436) / (0.033383304 + 2 * gasGravity ** 0.084226069 * Rs ** -0.272945957) ** 2
        pb = 1.86997927 * (Rs ** 1.221486524 * 10 ** a / gasGravity ** 1.370508349 + 0.011688308)
    
    else:
        pb = "Select Method"

    calcBubblePoint = pb
    return calcBubblePoint

def calcOilFormationFactor(method, Rs, oilDensity, gasGravity, temperature, pressure=None, pb=None, Psp=None, Tsp=None):
#Rs: scf/stb    OilDensity: API  GasGravity: SG  Temperature: F    Pressure: psi
#Recommended methods are Glaso, Standing, Lasater, Vasquez & Beggs, Petrosky, Al-Marhoun, Velarde

    if method == "Standing":
        oilDensity = unit(oilDensity, "API", "liquid sg")
        B = 0.972 + 1.47 * 10 ** -4 * (Rs * (gasGravity / oilDensity) ** 0.5 + 1.25 * temperature) ** 1.175
    
    elif method == "Elam":
        oilDensity = unit(oilDensity, "API", "liquid sg")
        B = math.exp(-0.0355 + 3.55 * 10 ** -4 * Rs(gasGravity / oilDensity) ** 0.5 + 7.1 * 10 ** -4 * temperature)
    
    elif method == "Vazquez & Beggs":
        if oilDensity <= 30:
            a1 = 4.677 * 10 ** -4
            a2 = 1.751 * 10 ** -5
            a3 = -1.8106 * 10 ** -8
        else:
            a1 = 4.67 * 10 ** -4
            a2 = 1.1 * 10 ** -5
            a3 = 1.337 * 10 ** -9
        
        B = 1 + a1 * Rs + (temperature - 60) * (oilDensity / gasGravity) * (a2 + a3 * Rs)
    
    elif method == "Glaso":
        oilDensity = unit(oilDensity, "API", "liquid sg")
        x = Rs * (gasGravity / oilDensity) ** 0.526 + 0.968 * temperature
        B = 10 ** (-6.58511 + 2.91329 * math.log10(x) - 0.27683 * (math.log10(x)) ** 2) + 1
    
    elif method == "Labedi":
        Bob = 0.9976 + 5.273 * 10 ** -4 * Rs + 2.6636 * 10 ** -8 * (temperature - 60) * (oilDensity * Psp) + 1.6982 * 10 ** -5 * oilDensity * (temperature - 60)
        if Bob <= 1.758:
            x = (3.61 * Rs ** 0.4625 * Bob ** 13.398 * Psp ** 0.0775) / 10 ** (3.231 * Bob)
            B = Bob - x * (1 - pressure / pb)
        else:
            x = 1.6339 - 9.152 * 10 ** -4 * Rs + 1.584 * 10 ** -7 * Rs ** 2
            B = Bob - (1 - pressure / pb) ** x
        
    elif method == "Owolabi":
        oilDensity = unit(oilDensity, "API", "liquid sg")
        #B = 0.9871 + 4.0689 * 10 ** -4 * (Rs * (GasGravity / OilDensity) ** 0.526 + 1.25 * Temperature) #old solution
        B = 0.9957 + 3.7921 * 10 ** -4 * (Rs * (gasGravity / oilDensity) ** 0.526 + 1.25 * temperature)
    
    elif method == "Al-Marhoun 1985": #recommended
        oilDensity = unit(oilDensity, "API", "liquid sg")
        x = Rs ** 0.501538 * gasGravity ** -0.145526 * oilDensity ** -5.220726
        B = 0.574095 + 7.723532 * 10 ** -4 * (temperature + 459.67) + 2.454005 * 10 ** -3 * x + 3.727676 * 10 ** -5 * x ** 2
    
    elif method == "Obomanu & Okpobiri":
        oilDensity = unit(oilDensity, "API", "liquid sg")
        if oilDensity < 0.876:
            B = 0.3321 + 1.404154 * 10 ** -4 * Rs + 4.1588128 * 10 ** -4 * Rs * gasGravity / oilDensity + 1.15861 * 10 ** -5 * (temperature + 459.67)
        else:
            B = 1.0232 + 2.725 * 10 ** -5 * (Rs * (gasGravity / oilDensity + temperature)) ** 0.79
        
    elif method == "Al-Marhoun 1988": #recommended
        oilDensity = unit(oilDensity, "API", "liquid sg")
        x = Rs ** 0.74239 * gasGravity ** 0.323294 * oilDensity ** -1.20204
        B = 0.497069 + 8.62963 * 10 ** -4 * (temperature + 459.67) + 1.82594 * 10 ** -3 * x + 3.18099 * 10 ** -6 * x ** 2
    
    elif method == "Asgarpour-Viking":
        temperature = unit(temperature, "F", "C")
        a1 = 0.1203
        a2 = 0.0645
        a3 = 0.2452
        a4 = 0.1118
        a5 = 0.2321
        a6 = 1.342
        a7 = -0.5811
        a8 = 0.183
        a9 = 0.1656
        x = a1 * gasGravity ** a2 * oilDensity ** a3 * temperature ** a4 * (Rs / 5.614583) ** a5
        B = a6 + a7 * x + a8 * x ** 2 + a9 * x ** 3
    
    elif method == "Asgarpour-Nisku":
        temperature = unit(temperature, "F", "C")
        a1 = 0.251
        a2 = 0.0724
        a3 = 0.00275
        a4 = 0.1538
        a5 = 0.2235
        a6 = -2.3211
        a7 = 7.039
        a8 = -5.006
        a9 = 1.33
        x = a1 * gasGravity ** a2 * oilDensity ** a3 * temperature ** a4 * (Rs / 5.614583) ** a5
        B = a6 + a7 * x + a8 * x ** 2 + a9 * x ** 3
    
    elif method == "Asgarpour-Leduc":
        temperature = unit(temperature, "F", "C")
        a1 = 0.1941
        a2 = 0.0136
        a3 = 0.0912
        a4 = 0.159
        a5 = 0.211
        a6 = 0.8603
        a7 = 0.7341
        a8 = -0.9378
        a9 = 0.4686
        x = a1 * gasGravity ** a2 * oilDensity ** a3 * temperature ** a4 * (Rs / 5.614583) ** a5
        B = a6 + a7 * x + a8 * x ** 2 + a9 * x ** 3
    
    elif method == "Al-Najjar":
        oilDensity = unit(oilDensity, "API", "liquid sg")
        B = 0.96325 + 4.9 * 10 ** -4 * (Rs * (gasGravity / oilDensity) ** 0.5 + 1.25 * temperature)
    
    elif method == "Ahmed":
        B = -0.12869353 + Rs ** 0.023484894 * (oilDensity ** 0.015966573 / gasGravity ** 0.021946351) - 4.5243973 * 10 ** -4 * temperature + \
            3.9063637 * 10 ** -6 * temperature ** 2 - 5.5542509 / temperature - 5.760322 * 10 ** -6 * pressure - \
            3.9528992 * 10 ** -9 * pressure ** 2 + 16.289473 / pressure + 3.8718887 * 10 ** -4 * Rs + \
            7.0703685 * 10 ** -8 * Rs ** 2 - 1.4358395 / Rs
    
    elif method == "Abdul-Majeed & Salman":
        oilDensity = unit(oilDensity, "API", "liquid sg")
        x = Rs ** 1.2 / (gasGravity ** 0.147 * oilDensity ** 5.222)
        B = 0.9657876 + 4.8141 * 10 ** -5 * x - 6.8987 * 10 ** -10 * x ** 2 + 7.73 * 10 ** -4 * temperature
    
    elif method == "Dokla & Osman":
        oilDensity = unit(oilDensity, "API", "liquid sg")
        x = Rs ** 0.773572 * (gasGravity ** 0.40402 / oilDensity ** 0.882605)
        B = 4.31935 * 10 ** -2 + 1.56667 * 10 ** -3 * (temperature + 459.67) + 1.39775 * 10 ** -3 * x + 3.80525 * 10 ** -6 * x ** 2
    
    elif method == "Petrosky":
        oilDensity = unit(oilDensity, "API", "liquid sg")
        B = 1.0113 + 7.2046 * 10 ** -5 * (Rs ** 0.3738 * gasGravity ** 0.2914 / oilDensity ** 0.6265 + 0.24626 * temperature ** 0.5371) ** 3.0936
    
    elif method == "Kartoatmodjo & Schmidt":   #recommended
        gasGravityc = gasGravity * (1 + 0.1595 * oilDensity ** 0.4078 / Tsp ** 0.2466 * math.log10(Psp / 114.7))
        oilDensity = unit(oilDensity, "API", "liquid sg")
        B = 0.98496 + 1 * 10 ** -4 * (Rs ** 0.755 * gasGravityc ** 0.25 / oilDensity ** 1.5 + 0.45 * temperature) ** 1.5
    
    elif method == "Al-Marhoun 1992":          #recommended
        oilDensity = unit(oilDensity, "API", "liquid sg")
        B = 1 + 1.77342 * 10 ** -4 * Rs + 2.20163 * 10 ** -4 * Rs * gasGravity / oilDensity + 4.29258 * 10 ** -6 * Rs * (temperature - 60) * (1 - oilDensity) + 5.28707 * 10 ** -4 * (temperature - 60)
    
    elif method == "Farshad":                  #recommended
        oilDensity = unit(oilDensity, "API", "liquid sg")
        x = Rs ** 0.5956 * (gasGravity ** 0.2369 / oilDensity ** 1.3282) + 0.0976 * temperature
        B = 1 + 10 ** (-2.6541 + 0.557 * math.log10(x) + 0.3331 * math.log10(x) ** 2)
    
    elif method == "Macary":
        oilDensity = unit(oilDensity, "API", "liquid sg")
        B = (1.0031 + 0.0008 * temperature) * math.exp(0.0004 * Rs + 0.0006 * oilDensity / gasGravity)
    
    elif method == "Omar & Todd":
        x = 1.1663 + 7.62 * 10 ** -4 * oilDensity / gasGravity - 0.0339 * gasGravity
        oilDensity = unit(oilDensity, "API", "liquid sg")
        B = 0.972 + 1.47 * 10 ** -4 * (Rs * (gasGravity / oilDensity) ** 0.5 + 1.25 * temperature) ** x
    
    elif method == "Almedhaideb":
        oilDensity = unit(oilDensity, "API", "liquid sg")
        B = 1.122018 + 1.41 * 10 ** -6 * Rs * temperature / oilDensity ** 2
    
    elif method == "Elsharkawy":
        oilDensity = unit(oilDensity, "API", "liquid sg")
        B = 1 + 4.0428 * 10 ** -4 * Rs + 6.3802 * 10 ** -4 * (temperature - 60) + 7.8 * 10 ** -7 * Rs * (temperature - 60) * gasGravity / oilDensity
    
    elif method == "Khairy":
        B = 0.773413 + 7.05341 * 10 ** -4 * Rs + 0.18669 * gasGravity - 9.2589 * 10 ** -4 * oilDensity + 4.41 * 10 ** -4 * temperature
    
    elif method == "Al-Shammasi":              #recommended
        oilDensity = unit(oilDensity, "API", "liquid sg")
        B = 1 + 5.53 * 10 ** -7 * Rs * (temperature - 60) + 1.81 * 10 ** -4 * Rs / oilDensity + 4.49 * 10 ** -4 * (temperature - 60) / oilDensity + 2.06 * 10 ** -4 * Rs * gasGravity / oilDensity
    
    elif method == "Levitan & Murtha":
        oilDensity = unit(oilDensity, "API", "liquid sg")
        B = 1 + 0.0005 * Rs * (gasGravity / oilDensity) ** 0.25 + 0.2 * (temperature - 60) / (519.67 * gasGravity * oilDensity)
    
    #Velarde - requires iterations
    
    elif method == "Dindoruk & Christman":
        oilDensity = unit(oilDensity, "API", "liquid sg")
        x = (Rs ** 2.510755 / (gasGravity ** 4.852538 * oilDensity ** 11.835) + 1.365428 * 10 ** 5 * (temperature - 60) ** 2.25288 + 10.0719 * Rs) ** 0.4450849 /  (5.352624 + 2 * Rs ** -0.6309052 * gasGravity ** -0.9000749 * (temperature - 60)) ** 2
        gco = unit(oilDensity, "liquid sg", "API")
        B = 0.9871766 + 7.865146 * 10 ** -4 * x + 2.689173 * 10 ** -6 * x ** 2 + 1.100001 * 10 ** -5 * (temperature - 60) * gco / gasGravity
    
    else:
        B = "Select method"
    
    calcOilFormationFactor = B
    return calcOilFormationFactor

def calcOilViscosity(method, oilDensity, temperature, pb=None, Rs=None):     #Dead oil
#OilDensity: API  GasGravity: SG  Temperature: F    Pressure: psi
#Recommended are Beal, Beggs, Petrosky, Egbogah, Bergman-Sutton


    if method == "Andrade":  #need to research for A & B
        B = 0
        
    elif method == "Beal":
        x = math.exp(2.302585 * (0.43 + 8.33 / oilDensity))
        m = (0.32 + 1.8 * 10 ** 7 / oilDensity ** 4.53) * (360 / (temperature + 200)) ** x
    
    elif method == "Beggs & Robinson":
        x = 10 ** (3.0324 - 0.02023 * oilDensity) * temperature ** -1.163
        m = 10 ** x - 1
    
    elif method == "Glaso":
        m = (3.141 * 10 ** 10 / temperature ** 3.444) * (math.log10(oilDensity)) ** (10.313 * math.log10(temperature) - 36.447)
    
    elif method == "Labedi":
        m = 10 ** 9.224 * oilDensity ** -4.7013 * temperature ** -0.6739
    
    elif method == "Egbogah & Ng":
        x = 10 ** (1.8653 - 2.5086 * 10 ** -2 * oilDensity - 0.56411 * math.log10(temperature))
        m = 10 ** x - 1
    
    #elif method == "Al-Khafaji":      need review of source
    #    m = 10 ** (4.9563 - 0.00488 * Temperature * (OilDensity + Temperature / 30 - 14.29) ** 2.709)
    
    elif method == "Petrosky":
        m = 2.3511 * 10 ** 7 / temperature ** 2.10255 * math.log10(oilDensity) ** (4.59388 * math.log10(temperature) - 22.82792)
    
    elif method == "Kartoatmodjo & Schmidt":
        m = 1.6 * 10 ** 9 / temperature ** 2.8177 * math.log10(oilDensity) ** (5.7526 * math.log10(temperature) - 26.9718)
    
    elif method == "De Ghetto":
        if oilDensity <= 10:
            x = 10 ** (1.90296 - 1.2619 * 10 ** -2 * oilDensity - 0.61748 * math.log10(temperature))
            m = 10 ** x - 1
        elif oilDensity > 10 and oilDensity <= 22.3:
            x = 10 ** (2.06492 - 1.79 * 10 ** -2 * oilDensity - 0.70226 * math.log10(temperature))
            m = 10 ** x - 1
        elif oilDensity > 22.3 and oilDensity <= 31.1:
            m = 220.15 * 10 ** 9 / temperature ** 3.556 * math.log10(oilDensity) ** (12.5428 * math.log10(temperature) - 45.7874)
        elif oilDensity > 31.1:
            x = 10 ** (1.67083 - 1.7628 * 10 ** -2 * oilDensity - 0.61304 * math.log10(temperature))
            m = 10 ** x - 1
        
    elif method == "Agip":
        x = 10 ** (1.8513 - 2.5548 * 10 ** -2 * oilDensity - 0.56238 * math.log10(temperature))
        m = 10 ** x - 1
    
    elif method == "Fitzgerald":
        B = 0
        
    elif method == "Bennison":
        m = 10 ** (-0.8021 * oilDensity + 23.8765) * temperature ** (0.31458 * oilDensity - 9.21592)
    
    elif method == "Elsharkawy":
        x = 10 ** (2.16924 - 0.02525 * oilDensity - 0.68875 * math.log10(temperature))
        m = 10 ** x - 1
    
    elif method == "Bergman":
        x = math.exp(22.33 - 0.194 * oilDensity + 0.00033 * oilDensity ** 2 - (3.2 - 0.0185 * oilDensity) * math.log(temperature + 310))
        m = math.exp(x) - 1
    
    elif method == "Standing":
    #    OilDensity = Unit(OilDensity, "API", "liquid sg")
        x1 = 1 + 8.69 * math.log10((temperature + 459.67) / 559.67)
        x2 = 1 + 0.544 * math.log10((temperature + 459.67) / 559.67)
        x3 = -0.1285 * (2.87 * x1 - 1) * unit(oilDensity, "API", "liquid sg") / (2.87 * x1 - unit(oilDensity, "API", "liquid sg"))
        ro = unit(oilDensity, "API", "liquid sg") / (1 + 0.000321 * (temperature - 60) * 10 ** (0.00462 * oilDensity))
        kw = calcCharacterizationFactor(oilDensity)
        m = ro * 10 ** (1 / (x3 * (kw - 8.24 / unit(oilDensity, "API", "liquid sg")) + 1.639 * x2 - 1.059) - 2.17)
    
    elif method == "Dindoruk & Christman":
        a1 = 14.505357625
        a2 = -44.868655416
        a3 = 9.36579 * 10 ** 9
        a4 = -4.194017808
        a5 = -3.1461171 * 10 ** -9
        a6 = 1.517652716
        a7 = 0.010433654
        a8 = -0.00077688
        m = (a3 * temperature ** a4 * math.log10(oilDensity) ** (a1 * math.log10(temperature) + a2)) / (a5 * pb ** a6 + a7 * Rs ** a8)
    else:
        m = "Select Method"

    calcOilViscosity = m
    return calcOilViscosity

def calcBubblePointOilViscosity(method, MD, Rs):  #Gas saturated oil
#md: dead oil viscosity, cP
#Recommended Chew & Connally, Beggs & Robinson
    
    if method == "Chew & Connally":
        a = 0.2 + 0.8 / 10 ** (0.00081 * Rs)
        B = 0.43 + 0.57 / 10 ** (0.00072 * Rs)
    
    elif method == "Beggs & Robinson":
        a = 10.715 / ((Rs + 100) ** 0.515)
        B = 5.44 / ((Rs + 150) ** 0.338)
    
    calcBubblePointOilViscosity = a * MD ** B
    return calcBubblePointOilViscosity

def calcUndersaturatedOilViscosity(method, mb, pressure, pb):   #Undersaturated oil
#mb: Gas-Sat oil viscosity, cP  Pressure: pressure, psi    pb: bubble-point pressure, psi
#Recommended Beal, Kouzel, Vazquez & Beggs
    if method == "Beal":
        calcUndersaturatedOilViscosity = mb + (0.001 * (pressure - pb)) * (0.024 * mb ** 1.6 + 0.038 * mb ** 0.56)
    
    elif method == "Vazquez & Beggs":
        calcUndersaturatedOilViscosity = mb * (pressure / pb) ** (2.6 * pressure ** 1.187 * 10 ** (-3.9 * 10 ** -5 * pressure - 5))
    return calcUndersaturatedOilViscosity

def calcGasOilSurfaceTension():
    return "no data"    


def calcWaterHydrocarbonSurfaceTension():
    return "no data" 
    


def calcIsothermalCompressibility():
#Calhoun
#Standing
#Vazquez & Beggs
#Labedi
#Ahmed
#Petrosky - recommended
#Kartoatmodjo & Schmidt
#Al-Mahoun
#Farshad - recommended
#De Ghetto
#Almedhaideb
#Elsharkawy
#Dindoruk & Christman - recommended
    return "no data" 

def calcRelativePerm(Swi, Sor, krw, kro, nw, no, Sw, state):
    if state == 0:
        calcRelativePerm = krw * ((Sw - Swi) / (1 - Swi - Sor)) ** nw
    else:
        calcRelativePerm = kro * ((1 - Sor - Sw) / (1 - Swi - Sor)) ** no
    return calcRelativePerm
    

def calcGOR(gasGravity, oilDensity, pressure, temperature):
#GasGravity:SG      OilDensity:API      Pressure:psi        Temperature:F
    calcGOR = gasGravity * (pressure / 18 * 10 ** (0.0125 * oilDensity - 0.00091 * temperature)) ** 1.2048
    return calcGOR

def calcCharacterizationFactor(oilDensity):
    Mo = 6084 / (oilDensity - 5.9)
    calcCharacterizationFactor = 4.5579 * Mo ** 0.15178 * unit(oilDensity, "API", "liquid sg") ** -0.84573
    return calcCharacterizationFactor

def calcNonHCEffect(oilDensity, temperature, N2, CO2, H2S):
#OilDensity:API     Temperature:F
    pbN2 = 1.1585 + 2.86 * N2 - 1.07 * 10 ** -3
    #pbN2 = 1 + ((-2.65 * 10 ** -4 * OilDensity + 5.5 * 10 ** -3) * Temperature + (0.0931 * OilDensity - 0.8295)) * N2 + ((1.954 * 10 ** -11 * OilDensity ** 4.699) * Temperature + (0.027 * OilDensity - 2.366)) * N2 ** 2
    pbCO2 = 1 - 693.8 * CO2 * temperature ** -1.553
    pbH2S = 1 - (0.9035 + 0.0015 * oilDensity) * H2S + 0.019 * (45 - oilDensity) * H2S ** 2
    calcNonHCEffect = pbN2 * pbCO2 * pbH2S
    return calcNonHCEffect

def calcGasFormationFactor(pressure, temperature, ZFactor):
#output:cf/scf
    calcGasFormationFactor = 14.7 / unit(60, "F", "R") * ZFactor * unit(temperature, "F", "R") / pressure
    return calcGasFormationFactor

def calcGasGravity(Properties, Comp):
    i = 1
    x = 0
    for i in range(12):
        x = x + Properties[i, 1] * Comp[i]
    calcGasGravity = x / 28.97
    return calcGasGravity

def calcGasDensity(gasGravity, pressure, temperature, ZFactor):
#GasDensity:lbm/ft3      Pressure:psi    Temperature:F
    calcGasDensity = 28.967 / 10.732 * gasGravity * pressure / (ZFactor * unit(temperature, "F", "R"))
    return calcGasDensity

def calcViscosityCarr(gasGravity, pressure, temperature, N2, CO2, H2S):
#viscosity calculations by Carr - Dempsey
#temperature:degF
    temperature = unit(temperature, "F", "R")
    
    Pc = 709.604 - 58.718 * gasGravity
    Tc = 170.491 + 307.344 * gasGravity
    Tr = temperature / Tc
    pr = pressure / Pc
    
    #constants
    a0 = -2.4621182
    a1 = 2.97054714 * pr
    a2 = -0.28624054 * pr ** 2
    a3 = 0.00805422 * pr ** 3
    a4 = 2.80860949
    a5 = -3.49803305 * pr
    a6 = 0.36037302 * pr ** 2
    a7 = -0.0104432413 * pr ** 3
    a8 = -0.793385684
    a9 = 1.39643306 * pr
    a10 = -0.149144925 * pr ** 2
    a11 = 0.00441015512 * pr ** 3
    a12 = 0.0839387176
    a13 = -0.186408848 * pr
    a14 = 0.0203367881 * pr ** 2
    a15 = -0.000609579263 * pr ** 3
    
    ViscUncorrected = (temperature - 460) * (1.709 * 10 ** -5 - (2.062 * 10 ** -6) * gasGravity) + (8.188 * 10 ** -3 - (6.15 * 10 ** -3) * math.log(gasGravity))
    
    H2SCorrection = H2S * (8.49 * 10 ** -3 * math.log(gasGravity) + 3.73 * 10 ** -3)
    CO2Correction = CO2 * (9.08 * 10 ** -3 * math.log(gasGravity) + 6.24 * 10 ** -3)
    N2Correction = N2 * (8.48 * 10 ** -3 * math.log(gasGravity) + 9.59 * 10 ** -3)
    Visc1 = ViscUncorrected + H2SCorrection + CO2Correction + N2Correction
    
    viscosity = Visc1 / Tr * math.exp((a0 + a1 + a2 + a3) + Tr * (a4 + a5 + a6 + a7) + Tr ** 2 * (a8 + a9 + a10 + a11) + Tr ** 3 * (a12 + a13 + a14 + a15))
    
    #ViscosityCarr = Unit(viscosity, "cP", Munit)
    calcViscosityCarr = viscosity
    return calcViscosityCarr

def calcZFactor(gasGravity, pressure, temperature, N2, CO2, H2S):
#gasGravity:SG      Pressure:psi        Temperature:F
#Z factor calculations by Hall-Yarborough with other gases present
    temperature = unit(temperature, "F", "R")
    
    Pc = 678 - 50 * (gasGravity - 0.5) - 206.7 * N2 + 440 * CO2 + 606.7 * H2S
    Tc = 326 + 315.7 * (gasGravity - 0.5) - 240 * N2 - 83.3 * CO2 + 133.3 * H2S
    pr = pressure / Pc
    tr = Tc / temperature
    y = 0.0001
    f = 1
    i = 0
    
    while abs(f) > 0.0001 and i < 1000:    #Newton Raphson solution
        i = i + 1   #keep from getting stuck in loop
        Z = (0.06125 * pr * tr * math.exp(-1.2 * (1 - tr) ** 2)) / y
        f = -0.06125 * pr * tr * math.exp(-1.2 * (1 - tr) ** 2) + \
        (y + y ** 2 + y ** 3 - y ** 4) / (1 - y) ** 3 - \
        (14.76 * tr - 9.76 * tr ** 2 + 4.58 * tr ** 3) * y ** 2 + \
        (90.7 * tr - 242.2 * tr ** 2 + 42.4 * tr ** 3) * y ** (2.18 + 2.82 * tr)
        
        dF = (1 + 4 * y + 4 * y ** 2 - 4 * y ** 3 + y ** 4) / (1 - y) ** 4 - \
        (29.52 * tr - 19.52 * tr ** 2 + 9.16 * tr ** 3) * y + \
        (2.18 + 2.82 * tr) * (90.7 * tr - 242.2 * tr ** 2 + 42.2 * tr ** 3) * y ** (1.18 + 2.82 * tr)       
        y = y - f / dF
    calcZFactor = Z
    return calcZFactor

def calcZFactorG(gasGravity, pressure, temperature):
#GasGravity:SG      Pressure:psi        Temperature:F
#Z factor calculations by Hall-Yarborough with no additional gases
    temperature = unit(temperature, "F", "R")
    
    Pc = 709.604 - 58.718 * gasGravity
    Tc = 170.491 + 307.344 * gasGravity
    pr = pressure / Pc
    tr = Tc / temperature
    y = 0.0001
    f = 1
    i = 0
    
    while abs(f) > 0.0001 and i < 1000:    #Newton Raphson solution
        i = i + 1   #keep from getting stuck in loop
        Z = (0.06125 * pr * tr * math.exp(-1.2 * (1 - tr) ** 2)) / y
        f = -0.06125 * pr * tr * math.exp(-1.2 * (1 - tr) ** 2) + \
        (y + y ** 2 + y ** 3 - y ** 4) / (1 - y) ** 3 - \
        (14.76 * tr - 9.76 * tr ** 2 + 4.58 * tr ** 3) * y ** 2 + (90.7 * tr - 242.2 * tr ** 2 + 42.4 * tr ** 3) * y ** (2.18 + 2.82 * tr)
        
        dF = (1 + 4 * y + 4 * y ** 2 - 4 * y ** 3 + y ** 4) / (1 - y) ** 4 - \
        (29.52 * tr - 19.52 * tr ** 2 + 9.16 * tr ** 3) * y + \
        (2.18 + 2.82 * tr) * (90.7 * tr - 242.2 * tr ** 2 + 42.2 * tr ** 3) * y ** (1.18 + 2.82 * tr)     
        y = y - f / dF
    calcZFactorG = Z
    return calcZFactorG






#*************************************************************************
#OTHER EQUATIONS
#*************************************************************************
def calc_area(outer_diameter, inner_diameter):
    #outer_diameter,inner_diameter:consistent length    output:length**2
    calc_area = math.pi / 4 * (outer_diameter ** 2 - inner_diameter ** 2)
    return calc_area

def calc_fluid_compressibility(pressure, volume):
    #pressure:psi   volume:any unit     output:same as volume
    #compressibility = 293000   #15 psi, 32 degF
    #compressibility = 320000   #15 psi, 68 degF
    #compressibility = 333000   #15 psi, 120 degF
    #compressibility = 308000   #15 psi, 200 degF
    #
    #compressibility = 300000   #1500 psi, 32 degF
    #compressibility = 330000   #1500 psi, 68 degF
    #compressibility = 342000   #1500 psi, 120 degF
    #compressibility = 319000   #1500 psi, 200 degF
    #compressibility = 248000   #1500 psi, 300 degF
    #
    #compressibility = 317000   #4500 psi, 32 degF
    #compressibility = 348000   #4500 psi, 68 degF
    #compressibility = 362000   #4500 psi, 120 degF
    #compressibility = 338000   #4500 psi, 200 degF
    #compressibility = 271000   #4500 psi, 300 degF
    #
    #compressibility = 380000   #15,000 psi, 32 degF
    #compressibility = 410000   #15,000 psi, 68 degF
    #compressibility = 426000   #15,000 psi, 120 degF
    #compressibility = 405000   #15,000 psi, 200 degF
    #compressibility = 350000   #15,000 psi, 300 degF
    
    #compressibility = 312000#
    compressibility = 339000#

    calc_fluid_compressibility = pressure * volume / compressibility
    return calc_fluid_compressibility

def calc_water_hammer(od_tubing, id_tubing, solid_modulus, poissons_ratio, pipe_state, fluid_rate, fluid_density, fluid_modulus):
    fluid_velocity = calc_fluid_velocity(fluid_rate, 0, id_tubing)
    pipe_thickness = (od_tubing - id_tubing) / 2

    if id_tubing / pipe_thickness > 10:   #thin walled
        if pipe_state == 0:       #free to expand
            pipe_support_factor = 1
        elif pipe_state == 1:   #anchored upstream
            pipe_support_factor = 1 - 0.5 * poissons_ratio
        elif pipe_state == 2:   #anchored throughout
            pipe_support_factor = 1 - poissons_ratio ** 2
        elif pipe_state == 3:   #circular flow
            pipe_support_factor = 2 * pipe_thickness / id_tubing * (1 + poissons_ratio)
        
    else:                        #thick walled
        if pipe_state == 0:       #free to expand
            pipe_support_factor = 2 * pipe_thickness / id_tubing * (1 + poissons_ratio) + id_tubing / (id_tubing + pipe_thickness)
        elif pipe_state == 1:   #anchored upstream
            pipe_support_factor = 2 * pipe_thickness / id_tubing * (1 + poissons_ratio) + id_tubing / (id_tubing + pipe_thickness) * (1 - poissons_ratio / 2)
        elif pipe_state == 2:   #anchored throughout
            pipe_support_factor = 2 * pipe_thickness / id_tubing * (1 + poissons_ratio) + id_tubing / (id_tubing + pipe_thickness) * (1 - poissons_ratio ** 2)
        elif pipe_state == 3:   #circular flow
            pipe_support_factor = 2 * pipe_thickness / id_tubing * (1 + poissons_ratio)

    _CONST1 = 24.887
    _CONST2 = 1.615 * 10 ** -3
    sound_velocity = _CONST1 * math.sqrt(1 / (fluid_density * (1 / fluid_modulus + id_tubing * pipe_support_factor / (solid_modulus * pipe_thickness))))
    calc_water_hammer = _CONST2 * fluid_density * sound_velocity * fluid_velocity
    return calc_water_hammer






#*************************************************************************
#UNIT CONVERSION
#*************************************************************************
#*************************************************************************
def unit(value, x_unit, y_unit):   #value, original unit, new unit
#function will convert input unit to the top (usually SI) unit in each category, and will convert output unit from SI unit to output unit
#when units are inverted (1/x) simply switch xUnit and yUnit when calling function
    x = value   #avoids function accidentally changing the variable in other parts of the program

    #Multiplying variable by key-value to convert units to internal SI system. Inverse of dictionary will be used to convert back     
    unit_dictionary = {
        #Angle
        'degree': x,
        'radian': x / math.pi * 180,
        #Area
        'ft2': x,
        'in2': x / 144,
        'acre': x * 43560,
        #Capacity
        'm3/m': x,
        'l/m': x / 1000,
        'bbls/ft': x * 0.5216,
        'bbl/m': x * 0.159,
        'gal/ft': x * 0.0124,
        'gal/m': x * 0.0038,
        'liter/m': x / 1000,
        #Concentration
        'pptg': x,
        'gpml': x / 0.0001198,
        #Density
        'liquid sg': x,
        'sg': x,
        #Case 'gas sg'
        'ppg': x / 8.345,
        'ppga': x / 8.345,
        'API': 141.5 / (x + 131.5),
        'kg/m3': x / 1000,
        'lb/gal': x / 8.345,
        'lbs/gal': x / 8.345,
        'lb/ft3': x / 62.39431984,
        #Force
        'N': x,
        'kN': x * 1000,
        'lbf': x * 4.448,
        'klbf': x * 4.448 * 1000,
        #Length
        'm': x,
        'cm': x / 100,
        'mm': x / 1000,
        'ft': x * 0.3048,
        'in': x / 12 * 0.3048,
        'inch': x / 12 * 0.3048,
        'inches': x / 12 * 0.3048,
        #Linear mass
        'kg/m': x,
        'lbm/ft': x * 1.488,
        #Mass
        'kg': x,
        'g': x / 1000,
        'lbm': x / 2.205,
        'oz': x / 2.205 / 16,
        #Permeability
        'D': x,
        'mD': x / 1000,
        'md': x / 1000,
        #Pressure
        'kPa': x,
        'MPa': x * 1000,
        'psi': x / 0.145037744,
        'bar': x * 100,
        #Pressure Gradient
        'kPa/m': x,
        'psi/ft': x / 0.145037744 * 0.3048,
        'ppg_ECD': x * 1.1751,
        #Temperature
        'C': x,
        'K': x - 273,
        'F': (x - 32) * 5 / 9,
        'R': (x - 459.67 - 32) * 5 / 9,
        'deg C': x,
        'deg F': (x - 32) * 5 / 9,
        #Velocity
        'm/s': x,
        'ft/s': x * 0.3048,
        #Viscosity
        'cP': x,
        'poise': x * 100,
        'lbm/ft-s': x * 1488.1639,
        'lbf-s/ft2': x * 47880.259,
        'Pa-s': x * 1000,
        #Volume
        'm3': x,
        'bbl': x * 0.1589698,
        'ft3': x * 0.3048 ** 3,
        'stb': x * 0.1589698,
        'scf': x * 0.3048 ** 3,
        'Mscf': x * 0.3048 ** 3 * 1000,
        'MMscf': x * 0.3048 ** 3 * 1000000,
        'gal': x * 0.003785,
        'liter': x * 0.001,
        'Mgal': x * 3.785,
        #Volumetric Rate
        'm3/day': x,
        'm3/min': x * 1.44 * 10 ** 3,
        'bbl/day': x * 0.15897,
        'bpd': x * 0.15897,
        'bbl/min': x * 228.942,
        'bpm': x * 228.942,
        'stb/day': x * 0.15897,
        'stb/min': x * 228.942,
        'ft3/day': x * 0.028317,
        'scf/day': x * 0.028317,
        'ft3/min': x * 0.028317 * 24 * 60,
        'Mscf/day': x * 28.317,
        'MMscf/day': x * 28.317 * 10 ** 3,
        'gal/min': x * 5.451,
        'gpm': x * 5.451,
        'bbls/min': x * 228.942,
        'liter/min': x * 1.44,
    }
    
    #print(unit_dictionary[x_unit])
    #print(unit_dictionary[y_unit])
    unit = value * unit_dictionary[x_unit] / unit_dictionary[y_unit]
    return unit






#*************************************************************************
#OTHER MATH FUNCTIONS
#*************************************************************************
def calc_interpolate(x, x_vector, y_vector):

    for i in range(len(x_vector)):
        if x == x_vector[i]:     #exact solution
            y = y_vector[i]
            i = len(x_vector)
        elif x < x_vector[0]:     #if below range
            x_1 = x_vector[0]
            x_2 = x_vector[1]
            y_1 = y_vector[0]
            y_2 = y_vector[1]
            y = y_1 - (x_1 - x) / (x_2 - x_1) * (y_2 - y_1)
            i = len(x_vector)
        elif x > x_vector[len(x_vector)]:     #if above range
            x_1 = x_vector[-2]
            x_2 = x_vector[-1]
            y_1 = y_vector[-2]
            y_2 = y_vector[-1]
            y = (x - x_2) / (x_2 - x_1) * (y_2 - y_1) + y_2
            i = len(y_vector)
        elif x > x_vector[i] and x < x_vector[i + 1]:     #if between values
            x_1 = x_vector[i]
            x_2 = x_vector[i + 1]
            y_1 = y_vector[i]
            y_2 = y_vector[i + 1]
            y = (x - x_1) / (x_2 - x_1) * (y_2 - y_1) + y_1
            i = len(x_vector)  
    calcInterp = y
    """    
    for i in range(len(x_vector)):
        if x == x_vector[i]:     #exact solution
            y = y_vector[i]
            i = len(x_vector)
        elif x < x_vector[0]:     #if below range
            x_1 = x_vector[0]
            x_2 = x_vector[1]
            y_1 = y_vector[0]
            y_2 = y_vector[1]
            y = y_1 - (x_1 - x) / (x_2 - x_1) * (y_2 - y_1)
            i = len(x_vector)
        elif x > x_vector[len(x_vector)]:     #if above range
            x_1 = x_vector[len(x_vector) - 2]
            x_2 = x_vector[len(x_vector) - 1]
            y_1 = y_vector[len(y_vector) - 2]
            y_2 = y_vector[len(y_vector) - 1]
            y = (x - x_2) / (x_2 - x_1) * (y_2 - y_1) + y_2
            i = len(y_vector)
        elif x > x_vector[i] and x < x_vector[i + 1]:     #if between values
            x_1 = x_vector[i]
            x_2 = x_vector[i + 1]
            y_1 = y_vector[i]
            y_2 = y_vector[i + 1]
            y = (x - x_1) / (x_2 - x_1) * (y_2 - y_1) + y_1
            i = len(x_vector)  
    calcInterp = y
    """
    return calcInterp

# #*************************************************************************
# #from https://msdn.microsoft.com/en-us/library/csfk8t62(VS.85).aspx
# def Log10(x)   #log base 10 conversion
#     Log10 = Log(x) / Log(10)

# def Sec(x)     #Secant
#     Sec = 1 / Cos(x)

# def Cosec(x)   #Cosecant
#     Cosec = 1 / Sin(x)

# def Cotan(x)   #Cotangent
#     Cotan = 1 / Tan(x)

# def Arcsin(x)  #Inverse Sine
#     if x = 1:
#         Arcsin = pi / 2
#     elif x = -1:
#         Arcsin = -pi / 2
#     Else
#         Arcsin = Atn(x / sqrt(-x * x + 1))   

# def Arccos(x)  #Inverse Cosine
#     if x = -1:
#         Arccos = pi
#     elif x = 1:
#         Arccos = 0
#     Else
#         Arccos = Atn(-x / sqrt(-x * x + 1)) + 2 * Atn(1)   

# def Arcsec(x)  #Inverse Secant
#     Arcsec = Atn(x / sqrt(x * x - 1)) + Sgn((x) - 1) * (2 * Atn(1))

# def Arccosec(x)    #Inverse Cosecant
#     Arccosec = Atn(x / sqrt(x * x - 1)) + (Sgn(x) - 1) * (2 * Atn(1))

# def Arccotan(x)    #Inverse Cotangent
#     Arccotan = -Atn(x) + 2 * Atn(1)

# def HSin(x)    #Hyperbolic Sine
#     HSin = (Exp(x) - Exp(-x)) / 2

# def HCos(x)    #Hyperbolic Cosine
#     HCos = (Exp(x) + Exp(-x)) / 2

# def HTan(x)    #Hyperbolic Tangent
#     HTan = (Exp(x) - Exp(-x)) / (Exp(x) + Exp(-x))

# def HSec(x)    #Hyperbolic Secant
#     HSec = 2 / (Exp(x) + Exp(-x))

# def HCosec(x)  #Hyperbolic Cosecant
#     HCosec = 2 / (Exp(x) - Exp(-x))

# def HCotan(x)  #Hyperbolic Cotangent
#     HCotan = (Exp(x) + Exp(-x)) / (Exp(x) - Exp(-x))

# def HArcsin(x) #Inverse Hyperbolic Sine
#     HArcsin = Log(x + sqrt(x * x + 1))

# def HArccos(x) #Inverse Hyperbolic Cosine
#     HArccos = Log(x + sqrt(x * x - 1))

# def HArctan(x) #Inverse Hyperbolic Tangent
#     HArctan = Log((1 + x) / (1 - x)) / 2

# def HArcsec(x) #Inverse Hyperbolic Secant
#     HArcsec = Log((sqrt(-x * x + 1) + 1) / x)

# def HArccosec(x)   #Inverse Hyperbolic Cosecant
#     HArccosec = Log((Sgn(x) * sqrt(x * x + 1) + 1) / x)

# def HArccotan(x)   #Inverse Hyperbolic Cotangent
#     HArccotan = Log((x + 1) / (x - 1)) / 2

# def ATn2(ByVal DX As Double, ByVal DY As Double) As Double
#     if DY < 0:
#         ATn2 = -ATn2(DX, -DY)
#     elif DX < 0:
#         ATn2 = pi - Atn(-DY / DX)
#     elif DX > 0:
#         ATn2 = Atn(DY / DX)
#     elif DY <> 0:
#         ATn2 = pi / 2
#     Else
#         ATn2 = 1 / 0
    
# def Atn2v(ByVal DX As Double, ByVal DY As Double) As Double
#     if DY < 0:
#         Atn2v = -Atn2v(DX, -DY)
#     elif DX < 0:
#         Atn2v = pi - Atn(-DY / DX)
#     elif DX > 0:
#         Atn2v = Atn(DY / DX)
#     elif DY <> 0:
#         Atn2v = pi / 2
#     Else
#         Atn2v = CVErr(xlErrDiv0)
