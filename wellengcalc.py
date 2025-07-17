'''
@author: Jack Charles   https://jackcharlesconsulting.com/
'''

#*************************************************************************
#_PI = 3.14159265358979
GRAVITY_CONSTANT = 32.174    #gravitational constant, ft/s**2
#*************************************************************************
#VELOCITY, RHEOLOGY, FRICTION EQUATIONS
#DELTA-PRESSURE EQUATIONS
#SETTLING AND TRANSPORT EQUATIONS
#ALPHA BETA PACKING
#GRAVEL PACK SLURRY EQUATIONS
#SANDOUT EQUATIONS
#FRACTURE MECHANICS EQUATIONS
#SURVEY EQUATIONS
#TUBULAR STRESS ANALYSIS EQUATIONS
#PRODUCTION EQUATIONS
#PVT CALCULATIONS
#OTHER EQUATIONS
#MATH FUNCTIONS & USEFUL CONVERSIONS
#*************************************************************************

import math
import numpy as np
#import matplotlib.pyplot as plt
#import thermopy as tp
from util.unit import convert as ucon

#generic error handling for all functions
def err_handle(func):
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except ZeroDivisionError as e:
            print(f"Error calling {func.__name__}: {e}. Returning NaN.")
            return float('nan') # Example of returning a safe value
        except ValueError as e:
            print(f"Error calling {func.__name__}: {e}. Returning NaN.")
            return float('nan') # Example of returning a safe value
        except TypeError as e:
            print(f"Error calling {func.__name__}: {e}. Check input types.")
            raise # Re-raise the exception for further handling
        except OverflowError as e:
            print(f"Error calling {func.__name__}: {e}. Value too large.")
            return float('inf') if args[0] > 0 else float('-inf') # Example of handling overflow
        except Exception as e:
            print(f"An unexpected error occurred calling {func.__name__}: {e}")
            raise
    return wrapper


#*************************************************************************
#VELOCITY, RHEOLOGY, FRICTION EQUATIONS
#*************************************************************************
@err_handle
def calc_fluid_velocity(fluid_rate, diameter, inner_diameter=0.0):
    #flow_rate:bbl/min   diameter, inner_diameter:in   output:ft/s
    #this def only needs one diameter input, other can be 0. No preference on order
    _CONST = 13.475
    if abs(diameter ** 2 - inner_diameter ** 2) == 0: 
        fluid_velocity = 0
    else: 
        fluid_velocity = _CONST * 4 * fluid_rate / math.pi / abs(diameter ** 2 - inner_diameter ** 2)
    return fluid_velocity

def calc_shear_rate(nPrime, fluid_velocity, hydraulic_diameter, geometry: str):
    #fluid_velocity:ft/s    hydraulic_diameter:in   geometry:Pipe or Slot   output:1/s
    _CONST = 12
    if geometry == "Pipe": 
        shear_rate = ((3 * nPrime + 1) / (4 * nPrime)) * _CONST * 8 * fluid_velocity / hydraulic_diameter
    elif geometry == "Slot": 
        shear_rate = ((2 * nPrime + 1) / (3 * nPrime)) * _CONST * 6 * fluid_velocity / hydraulic_diameter
    else: 
        shear_rate = 0
    return shear_rate

def calc_power_viscosity_pipe(nPrime, kPrime, fluid_velocity, hydraulic_diameter):
    #kPrime:lbf-sn/ft2, generalized K, not geometry specific!    fluid_velocity:ft/s   hydraulic_diameter:in     output:cP
    _CONST = 12
    shear_rate = ((3 * nPrime + 1) / (4 * nPrime)) * _CONST * 8 * fluid_velocity / hydraulic_diameter
    fluid_viscosity = (kPrime * ((3 * nPrime + 1) / (4 * nPrime)) ** nPrime) * shear_rate ** (nPrime - 1)
    power_viscosity_pipe = ucon(fluid_viscosity, "lbf-sec/ft\u00b2", "cP")
    return power_viscosity_pipe

def calc_power_viscosity_slot(nPrime, kPrime, fluid_velocity, hydraulic_diameter):
    #kPrime:lbf-sn/ft2, generalized K, not geometry specific!    fluid_velocity:ft/s   hydraulic_diameter(slot width):in     output:cP
    _CONST = 12
    shear_rate = ((2 * nPrime + 1) / (3 * nPrime)) * _CONST * 6 * fluid_velocity / hydraulic_diameter
    fluid_viscosity = (kPrime * ((2 * nPrime + 1) / (3 * nPrime)) ** nPrime) * shear_rate ** (nPrime - 1)
    power_viscosity_slot = ucon(fluid_viscosity, "lbf-sec/ft\u00b2", "cP")
    return power_viscosity_slot

def calc_power_viscosity_apparent(nPrime, kPrime, shear_rate):      
    #kPrime:lbf-sn/ft2, generalized K, not geometry specific!    shearRate:1/s     output:cP
    #use with bob data
    _CONST = 1
    #fluid_viscosity = 47880 * kPrime / shearRate ** (1 - nPrime)   #no conversion needed, kPrime is geometry-dependent
    fluid_viscosity = _CONST * kPrime * ((3 * nPrime + 1) / (4 * nPrime)) ** nPrime * shear_rate ** (nPrime - 1)
    power_viscosity_apparent = ucon(fluid_viscosity, "lbf-sec/ft\u00b2", "cP")
    return power_viscosity_apparent

def calc_kPrime(nPrime, fluid_viscosity, shear_rate):
    #fluid_viscosity:cP    shear_rate:1/s     output:lbf-sn/ft2
    _CONST = 1
    kPrime = _CONST * ucon(fluid_viscosity, "cP", "lbf-sec/ft\u00b2") / (((3 * nPrime + 1) / (4 * nPrime)) ** nPrime * shear_rate ** (nPrime - 1))
    return kPrime

def calc_slurry_viscosity_keck(nPrime, fluid_velocity, hydraulic_diameter, fluid_density, solid_density, solid_loading):
    #fluid_velocity:ft/s   hydraulic_diameter:in    solid_density,fluid_density,solid_loading:ppg(a)  output:cP
    #SPE 19771 & SPE 104253 - multiply result by base fluid viscosity
    c = solid_loading / (solid_density + solid_loading)
    shear_rate = 12 * 8 * fluid_velocity / hydraulic_diameter
    relative_viscosity = (1 + (0.75 * (math.exp(1.5 * nPrime) - 1) * math.exp(-(1 - nPrime) * shear_rate / 1000)) * (1.25 * c / (1 - 1.5 * c))) ** 2
    relative_density = (1 + solid_loading / fluid_density) / (1 + solid_loading / solid_density)
    slurry_viscosity_keck = relative_viscosity ** 0.55 * relative_density ** 0.45
        
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
    #ff = 1 / 4 / (-2 * log10(e / 3.7065 - 5.0452 / NRe * log10(e ** 1.1098 / 2.8257 + 5.8506 / (NRe ** 0.8981)))) ** 2
    if NRe < 2100:
        friction_chen = 16 / NRe
    else:
        _e = roughness / hydraulic_diameter
        friction_chen = 1 / (-4 * math.log10(_e / 3.7065 - 5.0452 / NRe * math.log10(_e ** 1.1098 / 2.8257 + (7.149 / NRe) ** 0.8981))) ** 2
    return friction_chen


def calc_friction_colebrook(hydraulic_diameter, NRe, roughness):    
    #hydraulic_diameter, roughness:in
    #Fanning friction factor, Colebrook-White implicit solution. Note Fanning = 4xDarcy friction
    #ff = 1 / (f ** 0.5) + 4 * math.log10(e / 3.7065 + 1.2613 / (NRe * (f ** 0.5)))
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
            #ff1 = 1 / (f0 ** 0.5) + 2 * log10(_e / 3.7065 + 2.51 / (NRe * (f0 ** 0.5)))                    #Darcy friction factor
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
    
def calc_standoff(tbg_od, csg_id, centralizer_od, blade_count):
    angle = math.pi / blade_count
    standoff1 = (centralizer_od * math.cos(angle) - tbg_od) / 2
    standoff2 = csg_id / 2 - ((csg_id / 2) ** 2 - (centralizer_od / 2 * math.sin(angle)) ** 2) ** 0.5
    #alternative form of standoff2. http://mathworld.wolfram.com/CircularSegment.html
    #standoff2 = (2 * csg_id / 2 - ((2 * csg_id / 2) ** 2 - 4 * (centralizer_od / 2 * math.sin(angle)) ** 2) ** 0.5) / 2
    standoff = standoff1 + standoff2
    return standoff

def calc_eccentricity(tbg_od, csg_id, standoff_radius):
    #SPE 111514 - Haciislamoglu & Langlinais solution Eqn 1, 1990
    #range from 0 (fully concentric) to 1 (fully eccentric)
    #standoff_radius is used instead of centralizer_od because the standoff on a bladed centralizer will be less. See calc_standoff
    if standoff_radius < 0:
        eccentricity = 1
    else: 
        eccentricity = (csg_id / 2 - (standoff_radius + tbg_od / 2)) / (csg_id / 2 - tbg_od / 2)
        #eccentricity = (csg_id - centralizer_od / (csg_id - tbg_od)
    return eccentricity

def calc_eccentricity_factor(tbg_od, csg_id, eccentricity):
    #SPE 31147 - Hang       Two different solutions for eccentricity are provided in the versions of the papers
    #This is multiplied by hydraulic diameter, per paper Deff = S*Dh, which is then used to determine NRe and dP
    if eccentricity == 0:
        eccentricity_factor = 0.66667
    else:
        #eccentricity_factor = 4.33 + (tbg_od / csg_id) * (10.54 * tbg_od / csg_id - 12.26)
        eccentricity_factor = 1.004 + 1.873 * (tbg_od / csg_id) - 1.826 * (tbg_od / csg_id) ** 2 + 0.617 * (tbg_od / csg_id) ** 3
    return eccentricity_factor
    
def calc_eccentricity_factor_powerlaw(nPrime, NRe, tbg_od, csg_id, eccentricity):  
    #this is multiplied by dP of a concentric annulus. Can be valid for nPrime = 1
    diameter_ratio = tbg_od / csg_id
    if NRe < (3250 - 1150 * nPrime):    
        #SPE 111514 - Haciislamoglu & Langlinais solution Eqn 2, 1990
        #typically valid for n > 0.4
        eccentricity_factor_powerlaw = 1 - (0.072 * eccentricity / nPrime * diameter_ratio ** 0.8454) - (1.5 * eccentricity ** 2 * math.sqrt(nPrime) * diameter_ratio ** 0.1852) + (0.96 * eccentricity ** 3 * math.sqrt(nPrime) * diameter_ratio ** 0.2527)
        #Eqn 3 - valid 0.2 < diameter_ratio < 0.8, 0 < ecc < 0.98, 0.2 < nPrime < 1.0
        #eccentricity_factor_powerlaw = 1 - (0.1019 * eccentricity * nPrime * diameter_ratio ** -0.4675) - (1.6152 * eccentricity ** 2 * nPrime ** 0.085 * diameter_ratio ** 0.7875) + (1.1434 * eccentricity ** 3 * nPrime ** 0.0547 * diameter_ratio ** 1.1655)
    else:
        #also described in Advanced Drilling & Well Technology - Aadnoy, Cooper, Miska, Mitchel, & Payne, 2009. Basis for split between laminar and turbulent flow
        eccentricity_factor_powerlaw = 1 - (0.048 * eccentricity / nPrime * diameter_ratio ** 0.8454) - (0.67 * eccentricity ** 2 * math.sqrt(nPrime) * diameter_ratio ** 0.1852) + (0.28 * eccentricity ** 3 * math.sqrt(nPrime) * diameter_ratio ** 0.2527) 
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
    #pressure drop due to Fanning friction
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
    
    _CONST = 0.01938         #constant defined here to avoid interference with other functions above
    DPf_full = _CONST * 2 * fluid_friction * fluid_density * length * fluid_velocity ** 2 / hydraulic_diameter
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
    
    settling_velocity = (4 / 3 * GRAVITY_CONSTANT * hydraulic_diameter / Cd * (solid_density - fluid_density) / fluid_density) ** 0.5
    return settling_velocity

def calc_settling_Stokes(solid_density, fluid_density, hydraulic_diameter, fluid_viscosity):   #solution given by Stoke#s equation, Cd = 24/Re
    #solid_density,fluid_density:ppg      hydraulic_diameter:in     fluid_viscosity:cP    output:ft/s
    _CONST = 77.307
    settling_Stokes = _CONST * ((solid_density - fluid_density) * GRAVITY_CONSTANT * hydraulic_diameter ** 2) / (18 * fluid_viscosity)
    return settling_Stokes

def calc_settling_intermediate(solid_density, fluid_density, hydraulic_diameter, fluid_viscosity): #solution given by Cd = 18*Re**-0.6
    #solid_density,fluid_density:ppg      hydraulic_diameter:in     fluid_viscosity:cP    output:ft/s
    _CONST = 3.169
    settling_intermediate = _CONST * (2 * GRAVITY_CONSTANT / 27 * (solid_density - fluid_density) / fluid_density) ** (5 / 7) * (hydraulic_diameter ** (8 / 7)) / ((fluid_viscosity / fluid_density) ** (3 / 7))
    return settling_intermediate

def calc_settling_Newtonian(solid_density, fluid_density, hydraulic_diameter):       #solution given by Cd = 0.44
    #solid_density,fluid_density:ppg      hydraulic_diameter:in     output:ft/s
    _CONST = 0.289
    f = 0.44            #for Re > 1000
    settling_Newtonian = _CONST * (3 * GRAVITY_CONSTANT * hydraulic_diameter * (solid_density - fluid_density) / fluid_density) ** 0.5
    return settling_Newtonian

def calc_settling_power_law(solid_density, fluid_density, hydraulic_diameter, nPrime, kPrime):    #power law fluid settling
    #solid_density,fluid_density:ppg      hydraulic_diameter:in   output:ft/s
    #settling_power_law = (0.8667 * (solid_density - fluid_density) / kPrime) * (pipeDiameter ** (1 / nPrime) * (2 * nPrime + 1) * solid_density) / (108 * nPrime)       #Halliburton CT manual
    #settling_power_law = (((solid_density - fluid_density) * GRAVITY_CONSTANT * hydraulic_diameter ** (nPrime + 1)) / (2 * 3 ** (nPrime + 1) * kPrime)) ** (1 / nPrime)    #Economides Modern Hydraulic Fracturing
    settling_power_law = (((solid_density - fluid_density) * GRAVITY_CONSTANT * hydraulic_diameter ** (nPrime + 1)) / (18 * 3 ** (nPrime - 1) * kPrime)) ** (1 / nPrime)    #Economides Reservoir Stimulation
    return settling_power_law

def calc_settling_Budryck(solid_density, fluid_density, hydraulic_diameter): #intermediate solution with Budryck equation
    #solid_density,fluid_density:ppg      hydraulic_diameter:in     output:ft/s
    _CONST = 1 / 304.8
    settling_Budryck = _CONST * 8.925 / (hydraulic_diameter * 25.4) * (math.sqrt(1 + 95 * (hydraulic_diameter * 25.4) ** 3 * (solid_density - fluid_density) / fluid_density) - 1)
    #settling_Budryck = _CONST * 8.925 / (fydraulicDiameter) * (sqrt(1 + 95 * (hydraulic_diameter) ** 3 * (solid_density - fluid_density) / fluid_density) - 1)
    return settling_Budryck

def calc_settling_Rittinger(solid_density, fluid_density, hydraulic_diameter): #turbulent solution with Rittinger equation
    #solid_density,fluid_density:ppg      hydraulic_diameter:in     output:ft/s
    _CONST = 1 / 304.8
    settling_Rittinger = _CONST * 87 * ((solid_density - fluid_density) / fluid_density * hydraulic_diameter * 25.4) ** 0.5
    return settling_Rittinger

def calc_settling_McCabeSmith(solid_density, fluid_density, hydraulic_diameter, fluid_viscosity): #intermediate solution with McCabe-Smith equation, 2<Re<500. SPE 187498
    #solid_density,fluid_density:ppg      hydraulic_diameter:in     output:ft/min
    _CONST = 1 / 60
    settling_McCabeSmith = _CONST * (0.072 * GRAVITY_CONSTANT * ucon(hydraulic_diameter, "m", "in") ** 1.6 * (solid_density - fluid_density) / (fluid_density ** 0.4 * fluid_viscosity ** 0.6)) ** 0.71
    return settling_McCabeSmith

def calc_settling_turbulent(solid_density, fluid_density, hydraulic_diameter): 
    #turbulent solution with turbulent equation, 500<Re. SPE 187498
    #solid_density,fluid_density:ppg      hydraulic_diameter:in     output:ft/min
    _CONST = 1 / 60
    settling_turbulent = _CONST * 1.74 * (GRAVITY_CONSTANT * (solid_density - fluid_density) / fluid_density * ucon(hydraulic_diameter, "m", "in")) ** 0.5
    return settling_turbulent

def calc_void_fraction(bulk_density, solid_density):
    #bulk_density,solid_density in same units
    void_fraction = 1 - bulk_density / solid_density
    return void_fraction

def calc_solids_fraction(solid_loading, solid_density):
    #also known as solids concentration, c
    solids_fraction = solid_loading / (solid_loading + solid_density)
    return solids_fraction

def calc_fluidization(hydraulic_diameter, solid_density, fluid_density, fluid_viscosity, porosity, sphericity):
    #hydraulic_diameter:in   solid_density,fluid_density:ppg   fluid_viscosity:cP    porosity,sphericity:fraction    output:ft/s
    #only good for small Re<10
    #according to McCabe, 30x fluidization velocity is ok to produce at
    _CONST = 77.307
    fluidization = _CONST * ((solid_density - fluid_density) * GRAVITY_CONSTANT * hydraulic_diameter ** 2) / (150 * fluid_viscosity) * (porosity ** 3 * sphericity ** 2 / (1 - porosity))
    return fluidization

def calc_angle_correction(inclination):
    #angle_correction = 1 + 2 * inclination / 45                                        #Rubiandini model, for transport UP
    #angle_correction = 0.0342 * inclination - 0.000233 * inclination ** 2 - 0.213      #SPE 25872 Larsen, for transport UP
    angle_correction = 0.022 * inclination - 0.9793                                     #Hang model
    return angle_correction

def calc_horizontal_transport_Oroskar(csg_id, solid_diameter, solid_density, fluid_density, fluid_viscosity, c, 
                                      x = 0.95, _CONSTY = 1.85,  _CONSTn1 = 0.1536,  _CONSTn2 = 0.3564,  _CONSTn3 = -0.378 ,  _CONSTn4 = 0.09,  _CONSTn5 = 0.3):
    #solid_density,fluid_density:ppg      csg_id,solid_diameter:in     fluid_viscosity:cP    output:ft/s
    #c: solids concentration, loading/(loading+density)
    #for non-concentric flow, however can replace csg_id with hydraulic_diameter or equivalent_diameter of flow area for interesting results
    _CONST1 = math.sqrt(1/12)   #constant, sqrt(1/12)
    _CONST2 = 927.6866          #constant for modified NRe to multiply by ft/s. See calc_NRe_newton function

    horizontal_transport_Oroskar = _CONST1 * math.sqrt(GRAVITY_CONSTANT * solid_diameter * (solid_density / fluid_density - 1)) \
        * _CONSTY * c ** _CONSTn1 \
        * (1 - c) ** _CONSTn2 \
        * (solid_diameter / csg_id) ** _CONSTn3 \
        * (_CONST1 * _CONST2 * csg_id * fluid_density / fluid_viscosity * math.sqrt(GRAVITY_CONSTANT * solid_diameter * (solid_density / fluid_density - 1))) ** _CONSTn4 \
        * x ** _CONSTn5
    return horizontal_transport_Oroskar

def calc_horizontal_transport_OroskarMod(hydraulic_diameter, solid_diameter, solid_density, fluid_density, fluid_viscosity, c, 
                                      x = 1, _CONSTY = 3.52,  _CONSTn1 = -0.111,  _CONSTn2 = -2.97,  _CONSTn3 = -0.357 ,  _CONSTn4 = -0.0595,  _CONSTn5 = 0.3):
    #solid_density,fluid_density:ppg      hydraulic_diameter,solid_diameter:in     fluid_viscosity:cP    output:ft/s
    #c: solids concentration, loading/(loading+density)
    #modified Oroskar for concentric flow, updated constants from regression fit, uncomment to compare
    _CONST1 = math.sqrt(1/12)   #constant, sqrt(1/12)
    _CONST2 = 927.6866          #constant for modified NRe to multiply by ft/s. See calc_NRe_newton function
            
    horizontal_transport_OroskarMod = _CONST1 * math.sqrt(GRAVITY_CONSTANT * solid_diameter * (solid_density / fluid_density - 1)) \
        * _CONSTY * c ** _CONSTn1 \
        * (1 - c) ** _CONSTn2 \
        * (solid_diameter / hydraulic_diameter) ** _CONSTn3 \
        * (_CONST1 * _CONST2 * hydraulic_diameter * fluid_density / fluid_viscosity * math.sqrt(GRAVITY_CONSTANT * solid_diameter * (solid_density / fluid_density - 1))) ** _CONSTn4 \
        * x ** _CONSTn5
    return horizontal_transport_OroskarMod

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
#ALPHA BETA PACKING
#*************************************************************************
#Alpha-Beta calculations, Hang solution
#csg_id may be replaced with openhole ID as necessary
def calc_Cfunction(csg_id, dune_height):
    #csg_id,dune_height:ft      output:ft2
    Cfunction = 2 * math.sqrt(dune_height * (csg_id - dune_height))
    return Cfunction

def calc_Sfunction(csg_id, dune_height):
    #csg_id,dune_height:ft      output:ft2
    Sfunction = csg_id * math.acos((2 * dune_height - csg_id) / csg_id)
    return Sfunction

def calc_Gfunction(csg_id, dune_height):
    #csg_id,dune_height:ft      output:ft2
    Gfunction = 0.25 * (csg_id * calc_Sfunction(csg_id, dune_height) - (2 * dune_height - csg_id) * calc_Cfunction(csg_id, dune_height))
    return Gfunction

def calc_Pfunction(csg_id, dune_height, H1, H2, screen_OD):
    #csg_id,dune_height,H1,H2,inner_diameter:in      output:in
    if dune_height < H1: 
        Pfunction = calc_Cfunction(csg_id, dune_height) + calc_Sfunction(csg_id, dune_height) + math.pi * screen_OD
    elif dune_height > H2: 
        Pfunction = calc_Cfunction(csg_id, dune_height) + calc_Sfunction(csg_id, dune_height)
    else: 
        Pfunction = calc_Cfunction(csg_id, dune_height) - calc_Cfunction(screen_OD, dune_height - H1) + calc_Sfunction(csg_id, dune_height) + calc_Sfunction(screen_OD, dune_height - H1)
    return Pfunction

def calc_Afunction(csg_id, dune_height, H1, H2, screen_OD):
    #csg_id,dune_height,H1,H2,inner_diameter:in    output:in2
    if dune_height < H1:
        Afunction = calc_Gfunction(csg_id, dune_height) - math.pi / 4 * screen_OD ** 2
    elif dune_height > H2:
        Afunction = calc_Gfunction(csg_id, dune_height)
    else:
        Afunction = calc_Gfunction(csg_id, dune_height) - calc_Gfunction(screen_OD, dune_height - H1)
    return Afunction 

def calc_alphawave_DH_Hang(csg_id, screen_OD, centralizer_OD, dune_height_ratio):
    #csg_id, screen_OD, centralizer_OD:in    DHR:dimensionless       output:in
    _H1 = (centralizer_OD - screen_OD) / 2
    _H2 = _H1 + screen_OD
    dune_height = dune_height_ratio * csg_id
    alphawave_DH_Hang = 4 * calc_Afunction(csg_id, dune_height, _H1, _H2, screen_OD) / calc_Pfunction(csg_id, dune_height, _H1, _H2, screen_OD)
    return alphawave_DH_Hang

#****************************************************************
def calc_tbg_dune_height(diameter, dune_height_ratio):     
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

    return theta * 180 / math.pi, area, perimeter, hydraulic_diameter, equivalent_diameter

#Alpha-Beta Calculations
def calc_alphawave_dune_height(csg_id, screen_od, centralizer_od, dune_height_ratio):
    #csg_id, screen_od, centralizer_od,output:in
    h = dune_height_ratio * csg_id
    g_i = min((csg_id - screen_od), (centralizer_od - screen_od)) / 2
    
    h_o = min(max(h, 0), csg_id)
    h_i = min(max(h - g_i, 0), screen_od)
    
    theta_o = 2 * math.asin((csg_id - 2 * h_o) / csg_id) + math.pi
    theta_i = 2 * math.asin((screen_od - 2 * h_i) / screen_od) + math.pi

    area_o = csg_id ** 2 / 8 * (theta_o + math.sin(theta_o - math.pi))
    area_i = screen_od ** 2 / 8 * (theta_i + math.sin(theta_i - math.pi))

    perimeter_o = csg_id * theta_o / 2 + csg_id * math.cos((theta_o - math.pi) / 2)
    perimeter_i = screen_od * theta_i / 2 - screen_od * math.cos((theta_i - math.pi) / 2)

    width_o = csg_id * math.cos((theta_o - math.pi) / 2)
    width_i = screen_od * math.cos((theta_i - math.pi) / 2)

    hydraulic_diameter = 4 * (area_o - area_i) / (perimeter_o + perimeter_i)
    equivalent_diameter = math.sqrt(4 * (area_o - area_i) / math.pi)

    return hydraulic_diameter, equivalent_diameter, area_o, area_i, perimeter_o, perimeter_i, width_o, width_i





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
        output.append([proppant_concentration[i],sand_mass[i] / proppant_concentration[i] * (1 + proppant_concentration[i] * solid_absVol)])
    Nolte_schedule = output[:]
    return Nolte_schedule






#*************************************************************************
#SANDOUT EQUATIONS
#*************************************************************************
def calc_Darcy_sandout(csg_id, tbg_od, fluid_rate, fluid_viscosity, permeability, height):
    #csg_id, tbg_od:in     flowRate:bbl/min   fluid_viscosity:cP    permeability:D  height:ft   output:psi
    area = ucon(calc_area(csg_id, tbg_od), "in\u00b2", "ft\u00b2")
    Darcy_sandout = height * fluid_viscosity * fluid_rate / (0.000783 * permeability * area)
    return Darcy_sandout

def calc_Forcheimer_sandout(csg_id, tbg_od, fluid_rate, fluid_density, fluid_viscosity, permeability, height):
    #csg_id, tbg_od:in     fluid_rate:bbl/min   fluid_density:ppg  fluid_viscosity:cP    permeability:D  height:ft  output:psi
    area = ucon(calc_area(csg_id, tbg_od), "in\u00b2", "ft\u00b2")
    #Forcheimer_sandout = height * (1279 * fluid_viscosity * fluid_rate * area ** 2 + 2.82 * fluid_density * fluid_rate ** 2 * permeability ** 0.698) / (permeability * area ** 2)
    Forcheimer_sandout = height * (1279 * fluid_viscosity * fluid_rate * area + 2.82 * fluid_density * fluid_rate ** 2 * permeability ** 0.698) / (permeability * area ** 2)
    #Forcheimer_sandout = height * (1279 * fluid_viscosity * fluid_rate * area + 4.63 * fluid_density * fluid_rate ** 2 * permeability ** 0.45) / (permeability * area ** 2)
    return Forcheimer_sandout

def calc_Ergun_sandout(csg_id, tbg_od, solid_diameter, fluid_rate, fluid_viscosity, fluid_density, height, void_fraction, sphericity):
    #csg_id, tbg_od:in     fluid_rate:bbl/min   fluid_viscosity:cP    fluid_density:ppg  height:ft   output:psi
    _CONST1 = 1.954 * 10 ** -6
    _CONST2 = 1.697 * 10 ** -4
    area = ucon(calc_area(csg_id, tbg_od), "in\u00b2", "ft\u00b2")
    Ergun_sandout = height * (_CONST1 * 150 * (fluid_rate / area * fluid_viscosity * (1 - void_fraction) ** 2) / ((sphericity * solid_diameter) ** 2 * void_fraction ** 3) \
                              + (_CONST2 * 1.75 * (fluid_density * (fluid_rate / area) ** 2) / (sphericity * solid_diameter) * (1 - void_fraction) / void_fraction))
    return Ergun_sandout

def calc_sand_fill_rate(csg_id, tbg_od, fluid_rate, solid_loading, absVol, bulk_density, proppant_perf_weight):
#csg_id,tbg_od:in   fluid_rate:bpm    solid_loading:ppga    absVol:gal/lbm   bulk_density:lbm/ft3  proppant_perf_weight:lbm/ft in perfs
    area = ucon(calc_area(csg_id, tbg_od), "in\u00b2", "ft\u00b2")
    sand_rate = calc_proppant_mass(fluid_rate * 42, absVol, solid_loading)
    sand_volume = proppant_perf_weight + bulk_density * area
    sand_fill_rate = sand_rate / sand_volume
    return sand_fill_rate

def calc_settling_factor(solid_loading, bulk_density, absVol):   #BHI tech facts
    #solid_loading:ppg    bulk_density:lbm/ft3  absVol:gal/lbm   output:%
    #this gives the bulk settling factor. this is distinct from the volume concentration c, in that this uses (1 / bulkDensity)
    #refer to calc_solids_concentration
    _CONST = 7.4608
    settling_factor = _CONST * solid_loading * (1 / bulk_density) / (1 + absVol * solid_loading) * 100
    return settling_factor




#*************************************************************************
#FRACTURE MECHANICS EQUATIONS
#*************************************************************************
def calc_perm_estimate(solid_diameter, porosity):  #Blake-Kozeny, SPE 31141
    #solid_diameter:mm
    perm_estimate = solid_diameter ** 2 * porosity ** 3 / (150 * (1 - porosity) ** 2)
    return perm_estimate

def calc_fracture_rate(permeability, formation_height, frac_gradient, formation_TVD, fluid_density, fluid_viscosity, wellbore_ID, Bo, skin):
    #permeability:mD    formation_height,formation_TVD:ft  frac_gradient:psi/ft    fluid_density:ppg    fluid_viscosity:cP   wellbore_ID:in   Bo:stb/rb
    _CONST = 1 / 24 / 60
    fracture_rate = _CONST * calc_Darcy_IPR(permeability, formation_height, frac_gradient * formation_TVD, 0, fluid_density * formation_TVD * 0.052, Bo, fluid_viscosity, calc_shape_radial(500 * 12, wellbore_ID / 2, skin, "Pseudoradial"))
    return fracture_rate

def calc_injectivity_index(rate, pressure):
    #rate:bpm   pressure:psi    output:bpd/psi
    #pressure should be total pressure over BHP, e.g. surface annulus + overbalance
    _CONST = 1440
    injectivity_index = rate * _CONST / pressure
    return injectivity_index

def calc_Eaton_closure(overburden_pressure, pore_pressure, Poissons_ratio, biot_number):
    #overburden_pressure,pore_pressure:psi    output:psi      #psi/ft can be used also but must be consistent!
    Eaton_closure = Poissons_ratio / (1 - Poissons_ratio) * (overburden_pressure - biot_number * pore_pressure) + biot_number * pore_pressure
    return Eaton_closure
    
def calc_Eaton_frac_impermeable(overburden_pressure, pore_pressure, Poissons_ratio, biot_number):
    #overburden_pressure,pore_pressure:psi    output:psi      #psi/ft can be used also but must be consistent!
    Eaton_frac_impermeable = 2 * Poissons_ratio / (1 - Poissons_ratio) * (overburden_pressure - biot_number * pore_pressure) + biot_number * pore_pressure
    return Eaton_frac_impermeable

def calc_Eaton_frac_permeable(overburden_pressure, pore_pressure, Poissons_ratio, biot_number):
    #overburden_pressure,pore_pressure:psi    output:psi      #psi/ft can be used also but must be consistent!
    Eaton_frac_permeable = 2 * Poissons_ratio * (overburden_pressure - biot_number * pore_pressure) + biot_number * pore_pressure
    return Eaton_frac_permeable

def calc_SGS_closure(overburden_pressure, pore_pressure):
    #overburden_pressure,pore_pressure:psi
    SGS_closure = 0.488 * (overburden_pressure - pore_pressure) + pore_pressure
    return SGS_closure

def calc_SGS_frac(overburden_pressure, pore_pressure):
    #rden_pressure,pore_pressure:psi
    SGS_frac = 0.6934 * (overburden_pressure - pore_pressure) + pore_pressure
    return SGS_frac

def calc_overburden(seawater_density, rock_density, liquid_density, formation_depth, seabed_depth, phi_0):    #from Applied Drilling Engineering
    #seawater_density,rock_density,liquid_density:ppg     formation_depth,seabed_depth:ft     output:psi
    #typical constants: phi_0~0.40  seawater_density~1.029SG     rock_density~2.6SG   liquid_density~1.074SG
    #equation will provide overburden pressure in offshore environments. the resultant gradient by dividing by TVD will give gradient to surface
    bed_thickness = formation_depth - seabed_depth

    if bed_thickness // 1000 == 0:     #calculate to nearest rounded thousand
        bulk_density = 1.95
        phi = 0.43
    elif bed_thickness // 1000 == 1:
        bulk_density = 2.02
        phi = 0.38
    elif bed_thickness // 1000 == 2:
        bulk_density = 2.06
        phi = 0.35
    elif bed_thickness // 1000 == 3:
        bulk_density = 2.11
        phi = 0.32
    elif bed_thickness // 1000 == 4:
        bulk_density = 2.16
        phi = 0.29
    elif bed_thickness // 1000 == 5:
        bulk_density = 2.19
        phi = 0.27
    elif bed_thickness // 1000 == 6:
        bulk_density = 2.24
        phi = 0.24
    elif bed_thickness // 1000 == 7:
        bulk_density = 2.27
        phi = 0.22
    elif bed_thickness // 1000 == 8:
        bulk_density = 2.29
        phi = 0.2
    elif bed_thickness // 1000 == 9:
        bulk_density = 2.33
        phi = 0.18
    elif bed_thickness // 1000 == 10:
        bulk_density = 2.35
        phi = 0.16
    elif bed_thickness // 1000 == 11:
        bulk_density = 2.37
        phi = 0.15
    elif bed_thickness // 1000 == 12:
        bulk_density = 2.38
        phi = 0.14
    elif bed_thickness // 1000 == 13:
        bulk_density = 2.4
        phi = 0.13
    elif bed_thickness // 1000 == 14:
        bulk_density = 2.41
        phi = 0.12
    elif bed_thickness // 1000 == 15:
        bulk_density = 2.43
        phi = 0.11
    elif bed_thickness // 1000 == 16:
        bulk_density = 2.44
        phi = 0.1
    elif bed_thickness // 1000 == 17:
        bulk_density = 2.45
        phi = 0.098
    elif bed_thickness // 1000 == 18:
        bulk_density = 2.46
        phi = 0.092
    elif bed_thickness // 1000 == 19:
        bulk_density = 2.47
        phi = 0.085
    elif bed_thickness // 1000 == 20:
        bulk_density = 2.48
        phi = 0.079
    else:
        bulk_density = -1 * 10 ** 9 * (bed_thickness) ** 2 + 5 * 10 ** -5 * (bed_thickness) + 1.9676  #from polynomial curve fit
        phi = 9 * 10 ** -10 * (bed_thickness) ** 2 - 3 * 10 ** -5 * (bed_thickness) + 0.4164         #from polynomial curve fit
    
    if phi == 0 or phi_0 == 0:
        k = 0
    else:
        k = math.log(phi_0 / phi) / bed_thickness
    
    overburden = 0.052 * (seawater_density * seabed_depth + rock_density * bed_thickness) - (rock_density - liquid_density) * 0.052 * phi_0 / k * (1 - math.exp(-k * bed_thickness))
    return overburden

def calc_formation_mean_stress(overburden_pressure, closure_pressure, pore_pressure):      #from GWong
    #overburden_pressure,closure_pressure,pore_pressure:psi/ft    output:psi
    formation_mean_stress = (overburden_pressure + 2 * closure_pressure) / 3 - pore_pressure
    return formation_mean_stress

def calc_formation_Youngs_modulus(mean_stress):      #from GWong
    #mean_stress:psi     output:psi
    formation_Youngs_modulus = 85.053 * mean_stress + 515477
    return formation_Youngs_modulus

def calc_formation_toughness(Youngs_modulus, Poissons_ratio, gamma):
    #Youngs_modulus:psi      gamma:in*lbm/in2, estimates are 0.5-2.0      output:psi*sqrt(in)
    formation_toughness = math.sqrt(gamma * Youngs_modulus / (1 - Poissons_ratio ** 2))
    return formation_toughness

def calc_elastic_shear_modulus(Youngs_modulus, Poissons_ratio):      #from Economides, Petroleum Production Systems
    elastic_shear_modulus = Youngs_modulus / (2 * (1 + Poissons_ratio))
    return elastic_shear_modulus

def calc_plane_strain_modulus(Youngs_modulus, Poissons_ratio):      #from Economides, Petroleum Production Systems
    plane_strain_modulus = Youngs_modulus / (1 - Poissons_ratio ** 2)
    return plane_strain_modulus

def calc_PKN_width(fluid_rate, fluid_viscosity, Poissons_ratio, half_length, Youngs_modulus, width_type: str):     #from Economides, Petroleum Production Systems
    #fluid_rate:bpm   fluid_viscosity:cP    halfLength:ft   Youngs_modulus:psi   width_type:max or average    #output:in
    _CONST = 0.13
    if width_type == "max": gamma = 1
    else: gamma = math.pi / 4 * 0.75
    
    PKN_width = _CONST * 2.31 * (fluid_rate * fluid_viscosity * 2 * (1 - Poissons_ratio ** 2) * half_length / Youngs_modulus) ** 0.25 * gamma
    return PKN_width

def calc_PKN_width_powerlaw(fluid_rate, nPrime, kPrime, half_length, frac_height, Youngs_modulus, width_type: str):     #from Economides, Petroleum Production Systems
    #fluid_rate:bpm   half_length, frac_height:ft   Youngs_modulus:psi   width_type:max or average    #output:in
    _CONST = 1
    if width_type == "max": gamma = 1
    else: gamma = math.pi / 4 * 0.75
    
    PKN_width_powerlaw = _CONST * 12 * ((128 / 3 / math.pi) * (nPrime + 1) * ((2 * nPrime + 1) / nPrime) ** nPrime * (0.9775 / 144) * (5.61 / 60) ** nPrime) ** (1 / (2 * nPrime + 2)) * \
    (fluid_rate ** nPrime * kPrime * half_length * frac_height ** (1 - nPrime) / Youngs_modulus) ** (1 / (2 * nPrime + 2)) * gamma
    return PKN_width_powerlaw
    
def calc_KGD_width(fluid_rate, fluid_viscosity, Poissons_ratio, half_length, Youngs_modulus, frac_height, width_type: str):     #from Economides, Petroleum Production Systems
    #fluid_rate:bpm   fluid_viscosity:cP    halfLength, frac_height:ft   Youngs_modulus:psi     widthType:max or average    #output:in
    _CONST = 0.128
    if width_type == "max": gamma = 1
    else: gamma = math.pi / 4
    
    KGD_width = _CONST * 2.27 * (fluid_rate * fluid_viscosity * 2 * (1 - Poissons_ratio ** 2) * half_length ** 2 / (Youngs_modulus * frac_height)) ** 0.25 * gamma
    return KGD_width

def calc_penny_width(fluid_rate, fluid_viscosity, Poissons_ratio, frac_radius, Youngs_modulus, width_type: str):       #from Guo & Ghalambor
    #fluid_rate:bpm   fluid_viscosity:cP    frac_radius:ft   Youngs_modulus:psi   WidthType:max or average    #output:in
    if width_type == "max": _CONST = 0.85
    else: _CONST = 2.56 / 12
    
    penny_width = _CONST * (fluid_viscosity * fluid_rate * (1 - Poissons_ratio) * frac_radius / Youngs_modulus) ** 0.25
    return penny_width
    
def calc_net_pressure_width(Youngs_modulus, Poissons_ratio, characteristic_length, net_pressure, model: str, width_type: str): #from Economides, Reservoir Stimulation
    #Youngs_modulus,net_pressure:psi  characteristic_length:ft     output:in
    _CONST = 1 / 12
    
    if width_type == "max": gamma = 1
    else: gamma = math.pi / 4 * 0.75
    
    plane_strain_modulus = calc_plane_strain_modulus(Youngs_modulus, Poissons_ratio)
    
    if model == "PKN": fracture_stiffness = _CONST * 2 * plane_strain_modulus / math.pi / characteristic_length      #PKN, characteristic_length is height
    elif model == "KGD": fracture_stiffness = _CONST * plane_strain_modulus / math.pi / characteristic_length         #KGD, characteristic_length is length
    elif model == "Penny": fracture_stiffness = _CONST * 3 * math.pi * plane_strain_modulus / 16 / characteristic_length #Penny, characteristic_length is radius
    
    net_pressure_width = net_pressure / fracture_stiffness * gamma
    return net_pressure_width

def calc_FCD(fracture_width, half_length, formation_permeability, fracture_permeability, formation_height, fracture_height):
    #fracture_width:in   halfLength,fractureHeight,formationHeight:ft   formationPermability,fracturePermeability:mD
    #if using kf-w,: use width of 12 in
    _CONST = 1 / 12
    FCD = _CONST * (fracture_width * fracture_permeability) / (half_length * formation_permeability) * (fracture_height / formation_height)
    return FCD
    
def calc_fracture_skin(FCD, half_length, wellbore_radius):
    #halfLength,wellboreRadius:ft
    u = math.log(FCD)
    fracture_skin = (1.65 - 0.328 * u + 0.116 * u ** 2) / (1 + 0.18 * u + 0.064 * u ** 2 + 0.005 * u ** 3) - math.log(half_length / wellbore_radius)
    #fracture_skin = log(wellboreRadius / halfLength * (2 + (pi / 2) ** 2 / FCD))
    return fracture_skin

def calc_equivalent_wellbore_radius(wellbore_radius, skin):
    equivalent_wellbore_radius = wellbore_radius * math.exp(-skin)
    return equivalent_wellbore_radius

def calc_FOI(reservoir_radius, wellbore_radius, equivalent_wellbore_radius):
    FOI = math.log(reservoir_radius / wellbore_radius) / math.log(reservoir_radius / equivalent_wellbore_radius)
    return FOI

def calc_fracture_volume(fracture_width, half_length, fracture_height):
#fracture_width:in   half_length,fracture_height:ft    output:ft3
    _CONST = 1 / 12
    fracture_volume = 2 * _CONST * fracture_width * half_length * fracture_height
    return fracture_volume

def calc_fracture_mass(fracture_width, half_length, fracture_height, bulk_density):
    #fracture_width:in   half_length,fracture_height:ft    bulk_density:lbm/ft3      output:lbm proppant
    _CONST = 1 / 12
    fracture_mass = 2 * _CONST * fracture_width * half_length * fracture_height * bulk_density
    return fracture_mass

def calc_UFD(FCD, half_length, reservoir_length, proppant_volume, fracture_permeability, formation_permeability, formation_height):
        #half_length,reservoir_length,formation_height:ft      proppant_volume:ft3      fracture_permeability,formation_permeability:mD       output:[ft,in]
    _CONST = 12
    Ix = 2 * half_length / reservoir_length
    Nprop = FCD * Ix ** 2
    xfOpt = ((proppant_volume / 2 * fracture_permeability) / (FCD * formation_height * formation_permeability)) ** 0.5
    wOpt = _CONST * ((FCD * proppant_volume / 2 * formation_permeability) / (formation_height * fracture_permeability)) ** 0.5
    
    output = [Ix, Nprop, xfOpt, wOpt]
    UFD = output[:]
    return UFD

def calc_fracture_width(areal_proppant_concentration, solid_density, porosity):
    #areal_proppant_concentration:lbm/ft2   solid_density:ppg    porosity:fraction
    calc_fracture_width = areal_proppant_concentration * 1.604 / (solid_density * (1 - porosity))
    return calc_fracture_width

def calc_areal_proppant_concentration(fracture_width, solid_density, porosity):
    #fracture_width:in   solid_density:ppg    porosity:fraction
    areal_proppant_concentration = fracture_width * solid_density * (1 - porosity) / 1.604
    return areal_proppant_concentration

def calc_equivalent_time(dt, time_pumped):
    equivalent_time = dt * time_pumped / (dt + time_pumped)
    return equivalent_time

def calc_Gdef_time(time_total, time_pumped, alpha):
    dt = (time_total - time_pumped) / (time_pumped)
    if alpha == 0.5:
        g0 = math.pi / 2
        gdTd = (1 + dt) * math.asin((1 + dt) ** -0.5) + dt ** 0.5
    elif alpha == 1:
        g0 = 4 / 3
        gdTd = 4 / 3 * ((1 + dt) ** 1.5 - dt ** 1.5)
    else:
        g0 = 4 / 3 * ((1 + (0)) ** 1.5 - (0) ** 1.5 - 1)
        gdTd = 4 / 3 * ((1 + dt) ** 1.5 - dt ** 1.5 - 1)
    
    Gdef_time = 4 / math.pi * (gdTd - g0)
    return Gdef_time
    

    

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

def calc_API_tensile(tbg_od, tbg_id, yield_point):
    API_tensile = yield_point * math.pi * (tbg_od**2 - tbg_id**2) / 4
    return API_tensile

def calc_API_burst(tbg_od, tbg_id, yield_point, eccentricity=0.875):
    #yield_point,output:stress unit     tbg_od,tbg_id:length
    #Barlow's equation
    API_burst = yield_point * eccentricity * (tbg_od - tbg_id) / tbg_od
    return API_burst

def calc_Lame_stress(Ro, Ri, r, Po, Pi):
    #Ro,Ri,r:inch        Po,Pi,output:psi
    Lame_hoop_stress = -(Ri**2 * Ro**2 * (Po - Pi))/(Ro**2 - Ri**2)*(1/r**2) + (Pi*Ri**2 - Po*Ro**2)/(Ro**2 - Ri**2)
    Lame_radial_stress = (Ri**2 * Ro**2 * (Po - Pi))/(Ro**2 - Ri**2)*(1/r**2) + (Pi*Ri**2 - Po*Ro**2)/(Ro**2 - Ri**2)
    return Lame_hoop_stress, Lame_radial_stress

def calc_VM_stress_envelope(tbg_od, tbg_id, yield_point):
    #yield_point,output:psi      tbg_od,tbg_id:inch
    #outputs stress in the pipe
    hoop_stress = []
    axial_stress = []
    Ro = tbg_od / 2
    Ri = tbg_id / 2
    _C = (1 - Ro**2 / Ri**2) / (Ro**2 / Ri**2 + 1)
    hoop_stress.append(2 * yield_point / math.sqrt(3*(_C**2 - 2*_C + 1)))
    axial_stress.append(0.5*(_C + 1)*hoop_stress[-1] - 0.5*math.sqrt(max(0, 4*yield_point**2 - 3*(_C**2 - 2*_C + 1)*hoop_stress[-1]**2)))
    
    for _x in range(0, 100):
        hoop_stress.append(hoop_stress[-1] - hoop_stress[0]/50)
        axial_stress.append(0.5*(_C + 1)*hoop_stress[-1] - 0.5*math.sqrt(max(0, 4*yield_point**2 - 3*(_C**2 - 2*_C + 1)*hoop_stress[-1]**2)))
    for _x in range(0, 100):
        hoop_stress.append(hoop_stress[-1] + hoop_stress[0]/50)
        axial_stress.append(0.5*(_C + 1)*hoop_stress[-1] + 0.5*math.sqrt(max(0, 4*yield_point**2 - 3*(_C**2 - 2*_C + 1)*hoop_stress[-1]**2)))
    return hoop_stress, axial_stress

def calc_VM_envelope(tbg_od, tbg_id, yield_point, radius, eccentricity=0.875, temperature=75.0, DF_burst=1.0, DF_tension=1.0):
    #yield_point:psi      tbg_od,tbg_id,radius:inch        temperature:deg F      output:psi & lbf/1000
    #translates stresses to equivalent tension and well pressures
    hoop_stress = []
    axial_stress = []
    pressure = []
    tension = []
    Ro = tbg_od / 2
    Ri = tbg_id / 2
    _C = (1 - Ro**2 / radius**2) / (Ro**2 / radius**2 + 1)
    Ypt = yield_point * calc_pipe_temperature_derating(temperature)
    cross_section_area = math.pi * (Ro**2 - Ri**2)
    uniaxial_burst = calc_API_burst(tbg_od, tbg_id, Ypt, eccentricity)
    hoop_stress.append(2 * Ypt / math.sqrt(3*(_C**2 - 2*_C + 1)))
    #max() addresses floating point error for working near roots
    axial_stress.append(0.5*(_C + 1)*hoop_stress[-1] - 0.5 * math.sqrt(max(0, 4*Ypt**2 - 3*(_C**2 - 2*_C + 1)*hoop_stress[-1]**2)))
    
    for _x in range(0, 100):
        hoop_stress.append(hoop_stress[-1] - hoop_stress[0]/50)
        axial_stress.append(0.5*(_C + 1)*hoop_stress[-1] - 0.5 * math.sqrt(max(0, 4*Ypt**2 - 3*(_C**2 - 2*_C + 1)*hoop_stress[-1]**2)))
    for _x in range(0,100):
        hoop_stress.append(hoop_stress[-1] + hoop_stress[0]/50)
        axial_stress.append(0.5*(_C + 1)*hoop_stress[-1] + 0.5 * math.sqrt(max(0, 4*Ypt**2 - 3*(_C**2 - 2*_C + 1)*hoop_stress[-1]**2)))
    for _x in range(len(hoop_stress)):
        pressure.append(uniaxial_burst * hoop_stress[_x] / (Ypt / math.sqrt(_C**2 - _C + 1)) / DF_burst)
        tension.append(axial_stress[_x] * cross_section_area / DF_tension)
    return pressure, tension

def calc_API_collapse(tbg_od, tbg_id, yield_point, Youngs_modulus, Poissons_ratio, temperature=75.0, axial_stress=0.0):
    #Yp,axial_stress,output:psi      outer_diameter,inner_diameter:inch     temperature:deg F
    Dt = tbg_od / ((tbg_od - tbg_id) / 2) 
    Ypt = yield_point * calc_pipe_temperature_derating(temperature)
    Ypa = Ypt * (math.sqrt(1 - 0.75*(axial_stress/Ypt)**2) - 0.5*axial_stress/yield_point)
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
        API_collapse = (2*Youngs_modulus)/(1-Poissons_ratio**2) * 1/(Dt * (Dt - 1)**2)
    return API_collapse    

#Buckling and Tubemove Equations
def calc_dL_HookesLaw(length, force, Youngs_modulus, area):   #units of L
    dL_HookesLaw = -length * force / Youngs_modulus / area
    return dL_HookesLaw

def calc_force_HookesLaw(length, dL, Youngs_modulus, area):   #lbf
    force_HookesLaw = dL * Youngs_modulus * area / length
    return force_HookesLaw

def calc_buoyed_weight(pipeWt, tbg_fluid_density, IDArea, csg_fluid_density, ODArea): #lb/ft
    buoyed_weight = pipeWt + (tbg_fluid_density * IDArea - csg_fluid_density * ODArea) * 12 / 231
    return buoyed_weight

def calc_force_pasley(buoy_wt, angle, Youngs_modulus, moment_inertia, rc):        #lbf
    criticalAngle = math.asin((1.94 / 2) ** 2 * rc * ((buoy_wt / 12) / (Youngs_modulus * moment_inertia)) ** (1 / 3))
    if (angle) >= criticalAngle:
        force_pasley = math.sqrt(4 * (buoy_wt / 12) * math.sin(angle) * Youngs_modulus * moment_inertia / rc)
    else:
        force_pasley = 1.94 * (Youngs_modulus * moment_inertia * (buoy_wt / 12) ** 2) ** (1 / 3)
    return force_pasley
    
def calc_normal_buckling_force(rc, force, Youngs_modulus, moment_inertia):        #lbf/ft
    normal_buckling_force = rc * force ** 2 / (4 * Youngs_modulus * moment_inertia) * 12
    return normal_buckling_force

def calc_buckling_dogleg(force, pasley_force, Youngs_modulus, moment_inertia, buckling_state: str):
    if buckling_state == "lateral":
        buckling_dogleg = 1.1227 / math.sqrt(2 * Youngs_modulus * moment_inertia) * force ** 0.04 * (force - pasley_force) ** 0.46
    elif buckling_state == "helical":
        buckling_dogleg = math.sqrt(force / (2 * Youngs_modulus * moment_inertia))
    else:
        buckling_dogleg = 0
    return buckling_dogleg
    
def calc_buckling_strain(force, pasley_force, Youngs_modulus, moment_inertia, rc, buckling_state: str):  #unit strain
    if buckling_state == "lateral":
        buckling_strain = -0.7285 * rc ** 2 / (4 * Youngs_modulus * moment_inertia) * force ** 0.08 * (force - pasley_force) ** 0.92
    elif buckling_state == "helical":
        buckling_strain = -rc ** 2 / (4 * Youngs_modulus * moment_inertia) * force
    else:
        buckling_strain = 0
    return buckling_strain

def calc_dl_buckling(force, force1, axialBuoyWt, pasley_force, Youngs_modulus, moment_inertia, rc, buckling_state: str):  #from Mitchell, ft
    if buckling_state == "lateral":
        dl_buckling = -(rc ** 2 / (4 * Youngs_modulus * moment_inertia * axialBuoyWt) * (force - pasley_force) * (0.3771 * force - 0.3668 * pasley_force))
    elif buckling_state == "helical":
        dl_buckling = -(rc ** 2 / (8 * Youngs_modulus * moment_inertia * axialBuoyWt) * (force ** 2 - force1 ** 2))
    else:
        dl_buckling = 0
    return dl_buckling
    
def calc_dL_balloon(length, Poissons_ratio, Youngs_modulus, dPTbg, dDensityTbg, dPAnn, dDensityAnn, R, d):    #in
    #D is pressure drop of fluid through the pipe, r = od/id of tubing
    #dL_balloon = -1 * Length ** 2 * Poissons_ratio / Youngs_modulus * (dDensityTbg - r ** 2 * dDensityAnn - (1 + 2 * Poissons_ratio) * D / (2 * Poissons_ratio)) / (r ** 2 - 1) - 2 * Length * Poissons_ratio / Youngs_modulus * (dPTbg - r ** 2 * dPAnn) / (r ** 2 - 1)
    dL_balloon = -2 * length * Poissons_ratio / Youngs_modulus * (dPTbg - R ** 2 * dPAnn) / (R ** 2 - 1)
    return dL_balloon

def calc_dL_temperature(alpha, dTemp, length):   #ft
    dL_temperature = alpha * dTemp * length
    return dL_temperature

def calc_dF_piston_packer(Ap, Ai, Ao, dPi, dPo):   #lbf, negative is pushing down against the packer (tensile force)
    dF_piston_packer = (Ap - Ai) * dPi - (Ap - Ao) * dPo
    return dF_piston_packer

def calc_dF_piston_xo(Ai1, Ai2, Ao1, Ao2, dPi, dPo):    #lbf
    dF_piston_xo = (Ai1 - Ai2) * dPi - (Ao1 - Ao2) * dPo
    return dF_piston_xo

def calc_dF_bend(F2, ff, angle):  #angle in radians
    #dF_bend = F2 * (exp(ff * angle) - 1) #calculates dF across a bend
    dF_bend = 2 * F2 * math.sin(angle / 2)    #alternative for normal force calculation
    #dF_bend = 0
    return dF_bend

def calc_ws_slackoff_force(force, axial_buoy_wt, normal_buoy_wt, normal_buckling_wt, pulley_wt, friction, length):
    ws_slackoff_force = force - (axial_buoy_wt - (normal_buoy_wt + normal_buckling_wt + pulley_wt) * friction) * length
    return ws_slackoff_force

def calc_ws_pickup_force(force, axial_buoy_wt, normal_buoy_wt, normal_buckling_wt, pulley_wt, friction, length):
    ws_pickup_force = force - (axial_buoy_wt + (normal_buoy_wt + normal_buckling_wt + pulley_wt) * friction) * length
    return ws_pickup_force





#*************************************************************************
#PRODUCTION EQUATIONS
#*************************************************************************
def calc_shape_radial(re, rw, S, state: str):
    #re and rw same units
    if state == "Steady State":
        shape_radial = math.log(re / rw) + S          #steady state
    else:
        shape_radial = math.log(0.472 * re / rw) + S  #pseudosteady state
    return shape_radial

def calc_shape_horizontal(kh, kv, re, rw, h, L, S, state: str):
    #kh & kv same units, re & rw same units, h & L same units
    if state == "Steady State":
        St = 1          #steady state
    else:
        St = 0.472     #pseudosteady state
    Iani = math.sqrt(kh / kv)
    a = L / 2 * (0.5 + (0.25 + (re / (0.5 * L)) ** 4) ** 0.5) ** 0.5
    shape_horizontal = math.log((a + math.sqrt(a ** 2 - (L / 2) ** 2)) / (L / 2)) + Iani * h / L * (math.log(St * Iani * h / (rw * (Iani + 1))) + S)
    return shape_horizontal

def calc_shape_dietz(c, a, rw, S):
    #A and rw in same units
    gamma = 1.78    #constant
    shape_dietz = 0.5 * math.log(4 * a / (gamma * c * rw ** 2)) + S
    return shape_dietz

def calc_Darcy_IPR(k, h, pr, pb, pwf, B, m, shape):   #calculate oil well IPR using Darcy & Vogel correlation
    #k:mD    h:ft    pr,pb,pwf:psi   m:cP    B:bbl/stb   #output stb/day
    _CONST = 0.001127           #constant to convert units to stb/day

    if pwf < pb:  #Vogel solution
        qb = (2 * math.pi * _CONST * k * h * (pr - pb)) / (B * m * shape)
        qv = qb * (pb / 1.8) / (pr - pb)
        q = qb + qv * (1 - 0.2 * (pwf / pb) - 0.8 * (pwf / pb) ** 2)
    else:            #Darcy solution
        q = (2 * math.pi * _CONST * k * h * (pr - pwf)) / (B * m * shape)
    Darcy_IPR = q
    return Darcy_IPR

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
    Wiggins_IPR = q
    return Wiggins_IPR

def calc_gas_IPR(k, h, pr, pwf, m, T, Z, shape):  #calculate gas well IPR
    #input D, ft, psi, cP, F    #output Mscf/day
    _CONST = 0.111807        #constant to convert units

    q = (2 * math.pi * _CONST * k * h * (pr ** 2 - pwf ** 2)) / (T * Z * m * shape)
    gas_IPR = q
    return gas_IPR

def calc_perforating_skin(well_radius, phase_angle, perf_length, perf_height, perf_EHD, horizontal_perm, vertical_perm): 
    #calculate perforating skin
    if phase_angle == 0:
        a_theta = 0.25
        a1 = -2.091
        a2 = 0.0453
        b1 = 5.1313
        b2 = 1.8672
        c1 = 0.16
        c2 = 2.675
    elif phase_angle == 180:
        a_theta = 0.5
        a1 = -2.025
        a2 = 0.0943
        b1 = 3.0373
        b2 = 1.8115
        c1 = 0.026
        c2 = 4.532
    elif phase_angle == 120:
        a_theta = 0.648
        a1 = -2.018
        a2 = 0.0634
        b1 = 1.6136
        b2 = 1.777
        c1 = 0.0066
        c2 = 5.32
    elif phase_angle == 90:
        a_theta = 0.726
        a1 = -1.905
        a2 = 0.1038
        b1 = 1.5674
        b2 = 1.6935
        c1 = 0.0019
        c2 = 6.155
    elif phase_angle == 60:
        a_theta = 0.813
        a1 = -1.898
        a2 = 0.1023
        b1 = 1.3654
        b2 = 1.649
        c1 = 0.0003
        c2 = 7.509
    elif phase_angle == 45:
        a_theta = 0.86
        a1 = -1.788
        a2 = 0.2398
        b1 = 1.1915
        b2 = 1.6392
        c1 = 0.000046
        c2 = 8.791

    #SH, plane flow effect
    if phase_angle == 0:
        rwPrime = perf_length / 4
    else: 
        rwPrime = a_theta * (well_radius + perf_length)
    
    SH = math.log(well_radius / rwPrime)

    #SV, vertical converging effect
    hD = perf_height / perf_length * (horizontal_perm / vertical_perm) ** 0.5
    rD = perf_EHD / (4 * perf_height) * (1 + (vertical_perm / horizontal_perm) ** 0.5)
    a = a1 * math.log10(rD) + a2
    B = b1 * rD + b2
    SV = 10 ** a * hD ** (B - 1) * rD ** B

    #SWB, wellbore effect
    rwD = well_radius / (perf_length + well_radius)
    SWB = c1 * math.exp(c2 * rwD)
    perforating_skin = SH + SV + SWB
    return perforating_skin

def calc_VcVs(rate, formation_volume_factor, formation_perm, formation_height, fluid_density, fluid_viscosity, skin, csg_od, csg_id, perf_interval_length, perf_diameter, proppant_perm, beta_factor, min_ann_gap):
    _CONST = 141.2
    dP_skin_mechanical = rate * _CONST * formation_volume_factor * fluid_viscosity * skin / (formation_perm * formation_height)
   
    perf_length = (csg_od - csg_id) / 2
    perf_area = calc_area(perf_diameter, 0)
    B = 10646.7 * fluid_viscosity * perf_length / proppant_perm
    a = 1.57501 * 10 ** -9 * fluid_density * perf_length * beta_factor
    effectiveSPF = (formation_volume_factor * rate) / ((-B + (B ** 2 + 4 * a * dP_skin_mechanical) ** 0.5) / (2 * a) * perf_area) / perf_interval_length
    
    csg_velocity = ((-B + (B ** 2 + 4 * a * dP_skin_mechanical) ** 0.5) / (2 * a)) * 144 * 5.6146 / 24 / 60 / 60
    screen_velocity = csg_velocity / (1 + 5.1285 * min_ann_gap + 9.0898 * min_ann_gap ** 2 - 0.7236 * min_ann_gap ** 3)
    
    return dP_skin_mechanical, effectiveSPF, csg_velocity, screen_velocity






#****************************************************************
#PVT CALCULATIONS
#****************************************************************
def calc_bubble_point(method: str, Rs, oil_density, gas_gravity, temperature, pressure=0, Tsp=0, Psp=0):
#Rs:scf/stb    oil_density:API  gas_gravity:SG  temperature:F    pressure:psi
#Recommended methods are Glaso, Standing, Lasater, Vasquez & Beggs, Petrosky, Al-Marhoun, Velarde

    if method == "Standing":
        pb = 18.2 * ((Rs / gas_gravity) ** 0.83 * 10 ** (0.00091 * temperature - 0.0125 * oil_density) - 1.4)
    
    elif method == "Elam":
        pb = Rs ** 0.702 / gas_gravity ** 0.514 * math.exp(0.00348 * temperature - 0.0282 * oil_density + 3.58)
    
    elif method == "Lasater":
        Mo = 6084 / (oil_density - 5.9)
        xg = (1 + ucon(oil_density, "API", "sg (water)") / (0.000007521 * Rs * Mo)) ** -1
        pf = math.exp((xg - 0.15649) / 0.33705) - 0.59162
        pb = pf * (temperature + 459.67) / gas_gravity
    
    elif method == "Vazquez & Beggs":
        if oil_density <= 30:
            a = 27.64
            B = -11.172
            c = 0.9143
        else:
            a = 56.06
            B = -10.393
            c = 0.8425
        
        gas_gravityc = gas_gravity * (1 + 0.00005912 * oil_density * Tsp * math.log10(Psp / 114.7))
        pb = (a * (Rs / gas_gravityc) * 10 ** (B * oil_density / (temperature + 459.67))) ** c
    
    elif method == "Glaso - non-volatile":
        x = (Rs / gas_gravity) ** 0.816 * (temperature ** 0.172 / oil_density ** 0.989)    #non-volatile oils
        pb = 10 ** (1.7669 + 1.7447 * math.log10(x) - 0.30218 * (math.log10(x)) ** 2)
    
    elif method == "Glaso - volatile":
        x = (Rs / gas_gravity) ** 0.816 * (temperature ** 0.13 / oil_density ** 0.989)    #volatile oil
        pb = 10 ** (1.7669 + 1.7447 * math.log10(x) - 0.30218 * (math.log10(x)) ** 2)
    
    elif method == "Labedi":
        pb = 6.0001 / gas_gravity * (Rs ** 0.6714 * (temperature / oil_density) ** 0.7097 * Tsp ** 0.08929) / (10 ** (Rs * 7.995 * 10 ** -5))
    
    elif method == "Owolabi":
        #pb = 55 + 0.8643 * ((Rs / gas_gravity) ** 1.255 * (Temperature ** 0.172 / oil_density ** 0.178))   #old solution
        pb = -987.56359 + 179.58816 * ((Rs / gas_gravity) ** 0.48088266 * (temperature ** 0.09353815 / oil_density ** 0.16648326))
    
    elif method == "Al-Marhoun 1985":
        oil_density = ucon(oil_density, "API", "sg (water)")
        x = Rs ** 0.722569 * oil_density ** 3.04659 / gas_gravity ** 1.879109 * (temperature + 459.67) ** 1.302347
        pb = -64.13891 + 7.02362 * 10 ** -3 * x - 2.278475 * 10 ** -9 * x ** 2
    
    elif method == "Obomanu & Okpobiri":
        pb = ((Rs * temperature ** 0.497 * 10 ** 0.811) / (1.01136371 * gas_gravity ** 2.15 * oil_density ** 1.27)) ** 1.0787
    
    elif method == "Al-Marhoun 1988":
        oil_density = ucon(oil_density, "API", "sg (water)")
        pb = 5.38088 * 10 ** -3 * Rs ** 0.715082 * oil_density ** 3.1437 * (temperature + 459.67) ** 1.32657 / gas_gravity ** 1.87784
    
    elif method == "Asgarpour-Viking":
        temperature = ucon(temperature, "F", "C")
        a1 = 70.9815
        a2 = -0.0101
        a3 = -0.2514
        a4 = 0.4593
        a5 = 0.5093
        pb = a1 * gas_gravity ** a2 * oil_density ** a3 * temperature ** a4 * (Rs / 5.614583) ** a5
    
    elif method == "Asgarpour-Nisku":
        temperature = ucon(temperature, "F", "C")
        a1 = 83.7883
        a2 = -0.4114
        a3 = -0.2697
        a4 = 0.281
        a5 = 0.6259
        pb = a1 * gas_gravity ** a2 * oil_density ** a3 * temperature ** a4 * (Rs / 5.614583) ** a5
    
    elif method == "Asgarpour-Leduc":
        temperature = ucon(temperature, "F", "C")
        a1 = 193.77
        a2 = 0.0928
        a3 = -1.1369
        a4 = 0.7899
        a5 = 0.6519
        pb = a1 * gas_gravity ** a2 * oil_density ** a3 * temperature ** a4 * (Rs / 5.614583) ** a5
    
    elif method == "Al-Najjar, Al-Soof & Al-Khalisy":
        if oil_density <= 30:
            a = 7.92
            B = 1.025
            c = -24.244
        else:
            a = 30.91
            B = 0.816
            c = -19.748
        
        pb = a * (Rs / gas_gravity) ** B * math.exp(c * oil_density / (temperature + 459.67))
    
    elif method == "Dokla & Osma":
        oil_density = ucon(oil_density, "API", "sg (water)")
        pb = 8363.86 * Rs ** 0.724047 * oil_density ** 0.107991 / (gas_gravity ** 1.01049 * (temperature + 459.67) ** 0.952584)
    
    elif method == "Petrosky":
        pb = 112.727 * (Rs ** 0.577421 / (gas_gravity ** 0.8439 * 10 ** (7.916 * 10 ** -4 * oil_density ** 1.541 - 4.561 * 10 ** -5 * temperature ** 1.3911)) - 12.34)
    
    elif method == "Kartoatmodjo & Schmidt":
        if oil_density <= 30:
            a = 0.05958
            B = 0.7972
            c = 13.1405
            d = 0.998602
        else:
            a = 0.0315
            B = 0.7587
            c = 11.289
            d = 0.914328
        
        gas_gravityc = gas_gravity * (1 + 0.1595 * oil_density ** 0.4078 / Tsp ** 0.2466 * math.log10(Psp / 114.7))
        pb = (Rs / (a * gas_gravityc ** B * 10 ** (c * oil_density / (temperature + 459.67)))) ** d
    
    elif method == "Farshad":
        pb = 64.14 * (Rs ** 0.6343 / (gas_gravity ** 1.15036 * 10 ** (7.97 * 10 ** -3 * oil_density - 3.35 * 10 ** -4 * temperature)) - 7.2818)
    
    elif method == "Macary":
        pb = 204.257 * math.exp(7.7 * 10 ** -4 * temperature - 9.7 * 10 ** -3 * oil_density - 0.4003 * gas_gravity) * (Rs ** 0.51 - 4.7927)
    
    elif method == "Omar & Todd":
        pb = 0
    
    elif method == "Hasan":
        pb = 18.3 * ((Rs / gas_gravity) ** 0.83 * 10 ** (0.00091 * temperature - 0.0125 * oil_density) + 2.2)
    
    elif method == "De Ghetto":    #questionable if gas_gravity or gas_gravityc
        if oil_density <= 10:
            pb = 10.7025 * (Rs / gas_gravity) ** 0.8986 * 10 ** (0.00156 * temperature - 0.00169 * oil_density)
        elif (oil_density > 10) and (oil_density <= 22.3):
            pb = (56.434 * Rs / (gas_gravity * 10 ** (10.9267 * oil_density / (temperature + 459.67)))) ** 0.8294
        elif (oil_density > 22.3) and (oil_density <= 31.1):
            pb = (Rs / (0.10084 * gas_gravity ** 0.2556 * 10 ** (7.4576 * oil_density / (temperature + 459.67)))) ** 1.0134
        else:
            pb = (Rs / (0.1347 * gas_gravity ** 0.3873 * 10 ** (12.753 * oil_density / (temperature + 459.67)))) ** 0.8536
        
    
    elif method == "Agip":    #questionable if gas_gravity or gas_gravityc
        pb = (37.966 * Rs / (gas_gravity * 10 ** (9.441 * oil_density / (temperature + 459.67)))) ** 0.8669
    
    elif method == "Almedhaideb":
        Bo = 1.122018 + (1.41 * 10 ** -6 * Rs * temperature / oil_density ** 2)
        pb = -620.592 + 6.23087 * Rs * oil_density / (gas_gravity * Bo ** 1.38559) + 2.89868 * temperature
    
    elif method == "Elsharkawy":
        if oil_density <= 30:
            pb = (Rs / (gas_gravity ** 0.04439 * oil_density ** 1.1394 * 10 ** (8.392 * 10 ** -4 * temperature - 2.188))) ** 1.0551194
        else:
            pb = (Rs / (gas_gravity * 10 ** (0.4636 * oil_density / temperature - 1.2179))) ** 0.847271
        
    
    elif method == "Khairy":
        pb = 49.3647 * Rs ** 0.5774 * temperature ** 0.6641 / (gas_gravity ** 1.4676 * oil_density ** 1.0305)
    
    elif method == "Al-Shammasi":
        oil_density = ucon(oil_density, "API", "sg (water)")
        pb = oil_density ** 5.527215 * (gas_gravity * Rs * (temperature + 459.67)) ** 0.783716 / math.exp(1.841408 * oil_density * gas_gravity)
    
    elif method == "Levitan & Murtha":
        oil_density = ucon(oil_density, "API", "sg (water)")
        pb = 14.7 * (Rs / gas_gravity) ** 0.85 * oil_density ** 5 * ((temperature + 459.67) / 519.67) ** 1.5
    
    elif method == "Velarde":
        #A = 9.73 * 10 ** -7 * gas_gravity ** 1.1672608 * oil_density ** 0.92987 * Temperature ** 0.247235 * pb ** 1.056052
        #B = 0.022339 * gas_gravity ** -1.00475 * oil_density ** 0.337711 * Temperature ** 0.132795 * pb ** 0.302065
        #C = 0.725167 * gas_gravity ** -1.48548 * oil_density ** -0.164741 * Temperature ** -0.09133 * pb ** 0.047094
        #pr = Pressure / pb
        #Rs = Rsb * (A * pr ** B + (1 - A) * pr ** C)
        #use above to solve for Rs in future
        x = 0.013098 * temperature ** 0.282372 - 8.2 * 10 ** -6 * oil_density ** 2.176124
        pb = 1091.47 * (Rs ** 0.081465 * 10 ** x / gas_gravity ** 0.161488 - 0.740152) ** 5.354891
    
    elif method == "Dindoruk & Christman":
        a = (1.42828 * 10 ** -10 * temperature ** 2.844591797 - 6.74896 * 10 ** -4 * oil_density ** 1.225226436) / (0.033383304 + 2 * gas_gravity ** 0.084226069 * Rs ** -0.272945957) ** 2
        pb = 1.86997927 * (Rs ** 1.221486524 * 10 ** a / gas_gravity ** 1.370508349 + 0.011688308)
    
    else:
        pb = "Select Method"

    bubble_point = pb
    return bubble_point

def calc_oil_formation_factor(method: str, Rs, oil_density, gas_gravity, temperature, pressure=0, pb=0, Psp=0, Tsp=0):
#Rs: scf/stb    oil_density: API  gas_gravity: SG  Temperature: F    Pressure: psi
#Recommended methods are Glaso, Standing, Lasater, Vasquez & Beggs, Petrosky, Al-Marhoun, Velarde

    if method == "Standing":
        oil_density = ucon(oil_density, "API", "sg (water)")
        B = 0.972 + 1.47 * 10 ** -4 * (Rs * (gas_gravity / oil_density) ** 0.5 + 1.25 * temperature) ** 1.175
    
    elif method == "Elam":
        oil_density = ucon(oil_density, "API", "sg (water)")
        B = math.exp(-0.0355 + 3.55 * 10 ** -4 * Rs(gas_gravity / oil_density) ** 0.5 + 7.1 * 10 ** -4 * temperature)
    
    elif method == "Vazquez & Beggs":
        if oil_density <= 30:
            a1 = 4.677 * 10 ** -4
            a2 = 1.751 * 10 ** -5
            a3 = -1.8106 * 10 ** -8
        else:
            a1 = 4.67 * 10 ** -4
            a2 = 1.1 * 10 ** -5
            a3 = 1.337 * 10 ** -9
        
        B = 1 + a1 * Rs + (temperature - 60) * (oil_density / gas_gravity) * (a2 + a3 * Rs)
    
    elif method == "Glaso":
        oil_density = ucon(oil_density, "API", "sg (water)")
        x = Rs * (gas_gravity / oil_density) ** 0.526 + 0.968 * temperature
        B = 10 ** (-6.58511 + 2.91329 * math.log10(x) - 0.27683 * (math.log10(x)) ** 2) + 1
    
    elif method == "Labedi":
        Bob = 0.9976 + 5.273 * 10 ** -4 * Rs + 2.6636 * 10 ** -8 * (temperature - 60) * (oil_density * Psp) + 1.6982 * 10 ** -5 * oil_density * (temperature - 60)
        if Bob <= 1.758:
            x = (3.61 * Rs ** 0.4625 * Bob ** 13.398 * Psp ** 0.0775) / 10 ** (3.231 * Bob)
            B = Bob - x * (1 - pressure / pb)
        else:
            x = 1.6339 - 9.152 * 10 ** -4 * Rs + 1.584 * 10 ** -7 * Rs ** 2
            B = Bob - (1 - pressure / pb) ** x
        
    elif method == "Owolabi":
        oil_density = ucon(oil_density, "API", "sg (water)")
        #B = 0.9871 + 4.0689 * 10 ** -4 * (Rs * (gas_gravity / oil_density) ** 0.526 + 1.25 * Temperature) #old solution
        B = 0.9957 + 3.7921 * 10 ** -4 * (Rs * (gas_gravity / oil_density) ** 0.526 + 1.25 * temperature)
    
    elif method == "Al-Marhoun 1985": #recommended
        oil_density = ucon(oil_density, "API", "sg (water)")
        x = Rs ** 0.501538 * gas_gravity ** -0.145526 * oil_density ** -5.220726
        B = 0.574095 + 7.723532 * 10 ** -4 * (temperature + 459.67) + 2.454005 * 10 ** -3 * x + 3.727676 * 10 ** -5 * x ** 2
    
    elif method == "Obomanu & Okpobiri":
        oil_density = ucon(oil_density, "API", "sg (water)")
        if oil_density < 0.876:
            B = 0.3321 + 1.404154 * 10 ** -4 * Rs + 4.1588128 * 10 ** -4 * Rs * gas_gravity / oil_density + 1.15861 * 10 ** -5 * (temperature + 459.67)
        else:
            B = 1.0232 + 2.725 * 10 ** -5 * (Rs * (gas_gravity / oil_density + temperature)) ** 0.79
        
    elif method == "Al-Marhoun 1988": #recommended
        oil_density = ucon(oil_density, "API", "sg (water)")
        x = Rs ** 0.74239 * gas_gravity ** 0.323294 * oil_density ** -1.20204
        B = 0.497069 + 8.62963 * 10 ** -4 * (temperature + 459.67) + 1.82594 * 10 ** -3 * x + 3.18099 * 10 ** -6 * x ** 2
    
    elif method == "Asgarpour-Viking":
        temperature = ucon(temperature, "F", "C")
        a1 = 0.1203
        a2 = 0.0645
        a3 = 0.2452
        a4 = 0.1118
        a5 = 0.2321
        a6 = 1.342
        a7 = -0.5811
        a8 = 0.183
        a9 = 0.1656
        x = a1 * gas_gravity ** a2 * oil_density ** a3 * temperature ** a4 * (Rs / 5.614583) ** a5
        B = a6 + a7 * x + a8 * x ** 2 + a9 * x ** 3
    
    elif method == "Asgarpour-Nisku":
        temperature = ucon(temperature, "F", "C")
        a1 = 0.251
        a2 = 0.0724
        a3 = 0.00275
        a4 = 0.1538
        a5 = 0.2235
        a6 = -2.3211
        a7 = 7.039
        a8 = -5.006
        a9 = 1.33
        x = a1 * gas_gravity ** a2 * oil_density ** a3 * temperature ** a4 * (Rs / 5.614583) ** a5
        B = a6 + a7 * x + a8 * x ** 2 + a9 * x ** 3
    
    elif method == "Asgarpour-Leduc":
        temperature = ucon(temperature, "F", "C")
        a1 = 0.1941
        a2 = 0.0136
        a3 = 0.0912
        a4 = 0.159
        a5 = 0.211
        a6 = 0.8603
        a7 = 0.7341
        a8 = -0.9378
        a9 = 0.4686
        x = a1 * gas_gravity ** a2 * oil_density ** a3 * temperature ** a4 * (Rs / 5.614583) ** a5
        B = a6 + a7 * x + a8 * x ** 2 + a9 * x ** 3
    
    elif method == "Al-Najjar":
        oil_density = ucon(oil_density, "API", "sg (water)")
        B = 0.96325 + 4.9 * 10 ** -4 * (Rs * (gas_gravity / oil_density) ** 0.5 + 1.25 * temperature)
    
    elif method == "Ahmed":
        B = -0.12869353 + Rs ** 0.023484894 * (oil_density ** 0.015966573 / gas_gravity ** 0.021946351) - 4.5243973 * 10 ** -4 * temperature + \
            3.9063637 * 10 ** -6 * temperature ** 2 - 5.5542509 / temperature - 5.760322 * 10 ** -6 * pressure - \
            3.9528992 * 10 ** -9 * pressure ** 2 + 16.289473 / pressure + 3.8718887 * 10 ** -4 * Rs + \
            7.0703685 * 10 ** -8 * Rs ** 2 - 1.4358395 / Rs
    
    elif method == "Abdul-Majeed & Salman":
        oil_density = ucon(oil_density, "API", "sg (water)")
        x = Rs ** 1.2 / (gas_gravity ** 0.147 * oil_density ** 5.222)
        B = 0.9657876 + 4.8141 * 10 ** -5 * x - 6.8987 * 10 ** -10 * x ** 2 + 7.73 * 10 ** -4 * temperature
    
    elif method == "Dokla & Osman":
        oil_density = ucon(oil_density, "API", "sg (water)")
        x = Rs ** 0.773572 * (gas_gravity ** 0.40402 / oil_density ** 0.882605)
        B = 4.31935 * 10 ** -2 + 1.56667 * 10 ** -3 * (temperature + 459.67) + 1.39775 * 10 ** -3 * x + 3.80525 * 10 ** -6 * x ** 2
    
    elif method == "Petrosky":
        oil_density = ucon(oil_density, "API", "sg (water)")
        B = 1.0113 + 7.2046 * 10 ** -5 * (Rs ** 0.3738 * gas_gravity ** 0.2914 / oil_density ** 0.6265 + 0.24626 * temperature ** 0.5371) ** 3.0936
    
    elif method == "Kartoatmodjo & Schmidt":   #recommended
        gas_gravityc = gas_gravity * (1 + 0.1595 * oil_density ** 0.4078 / Tsp ** 0.2466 * math.log10(Psp / 114.7))
        oil_density = ucon(oil_density, "API", "sg (water)")
        B = 0.98496 + 1 * 10 ** -4 * (Rs ** 0.755 * gas_gravityc ** 0.25 / oil_density ** 1.5 + 0.45 * temperature) ** 1.5
    
    elif method == "Al-Marhoun 1992":          #recommended
        oil_density = ucon(oil_density, "API", "sg (water)")
        B = 1 + 1.77342 * 10 ** -4 * Rs + 2.20163 * 10 ** -4 * Rs * gas_gravity / oil_density + 4.29258 * 10 ** -6 * Rs * (temperature - 60) * (1 - oil_density) + 5.28707 * 10 ** -4 * (temperature - 60)
    
    elif method == "Farshad":                  #recommended
        oil_density = ucon(oil_density, "API", "sg (water)")
        x = Rs ** 0.5956 * (gas_gravity ** 0.2369 / oil_density ** 1.3282) + 0.0976 * temperature
        B = 1 + 10 ** (-2.6541 + 0.557 * math.log10(x) + 0.3331 * math.log10(x) ** 2)
    
    elif method == "Macary":
        oil_density = ucon(oil_density, "API", "sg (water)")
        B = (1.0031 + 0.0008 * temperature) * math.exp(0.0004 * Rs + 0.0006 * oil_density / gas_gravity)
    
    elif method == "Omar & Todd":
        x = 1.1663 + 7.62 * 10 ** -4 * oil_density / gas_gravity - 0.0339 * gas_gravity
        oil_density = ucon(oil_density, "API", "sg (water)")
        B = 0.972 + 1.47 * 10 ** -4 * (Rs * (gas_gravity / oil_density) ** 0.5 + 1.25 * temperature) ** x
    
    elif method == "Almedhaideb":
        oil_density = ucon(oil_density, "API", "sg (water)")
        B = 1.122018 + 1.41 * 10 ** -6 * Rs * temperature / oil_density ** 2
    
    elif method == "Elsharkawy":
        oil_density = ucon(oil_density, "API", "sg (water)")
        B = 1 + 4.0428 * 10 ** -4 * Rs + 6.3802 * 10 ** -4 * (temperature - 60) + 7.8 * 10 ** -7 * Rs * (temperature - 60) * gas_gravity / oil_density
    
    elif method == "Khairy":
        B = 0.773413 + 7.05341 * 10 ** -4 * Rs + 0.18669 * gas_gravity - 9.2589 * 10 ** -4 * oil_density + 4.41 * 10 ** -4 * temperature
    
    elif method == "Al-Shammasi":              #recommended
        oil_density = ucon(oil_density, "API", "sg (water)")
        B = 1 + 5.53 * 10 ** -7 * Rs * (temperature - 60) + 1.81 * 10 ** -4 * Rs / oil_density + 4.49 * 10 ** -4 * (temperature - 60) / oil_density + 2.06 * 10 ** -4 * Rs * gas_gravity / oil_density
    
    elif method == "Levitan & Murtha":
        oil_density = ucon(oil_density, "API", "sg (water)")
        B = 1 + 0.0005 * Rs * (gas_gravity / oil_density) ** 0.25 + 0.2 * (temperature - 60) / (519.67 * gas_gravity * oil_density)
    
    #Velarde - requires iterations
    
    elif method == "Dindoruk & Christman":
        oil_density = ucon(oil_density, "API", "sg (water)")
        x = (Rs ** 2.510755 / (gas_gravity ** 4.852538 * oil_density ** 11.835) + 1.365428 * 10 ** 5 * (temperature - 60) ** 2.25288 + 10.0719 * Rs) ** 0.4450849 /  (5.352624 + 2 * Rs ** -0.6309052 * gas_gravity ** -0.9000749 * (temperature - 60)) ** 2
        gco = ucon(oil_density, "sg (water)", "API")
        B = 0.9871766 + 7.865146 * 10 ** -4 * x + 2.689173 * 10 ** -6 * x ** 2 + 1.100001 * 10 ** -5 * (temperature - 60) * gco / gas_gravity
    
    else:
        B = "Select method"
    
    oil_formation_factor = B
    return oil_formation_factor

def calc_oil_viscosity(method: str, oil_density, temperature, pb=0, Rs=0):     #Dead oil
#oil_density: API  gas_gravity: SG  Temperature: F    Pressure: psi
#Recommended are Beal, Beggs, Petrosky, Egbogah, Bergman-Sutton


    if method == "Andrade":  #need to research for A & B
        B = 0
        
    elif method == "Beal":
        x = math.exp(2.302585 * (0.43 + 8.33 / oil_density))
        m = (0.32 + 1.8 * 10 ** 7 / oil_density ** 4.53) * (360 / (temperature + 200)) ** x
    
    elif method == "Beggs & Robinson":
        x = 10 ** (3.0324 - 0.02023 * oil_density) * temperature ** -1.163
        m = 10 ** x - 1
    
    elif method == "Glaso":
        m = (3.141 * 10 ** 10 / temperature ** 3.444) * (math.log10(oil_density)) ** (10.313 * math.log10(temperature) - 36.447)
    
    elif method == "Labedi":
        m = 10 ** 9.224 * oil_density ** -4.7013 * temperature ** -0.6739
    
    elif method == "Egbogah & Ng":
        x = 10 ** (1.8653 - 2.5086 * 10 ** -2 * oil_density - 0.56411 * math.log10(temperature))
        m = 10 ** x - 1
    
    #elif method == "Al-Khafaji":      need review of source
    #    m = 10 ** (4.9563 - 0.00488 * Temperature * (oil_density + Temperature / 30 - 14.29) ** 2.709)
    
    elif method == "Petrosky":
        m = 2.3511 * 10 ** 7 / temperature ** 2.10255 * math.log10(oil_density) ** (4.59388 * math.log10(temperature) - 22.82792)
    
    elif method == "Kartoatmodjo & Schmidt":
        m = 1.6 * 10 ** 9 / temperature ** 2.8177 * math.log10(oil_density) ** (5.7526 * math.log10(temperature) - 26.9718)
    
    elif method == "De Ghetto":
        if oil_density <= 10:
            x = 10 ** (1.90296 - 1.2619 * 10 ** -2 * oil_density - 0.61748 * math.log10(temperature))
            m = 10 ** x - 1
        elif oil_density > 10 and oil_density <= 22.3:
            x = 10 ** (2.06492 - 1.79 * 10 ** -2 * oil_density - 0.70226 * math.log10(temperature))
            m = 10 ** x - 1
        elif oil_density > 22.3 and oil_density <= 31.1:
            m = 220.15 * 10 ** 9 / temperature ** 3.556 * math.log10(oil_density) ** (12.5428 * math.log10(temperature) - 45.7874)
        elif oil_density > 31.1:
            x = 10 ** (1.67083 - 1.7628 * 10 ** -2 * oil_density - 0.61304 * math.log10(temperature))
            m = 10 ** x - 1
        
    elif method == "Agip":
        x = 10 ** (1.8513 - 2.5548 * 10 ** -2 * oil_density - 0.56238 * math.log10(temperature))
        m = 10 ** x - 1
    
    elif method == "Fitzgerald":
        B = 0
        
    elif method == "Bennison":
        m = 10 ** (-0.8021 * oil_density + 23.8765) * temperature ** (0.31458 * oil_density - 9.21592)
    
    elif method == "Elsharkawy":
        x = 10 ** (2.16924 - 0.02525 * oil_density - 0.68875 * math.log10(temperature))
        m = 10 ** x - 1
    
    elif method == "Bergman":
        x = math.exp(22.33 - 0.194 * oil_density + 0.00033 * oil_density ** 2 - (3.2 - 0.0185 * oil_density) * math.log(temperature + 310))
        m = math.exp(x) - 1
    
    elif method == "Standing":
    #    oil_density = ucon(oil_density, "API", "sg (water)")
        x1 = 1 + 8.69 * math.log10((temperature + 459.67) / 559.67)
        x2 = 1 + 0.544 * math.log10((temperature + 459.67) / 559.67)
        x3 = -0.1285 * (2.87 * x1 - 1) * ucon(oil_density, "API", "sg (water)") / (2.87 * x1 - ucon(oil_density, "API", "sg (water)"))
        ro = ucon(oil_density, "API", "sg (water)") / (1 + 0.000321 * (temperature - 60) * 10 ** (0.00462 * oil_density))
        kw = calc_characterization_factor(oil_density)
        m = ro * 10 ** (1 / (x3 * (kw - 8.24 / ucon(oil_density, "API", "sg (water)")) + 1.639 * x2 - 1.059) - 2.17)
    
    elif method == "Dindoruk & Christman":
        a1 = 14.505357625
        a2 = -44.868655416
        a3 = 9.36579 * 10 ** 9
        a4 = -4.194017808
        a5 = -3.1461171 * 10 ** -9
        a6 = 1.517652716
        a7 = 0.010433654
        a8 = -0.00077688
        m = (a3 * temperature ** a4 * math.log10(oil_density) ** (a1 * math.log10(temperature) + a2)) / (a5 * pb ** a6 + a7 * Rs ** a8)
    else:
        m = "Select Method"

    oil_viscosity = m
    return oil_viscosity

def calc_bubble_point_oil_viscosity(method: str, MD, Rs):  #Gas saturated oil
#MD: dead oil viscosity, cP
#Recommended Chew & Connally, Beggs & Robinson
    
    if method == "Chew & Connally":
        a = 0.2 + 0.8 / 10 ** (0.00081 * Rs)
        B = 0.43 + 0.57 / 10 ** (0.00072 * Rs)
    
    elif method == "Beggs & Robinson":
        a = 10.715 / ((Rs + 100) ** 0.515)
        B = 5.44 / ((Rs + 150) ** 0.338)
    
    bubble_point_oil_viscosity = a * MD ** B
    return bubble_point_oil_viscosity

def calc_undersaturated_oil_viscosity(method: str, mb, pressure, pb):   #Undersaturated oil
    #mb: Gas-Sat oil viscosity, cP  Pressure: pressure, psi    pb: bubble-point pressure, psi
    #Recommended Beal, Kouzel, Vazquez & Beggs
    if method == "Beal":
        undersaturated_oil_viscosity = mb + (0.001 * (pressure - pb)) * (0.024 * mb ** 1.6 + 0.038 * mb ** 0.56)
    
    elif method == "Vazquez & Beggs":
        undersaturated_oil_viscosity = mb * (pressure / pb) ** (2.6 * pressure ** 1.187 * 10 ** (-3.9 * 10 ** -5 * pressure - 5))
    return undersaturated_oil_viscosity

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

def calc_relative_perm(Swi, Sor, krw, kro, nw, no, Sw, state):
    if state == 0:
        relative_perm = krw * ((Sw - Swi) / (1 - Swi - Sor)) ** nw
    else:
        relative_perm = kro * ((1 - Sor - Sw) / (1 - Swi - Sor)) ** no
    return relative_perm
    

def calc_GOR(gas_gravity, oil_density, pressure, temperature):
#gas_gravity:SG      oil_density:API      Pressure:psi        Temperature:F
    GOR = gas_gravity * (pressure / 18 * 10 ** (0.0125 * oil_density - 0.00091 * temperature)) ** 1.2048
    return GOR

def calc_characterization_factor(oil_density):
    Mo = 6084 / (oil_density - 5.9)
    characterization_factor = 4.5579 * Mo ** 0.15178 * ucon(oil_density, "API", "sg (water)") ** -0.84573
    return characterization_factor

def calc_non_HC_effect(oil_density, temperature, N2, CO2, H2S):
    #oil_density:API     Temperature:F
    pbN2 = 1.1585 + 2.86 * N2 - 1.07 * 10 ** -3
    #pbN2 = 1 + ((-2.65 * 10 ** -4 * oil_density + 5.5 * 10 ** -3) * Temperature + (0.0931 * oil_density - 0.8295)) * N2 + ((1.954 * 10 ** -11 * oil_density ** 4.699) * Temperature + (0.027 * oil_density - 2.366)) * N2 ** 2
    pbCO2 = 1 - 693.8 * CO2 * temperature ** -1.553
    pbH2S = 1 - (0.9035 + 0.0015 * oil_density) * H2S + 0.019 * (45 - oil_density) * H2S ** 2
    non_HC_effect = pbN2 * pbCO2 * pbH2S
    return non_HC_effect

def calc_gas_formation_factor(pressure, temperature, Zfactor):
    #output:cf/scf
    gas_formation_factor = 14.7 / ucon(60, "F", "R") * Zfactor * ucon(temperature, "F", "R") / pressure
    return gas_formation_factor

def calc_gas_gravity(properties, comp):
    i = 1
    x = 0
    for i in range(12):
        x = x + properties[i, 1] * comp[i]
    gas_gravity = x / 28.97
    return gas_gravity

def calc_gas_density(gas_gravity, pressure, temperature, ZFactor):
    #gas_density:lbm/ft3      Pressure:psi    Temperature:F
    gas_density = 28.967 / 10.732 * gas_gravity * pressure / (ZFactor * ucon(temperature, "F", "R"))
    return gas_density

def calc_viscosity_Carr(gas_gravity, pressure, temperature, N2, CO2, H2S):
    #viscosity calculations by Carr - Dempsey
    #temperature:degF
    temperature = ucon(temperature, "F", "R")
    
    Pc = 709.604 - 58.718 * gas_gravity
    Tc = 170.491 + 307.344 * gas_gravity
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
    
    visc_uncorrected = (temperature - 460) * (1.709 * 10 ** -5 - (2.062 * 10 ** -6) * gas_gravity) + (8.188 * 10 ** -3 - (6.15 * 10 ** -3) * math.log(gas_gravity))
    
    H2S_correction = H2S * (8.49 * 10 ** -3 * math.log(gas_gravity) + 3.73 * 10 ** -3)
    CO2_correction = CO2 * (9.08 * 10 ** -3 * math.log(gas_gravity) + 6.24 * 10 ** -3)
    N2_correction = N2 * (8.48 * 10 ** -3 * math.log(gas_gravity) + 9.59 * 10 ** -3)
    visc_corrected = visc_uncorrected + H2S_correction + CO2_correction + N2_correction
    
    viscosity = visc_corrected / Tr * math.exp((a0 + a1 + a2 + a3) + Tr * (a4 + a5 + a6 + a7) + Tr ** 2 * (a8 + a9 + a10 + a11) + Tr ** 3 * (a12 + a13 + a14 + a15))
    
    #viscosity_Carr = ucon(viscosity, "cP", Munit)
    viscosity_Carr = viscosity
    return viscosity_Carr

def calc_Zfactor(gas_gravity, pressure, temperature, N2, CO2, H2S):
    #gas_gravity:SG      Pressure:psi        Temperature:F
    #Z factor calculations by Hall-Yarborough with other gases present
    temperature = ucon(temperature, "F", "R")
    
    Pc = 678 - 50 * (gas_gravity - 0.5) - 206.7 * N2 + 440 * CO2 + 606.7 * H2S
    Tc = 326 + 315.7 * (gas_gravity - 0.5) - 240 * N2 - 83.3 * CO2 + 133.3 * H2S
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
    Zfactor = Z
    return Zfactor

def calc_Zfactor_nogas(gas_gravity, pressure, temperature):
#gas_gravity:SG      Pressure:psi        Temperature:F
#Z factor calculations by Hall-Yarborough with no additional gases
    temperature = ucon(temperature, "F", "R")
    
    Pc = 709.604 - 58.718 * gas_gravity
    Tc = 170.491 + 307.344 * gas_gravity
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
    Zfactor_nogas = Z
    return Zfactor_nogas






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

    fluid_compressibility = pressure * volume / compressibility
    return fluid_compressibility

def calc_water_hammer(tbg_od, tbg_id, solid_modulus, Poissons_ratio, pipe_state: str, fluid_rate, fluid_density, fluid_modulus):
    fluid_velocity = calc_fluid_velocity(fluid_rate, 0, tbg_id)
    pipe_thickness = (tbg_od - tbg_id) / 2

    if tbg_id / pipe_thickness > 10:   #thin walled
        if pipe_state == "Free":
            pipe_support_factor = 1
        elif pipe_state == "Upstream Anchored":
            pipe_support_factor = 1 - 0.5 * Poissons_ratio
        elif pipe_state == "All Anchored":
            pipe_support_factor = 1 - Poissons_ratio ** 2
        elif pipe_state == "Circular Flow":
            pipe_support_factor = 2 * pipe_thickness / tbg_id * (1 + Poissons_ratio)
        
    else:                        #thick walled
        if pipe_state == "Free":
            pipe_support_factor = 2 * pipe_thickness / tbg_id * (1 + Poissons_ratio) + tbg_id / (tbg_id + pipe_thickness)
        elif pipe_state == "Upstream Anchored":
            pipe_support_factor = 2 * pipe_thickness / tbg_id * (1 + Poissons_ratio) + tbg_id / (tbg_id + pipe_thickness) * (1 - Poissons_ratio / 2)
        elif pipe_state == "All Anchored":
            pipe_support_factor = 2 * pipe_thickness / tbg_id * (1 + Poissons_ratio) + tbg_id / (tbg_id + pipe_thickness) * (1 - Poissons_ratio ** 2)
        elif pipe_state == "Circular Flow":
            pipe_support_factor = 2 * pipe_thickness / tbg_id * (1 + Poissons_ratio)

    _CONST1 = 24.887
    _CONST2 = 1.615 * 10 ** -3
    sound_velocity = _CONST1 * math.sqrt(1 / (fluid_density * (1 / fluid_modulus + tbg_id * pipe_support_factor / (solid_modulus * pipe_thickness))))
    water_hammer = _CONST2 * fluid_density * sound_velocity * fluid_velocity
    return water_hammer




#*************************************************************************
#MATH FUNCTIONS
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
    interp = y
    return interp

#*************************************************************************
#Useful quick unit conversions
#NRe with API units - ft/s/in * ppg/cP multiply by 927.6866
#5.615 ft**3 in 1 bbl
#7.481 gal in 1 ft**3



#*************************************************************************
# Definitions of math functions
# from https://msdn.microsoft.com/en-us/library/csfk8t62(VS.85).aspx

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