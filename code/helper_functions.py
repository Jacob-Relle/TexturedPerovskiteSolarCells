from typing import Literal
import numpy as np
import sys, os
import pandas as pd
from decimal import Decimal
from pathlib import Path



def get_spectrum(density: Literal['full','thin','single'] , wavelength = 500):
    assert density in ['full', 'thin', 'single']

    # Generate the AM1.5 Sunspectrum
    working_dir = Path.cwd().parent
    spectral_data = pd.read_excel(f'{working_dir}/data/AM1.5.xls')
    if density == 'full':
        spectrum = np.around(np.array(spectral_data.iloc[41:1042,0].tolist()) * 1e-9,15)
        irradiance = np.array(spectral_data.iloc[41:1042,2].tolist())

    elif density == 'thin':
        spectrum = np.arange(300,910,10)*1e-9
        irradiance = np.array(spectral_data.iloc[41:742,2].tolist())
        irradiance = np.append(irradiance[0:200:20],irradiance[200::10])

    elif density == 'single':
        if Decimal(wavelength).as_tuple().exponent < -1:
            wavelength *= 1e9
        assert wavelength >= 300 and wavelength <= 1500
        irradiance = spectral_data.loc[spectral_data['ASTM G173-03 Reference Spectra Derived from SMARTS v. 2.9.2'] == wavelength].iloc[0,2]
        spectrum = np.array([wavelength])*1e-9

    return spectrum, irradiance

def koch_curve(order):
    """
    Generate the coordinates for a Koch curve starting from a unit intervall.
    
    Parameters:
    order (int): The recursive depth of the fraktal.
    
    Returns:
    list: List of complex numbers representing the vertices of the curve.
    """
    if order == 0:
        return np.array([0, 1], dtype=complex)
    else:
        # Recursively generate the smaller curves
        prev_curve = koch_curve(order - 1)
        new_curve = []

        for i in range(len(prev_curve) - 1):
            p1 = prev_curve[i]
            p2 = prev_curve[i + 1]
            segment = p2 - p1
            new_curve.extend([
                p1,
                p1 + segment / 3,
                p1 + segment / 3 + segment * (np.exp(1j * np.pi / 3)) / 3,
                p1 + 2 * segment / 3
            ])
        
        new_curve.append(prev_curve[-1])
        return np.array(new_curve)

def fresnel_s(n_1, n_2, theta=0):

    x = n_2*np.sqrt(1 - (n_1/n_2 * np.sin(theta))**2)
    num =  n_1*np.cos(theta) - x
    denum =  n_1*np.cos(theta) + x
    R_s = np.abs( num / denum )**2

    return R_s 

def fresnel_p(n_1, n_2, theta=0):

    x = n_1*np.sqrt(1 - (n_1/n_2 * np.sin(theta))**2)
    num = x - n_2*np.cos(theta)
    denum = x + n_2*np.cos(theta)
    R_p = np.abs( num / denum )**2

    return R_p
