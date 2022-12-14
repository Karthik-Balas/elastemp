"""
Elastemp is an open-source python package,distributed under MIT license, to compute quasi-harmonic approximated temperature dependent elastic constants of materials. 
author: Karthik Balasubramanian
"""

import argparse
from elastemp.input.make_elastic_input import make_elastic_wrapper,get_elastic_wrapper
from elastemp.input.make_dynamic_input import make_dynamic_wrapper, get_dynamic_wrapper,make_T_Zpe,get_T_Zpe
from elastemp.input.make_thermal_input import get_thermal_wrapper, get_thermal_constants, get_elastic_temp_plots

if __name__ == "__main__":
    
    symprec = 0.05
    angle_tolerance = 10

    parser = argparse.ArgumentParser()
    parser.add_argument("--operations", help="operations",\
                        choices=['make_elastic','get_elastic','make_bulk_dynamic','get_bulk_dynamic','make_temp_static','get_temp_static','make_dynamic','get_dynamic','get_thermal'])
    args = parser.parse_args()
    
    if args.operations=='make_elastic':
        make_elastic_wrapper(symprec=symprec,angle_tolerance=angle_tolerance)
    elif args.operations=='get_elastic':
        get_elastic_wrapper(make_plots=True)
    elif args.operations=='make_bulk_dynamic':
        make_dynamic_wrapper(calc_bulk_dynamic=True)
    elif args.operations=='get_bulk_dynamic':
        get_dynamic_wrapper(calc_bulk_dynamic=True)
    elif args.operations=='make_temp_static':
        make_T_Zpe()
    elif args.operations=='get_temp_static':
        get_T_Zpe()
    elif args.operations=='make_dynamic':
        make_dynamic_wrapper(calc_bulk_dynamic=False)
    elif args.operations=='get_dynamic':
        get_dynamic_wrapper(calc_bulk_dynamic=False)
    elif args.operations=='get_thermal':
        get_thermal_wrapper()
        get_thermal_constants()
        get_elastic_temp_plots()
    else:
        print('Operation not in choices')
