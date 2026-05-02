#!/usr/bin/env python3
"""
Comprehensive HDF5 to Text Converter for GR1D

Converts xg.h5 and dat.h5 back to original output.F90 text format.
Uses neutrino_energies from HDF5 as energy coordinate for spectra.
Uses radius from HDF5 as spatial coordinate for grid data.
"""

import h5py
import numpy as np
import sys
from pathlib import Path


def format_fortran_e18_9(value):
    """Format a single value in Fortran 1P20E18.9 format (1P E18.9)
    
    Fortran format 1P20E18.9 means:
    - 1P: scale by 10^0 with leading significant digit
    - E18.9: E notation with 18 total chars, 9 decimal places
    - Example: 1.234567E+00 (total 18 chars)
    """
    if abs(value) < 1.0e-90:
        value = 0.0
    # Format to match Fortran 1PE18.9: E notation with 18 chars total, 9 decimals
    # This produces output like: " 1.234567890E+00"
    return f"{value:18.9E}"


def write_xg_line(out, coord, value):
    """Write a single (coordinate, value) line in output_single format
    
    Matches Fortran: write(666,"(1P20E18.9)") coordinate, value
    """
    if abs(coord) < 1.0e-90:
        coord = 0.0
    if abs(value) < 1.0e-90:
        value = 0.0
    out.write(f"{format_fortran_e18_9(coord)}{format_fortran_e18_9(value)}\n")


def write_spectrum_line(out, energy, value):
    """Write a single (energy, value) line in output_spectra format
    
    Matches Fortran: write(666,"(1P20E18.9)") energy, value
    """
    if abs(energy) < 1.0e-90:
        energy = 0.0
    if abs(value) < 1.0e-90:
        value = 0.0
    out.write(f"{format_fortran_e18_9(energy)}{format_fortran_e18_9(value)}\n")


def write_scalar_line(out, time, *values):
    """Write a time + values line in output_scalar/output_many_scalars format
    
    Matches Fortran: write(666,"(1P20E18.9)") time, var(1:n) OR
                     write(666,"(1P256E18.9)") time, var(1:n)
    """
    if abs(time) < 1.0e-90:
        time = 0.0
    line = format_fortran_e18_9(time)
    for val in values:
        if abs(val) < 1.0e-90:
            val = 0.0
        line += format_fortran_e18_9(val)
    out.write(line + "\n")


def convert_xg_file(hdf5_file, output_dir):
    """
    Convert xg.h5 to individual .xg text files.
    
    For each variable in output_#### groups, create .xg file with:
    - Time header: "Time = <value>
    - Data rows: radius value (or mass value if using mass coordinate)
    - Blank lines between timesteps
    """
    print(f"Converting {hdf5_file.name} to text format...")
    
    with h5py.File(hdf5_file, 'r') as f:
        # Get all output groups
        output_groups = sorted([k for k in f.keys() if k.startswith('output_')])
        
        if not output_groups:
            print("No output groups found in xg.h5")
            return
        
        # Get metadata to determine coordinate system
        first_group = f[output_groups[0]]
        
        # Check for coordinate datasets at root level (new structure)
        has_radius = 'radius' in f
        has_neutrino_energies = 'neutrino_energies' in f
        has_mass = 'mass' in first_group
        
        # Get variable names (exclude coordinate/metadata datasets)
        exclude_keys = {'time', 'radius', 'mass', 'volume', 'neutrino_energies'}
        variables = sorted([k for k in first_group.keys() if k not in exclude_keys and not 
                            'bounce' in k])
        
        print(f"Found {len(output_groups)} timesteps with {len(variables)} variables")
        print(f"Coordinate system: {'radius' if has_radius else 'mass' if has_mass else 'grid index'}")
        
        # Read neutrino_energies from root dataset if it exists (new structure, for spectra conversion)
        neutrino_energies = None
        if has_neutrino_energies:
            neutrino_energies = f['neutrino_energies'][:]
            print(f"Found neutrino_energies coordinate ({len(neutrino_energies)} groups)")
        
        # Get radius for grid data (from root level in new structure)
        radius = None
        if has_radius:
            radius = f['radius'][:]
        
        # Classify variables by checking dimension size
        # Grid variables should match radius/mass size; spectrum variables use number_groups size
        grid_size = len(radius) if radius is not None else None
        spectrum_size = len(neutrino_energies) if neutrino_energies is not None else None
        
        # Classify each variable based on its actual data size
        grid_vars = []
        spectrum_vars = []
        
        for var in variables:
            first_group = f[output_groups[0]]
            var_data = first_group[var][:]
            var_size = len(var_data)
            
            # Check if it matches neutrino_energies dimension (spectrum)
            if spectrum_size is not None and var_size == spectrum_size:
                spectrum_vars.append(var)
            # Otherwise treat as grid data
            else:
                grid_vars.append(var)
        
        print(f"  {len(grid_vars)} grid variables, {len(spectrum_vars)} spectrum variables")
        
        # Process bounce variables
        bounce_vars = [var for var in variables if 'at_bounce' in var]
        if bounce_vars:
            tbounce = f['tbounce'][0]
            for var in bounce_vars:
                xg_file = output_dir / f"{var}.xg"
                print(f"  Writing {var}.xg...", end='', flush=True)
                
                with open(xg_file, 'w') as out:
                    group = f[output_group]
                    
                    # Write time header in Fortran list-directed format
                    out.write(f' "Time = {tbounce:25.16E}\n')

                    # Get variable data
                    data = group[var][:]
                    
                    # Get coordinate system (radius is at root level in new structure)
                    if 'mass' in var:
                        coords = radius
                    else:
                        coords = group['mass'][:]
                    
                    # Write data in columns
                    for i in range(len(data)):
                        write_xg_line(out, coords[i], data[i])
                    
                    # Write blank lines between outputs
                    out.write("\n\n")
                
                print(" done")

        # Process grid data variables
        if grid_vars:
            for var in grid_vars:
                xg_file = output_dir / f"{var}.xg"
                print(f"  Writing {var}.xg...", end='', flush=True)
                
                with open(xg_file, 'w') as out:
                    for output_group in output_groups:
                        group = f[output_group]
                        time = group['time'][0]
                        
                        # Write time header in Fortran list-directed format
                        if time == 0:
                            out.write(f' "Time =    0.0000000000000000     \n')
                        else:
                            out.write(f' "Time = {time:25.16E}\n')

                        # Get variable data
                        data = group[var][:]
                        
                        # Get coordinate system (radius is at root level in new structure)
                        if has_radius:
                            coords = radius
                        elif has_mass and 'mass' in group:
                            coords = group['mass'][:]
                        else:
                            coords = np.arange(len(data))
                        
                        # Write data in columns
                        for i in range(len(data)):
                            write_xg_line(out, coords[i], data[i])
                        
                        # Write blank lines between outputs
                        out.write("\n\n")
                
                print(" done")
        
        # Process spectrum variables
        if spectrum_vars:
            for var in spectrum_vars:
                xg_file = output_dir / f"{var}.xg"
                print(f"  Writing {var}.xg...", end='', flush=True)
                
                with open(xg_file, 'w') as out:
                    for output_group in output_groups:
                        group = f[output_group]
                        time = group.attrs.get('time', 0.0)
                        
                        # Write time header in Fortran list-directed format
                        out.write(f' "Time = " {time}\n')
                        
                        # Get spectrum data
                        spectrum_data = group[var][:]
                        
                        # Use neutrino_energies as energy coordinate if available
                        if neutrino_energies is not None:
                            for i in range(len(neutrino_energies)):
                                write_spectrum_line(out, neutrino_energies[i], spectrum_data[i])
                        else:
                            # Fallback to index-based coordinates
                            for i in range(len(spectrum_data)):
                                write_spectrum_line(out, float(i), spectrum_data[i])
                        
                        # Write blank lines between outputs
                        out.write("\n\n")
                
                print(" done")
        
        # Process root-level coordinate datasets (volume, etc.)
        print(f"\nProcessing root-level datasets...")
        
        if 'volume' in f:
            print(f"  Writing volume.xg...", end='', flush=True)
            volume = f['volume'][:]
            xg_file = output_dir / "volume.xg"
            
            with open(xg_file, 'w') as out:
                # Write volume only ONCE, since it's a static root-level dataset
                # We pull the first output group just to evaluate coordinates if needed
                first_group = f[output_groups[0]]
                
                # Write the static time header using your 0.0 format
                out.write(' "Time =    0.0000000000000000    \n')
                
                # Get coordinate system
                if has_radius:
                    coords = radius
                elif has_mass and 'mass' in first_group:
                    coords = first_group['mass'][:]
                else:
                    coords = np.arange(len(volume))
                
                # Write data in columns
                for i in range(len(volume)):
                    write_xg_line(out, coords[i], volume[i])
                
                # Write final blank lines
                out.write("\n\n")
            
            print(" done")
        
        print(f"Successfully converted {len(variables)} variables to .xg files")

def convert_dat_file(hdf5_file, output_dir):
    """
    Convert dat.h5 to individual .dat text files.
    
    For each scalar time series in /scalars group:
    - Create .dat file with: time + value(s) per row
    - Handle both 1D (single value) and 2D (multiple values) data
    - Time is stored in time_c or similar dataset that matches data length
    - Automatically adjusts time arrays for post-bounce variables
    - Writes tbounce as a single number in tbounce.dat
    """
    print(f"Converting {hdf5_file.name} to text format...")
    
    with h5py.File(hdf5_file, 'r') as f:
        if 'scalars' not in f:
            print("No /scalars group found in dat.h5")
            return
        
        scalars_group = f['scalars']
        
        # Get all variables except 'metadata' and explicitly exclude 'accretion_radii'
        exclude_keys = {'metadata', 'accretion_radii'}
        all_variables = sorted([k for k in scalars_group.keys() if k not in exclude_keys])
        
        # Find time array - look for time_c or any variable that could be time
        time_var_name = 'time'
        if time_var_name not in scalars_group:
            print(f"Error: No time dataset found. Available datasets: {all_variables}")
            return
            
        time = scalars_group[time_var_name][:]
        print(f"Found time in '{time_var_name}' with {len(time)} timesteps")
        
        # --- NEW: Extract tbounce, write to tbounce.dat, and exclude it from the loop ---
        tbounce_val = None
        if 'tbounce' in scalars_group:
            tbounce_val = scalars_group['tbounce'][0]
        elif 'tbounce' in f:
            tbounce_val = f['tbounce'][0]
            
        if tbounce_val is not None:
            print(f"Found tbounce = {tbounce_val}. Writing tbounce.dat and preparing to adjust arrays.")
            with open(output_dir / "tbounce.dat", 'w') as tb_out:
                tb_out.write(f"{tbounce_val:18.9E}\n")
            
            # Make sure we don't try to process tbounce as a time series below
            exclude_keys.add('tbounce')

        # Get data variables (exclude time and already excluded keys)
        exclude_keys.add(time_var_name)
        variables = sorted([k for k in all_variables if k not in exclude_keys])
        
        print(f"Found {len(variables)} scalar time series")
        
        # For each variable, create a .dat file
        for var in variables:
            dat_file = output_dir / f"{var}.dat"
            print(f"  Writing {var}.dat...", end='', flush=True)
            
            with open(dat_file, 'w') as out:
                data = scalars_group[var][:]
                
                # Determine data length along time axis
                data_len = len(data) if data.ndim == 1 else data.shape[1]
                current_time = time
                
                # Check for mismatch and apply post-bounce time slicing
                if data_len != len(time):
                    if tbounce_val is not None:
                        post_bounce_time = time[time >= tbounce_val]
                        
                        if data_len == len(post_bounce_time):
                            current_time = post_bounce_time
                        elif data_len < len(time):
                            # Fallback: Slicing from the tail is bulletproof for post-bounce arrays
                            current_time = time[-data_len:]
                        else:
                            print(f"\nWarning: Length mismatch for {var}: data has {data_len} entries, time has {len(time)}")
                            continue
                    else:
                        print(f"\nWarning: Length mismatch for {var}: data has {data_len} entries, time has {len(time)} (No tbounce found to adjust)")
                        continue
                
                # Handle 1D and 2D data using the properly sized current_time array
                if data.ndim == 1:
                    # 1D scalar time series (single value per timestep)
                    for t, val in zip(current_time, data):
                        write_scalar_line(out, t, val)
                        
                elif data.ndim == 2:
                    # 2D array stored as (components, time) in HDF5
                    # Need to transpose to (time, components) for output
                    for t_idx, t in enumerate(current_time):
                        values = data[:, t_idx]  # Get all components for this timestep
                        write_scalar_line(out, t, *values)
                        
                else:
                    print(f"\nWarning: Unexpected data dimensionality for {var}: {data.ndim}")
                    continue
            
            print(" done")
        
        # Handle special cases with headers
        # Read accretion_radii (Fortran implies root dataset, but checking scalars group too just in case)
        radii_data = None
        if 'accretion_radii' in f:
            radii_data = f['accretion_radii'][:]
        elif 'accretion_radii' in scalars_group:
            radii_data = scalars_group['accretion_radii'][:]
            
        if radii_data is not None:
            # Flatten the 2D (11, 1) array into a 1D sequence and cast to float
            header_str = "#Radii: " + "".join(f"{float(r):18.9E}" for r in radii_data.flatten()) + "\n"
            
            # Prepend the header to specific files
            target_files = ["accretion_rates.dat", "accreted_mass.dat"]
            for target in target_files:
                target_path = output_dir / target
                if target_path.exists():
                    print(f"  Prepending header to {target}...")
                    # Read existing content
                    with open(target_path, 'r') as file:
                        original_content = file.read()
                    # Rewrite with header
                    with open(target_path, 'w') as file:
                        file.write(header_str)
                        file.write(original_content)
                else:
                    print(f"  Warning: {target} was not generated; skipping header injection.")

        print(f"Successfully converted {len(variables)} scalar time series to .dat files")

def write_scalar_line(file_obj, time_val, *values):
    """Helper to write time and values to file in standard format"""
    # Assuming standard scientific notation formatting for the columns
    line = f"{time_val:18.9E} " + " ".join(f"{v:18.9E}" for v in values) + "\n"
    file_obj.write(line)

def main():
    if len(sys.argv) != 2:
        print("Usage: python hdf5_to_text.py <full_path_to_hdf5_file>")
        print("  Converts xg.h5 or dat.h5 to original output.F90 text format")
        print("  Output files are saved in the same directory as the input file")
        sys.exit(1)
    
    # Resolve creates an absolute path
    input_file = Path(sys.argv[1]).resolve()
    output_dir = input_file.parent
    
    # Check if file exists
    if not input_file.exists():
        print(f"Error: File {input_file} not found")
        sys.exit(1)
    
    # Ensure the path points to a file, not a directory
    if not input_file.is_file():
        print(f"Error: {input_file} is a directory. Please provide the path to a file.")
        sys.exit(1)
    
    # Check if HDF5 file
    if not h5py.is_hdf5(input_file):
        print(f"Error: {input_file} is not a valid HDF5 file")
        sys.exit(1)
    
    # Convert based on filename
    if input_file.name == 'xg.h5':
        convert_xg_file(input_file, output_dir)
    elif input_file.name == 'dat.h5':
        convert_dat_file(input_file, output_dir)
    else:
        print(f"Error: Unrecognized file '{input_file.name}'. File must be named exactly 'xg.h5' or 'dat.h5'")
        sys.exit(1)
    
    print(f"\nConversion complete! Files written to {output_dir}")


if __name__ == '__main__':
    main()
