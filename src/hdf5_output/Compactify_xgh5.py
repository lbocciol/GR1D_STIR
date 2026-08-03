import h5py
import pandas as pd
import numpy as np
import os

# geometrized length unit (cm -> code units), same value as GR1D_module.F90
length_gf = 6.77140812e-06

exp_dictionary = {'Lambda_MLT': 0, \
              'M1_anue_aveenergy_fluid_rad': 0, \
              'M1_anue_enden_lab_rad': 0, \
              'M1_anue_enspectra_cen' : -35, \
              'M1_anue_enspectra_out' : -35, \
              'M1_anue_fluxden_lab_rad' : 0, \
              'M1_anue_fluxspectra_cen' : 0, \
              'M1_anue_fluxspectra_out' : -21, \
              'M1_anue_luminosity_fluid_rad' : 21, \
              'M1_anue_luminosity_lab_rad' : 21, \
              'M1_anue_numluminosity_fluid_rad' : 25, \
              'M1_anue_rmsenergy_fluid_rad' : 0, \
              'M1_eddingtonfactor_enweighted_anue' : 0, \
              'M1_eddingtonfactor_enweighted_nue' : 0, \
              'M1_eddingtonfactor_enweighted_nux' : 0, \
              'M1_eddingtonfactor_fluxweighted_anue' : 0, \
              'M1_eddingtonfactor_fluxweighted_nue' : 0, \
              'M1_eddingtonfactor_fluxweighted_nux' : 0, \
              'M1_fluxfactor_enweighted_anue' : 0, \
              'M1_fluxfactor_enweighted_nue' : 0, \
              'M1_fluxfactor_enweighted_nux' : 0, \
              'M1_fluxfactor_fluxweighted_anue' : 0, \
              'M1_fluxfactor_fluxweighted_nue' : 0, \
              'M1_fluxfactor_fluxweighted_nux' : 0, \
              'M1_nue_aveenergy_fluid_rad' : 0, \
              'M1_nue_enden_lab_rad' : 31, \
              'M1_nue_enspectra_cen' : -35, \
              'M1_nue_enspectra_out' : -35, \
              'M1_nue_fluxden_lab_rad' : 31, \
              'M1_nue_fluxspectra_cen' : -35, \
              'M1_nue_fluxspectra_out' : -35, \
              'M1_nue_luminosity_fluid_rad' : 21, \
              'M1_nue_luminosity_lab_rad' : 21, \
              'M1_nue_ng1_rad' : -15, \
              'M1_nue_numluminosity_fluid_rad' : 25, \
              'M1_nue_rmsenergy_fluid_rad' : 0, \
              'M1_nux_aveenergy_fluid_rad' : 0, \
              'M1_nux_enden_lab_rad' : 0, \
              'M1_nux_enspectra_cen' : -35, \
              'M1_nux_enspectra_out' : -35, \
              'M1_nux_fluxden_lab_rad' : 0, \
              'M1_nux_fluxspectra_cen' : 0, \
              'M1_nux_fluxspectra_out' : 0, \
              'M1_nux_luminosity_fluid_rad' : 21, \
              'M1_nux_luminosity_lab_rad' : 21, \
              'M1_nux_numluminosity_fluid_rad' : 25, \
              'M1_nux_rmsenergy_fluid_rad' : 0, \
              'Spectra/Energy_groups' : 0, \
              'Spectra/time_xg_spectra' : 0, \
              'W' : 0, \
              'X' : 0, \
              'alpha' : 0, \
              'alpha' : 0, \
              'buoyancy_turb_eps' : 0, \
              'bounce/entropy_at_bounce' : 0, \
              'bounce/mass_bary_at_bounce' : 33, \
              'bounce/rho_at_bounce' : 0, \
              'bounce/temperature_at_bounce' : 0, \
              'bounce/v_at_bounce' : 0, \
              'bounce/ye_at_bounce' : 0, \
              'capturing_factors' : 0, \
              'cs' : 0, \
              'depsdt' : 21, \
              'dissipated_turb_eps' : 0, \
              'dnupdr' : 0, \
              'dyedt_hydro' : 0, \
              'dyedt_neutrino' : 0, \
              'energy_nu' : 33, \
              'entropy' : 0, \
              'eps' : 0, \
              'eps_kin' : 0, \
              'mass_bary' : 33, \
              'mass_grav' : 33, \
              'nuchem' : 0, \
              'omega2_BV' : 0, \
              'press' : 34, \
              'press_nu' : 33, \
              'r' : 0, \
              'rho' : 0, \
              'shear_turb_eps': 0, \
              'temperature' : 0, \
              'time_xg' : 0, \
              'v' : 0, \
              'v1' : 0, \
              'v_turb' : 0, \
              'volume' : 26, \
              'xa' : 0, \
              'xabar' : 0, \
              'xh' : 0, \
              'xn' : 0, \
              'xp' : 0, \
              'xzbar' : 0, \
              'ye' : 0, \
              'ynu' : 0, \
              'omega' : 0, \
              'ToverW' : 0, \
              'vphi' : 0, \
              'vphi1' : 0, \
                }

spectra_xg = ['fluxspectra', 'enspectra', 'capturing_factors']

def Convert_xg_to_32bit_stacked(This_dir, nghosts=4):
    """Convert xg.h5 (groups /hydro and /M1 of (nspace, ntime) datasets, shared
    /time axis at the root, one-shot root datasets) to the flat 32-bit stacked
    layout of xg_files_32.h5.

    Output layout (unchanged from the pre-restructure version of this script):
      r, volume, <var> ...                 spatial variables, ghosts stripped
      Spectra/Energy_groups, Spectra/<spectrum>, Spectra/time_xg_spectra
      bounce/<var>_at_bounce, tbounce
      time_xg
    All values are scaled by 10**exp_dictionary[...] and stored as float32.
    """
    print(f'{This_dir} in progress')

    input_h5 = os.path.join(This_dir, 'xg.h5')
    output_h5 = os.path.join(This_dir, 'xg_files_32.h5')

    try:
        with h5py.File(input_h5, 'r') as f_in, h5py.File(output_h5, 'w') as f_out:

            if 'time' not in f_in:
                print(f"No /time dataset found in {input_h5}")
                return

            time = f_in['time'][:]
            nT = len(time)

            f_out.create_group('Spectra')

            # ==========================================
            # 1. TIME AXIS (shared by grid and spectra variables)
            # ==========================================
            for out_path in ['time_xg', 'Spectra/time_xg_spectra']:
                exponent = exp_dictionary.get(out_path, 0)
                f_out.create_dataset(out_path, data=(time / 10.0**exponent).astype(np.float32))

            # ==========================================
            # 2. ROOT-LEVEL ONE-SHOT DATASETS
            # ==========================================
            group_names = [k for k in f_in.keys() if isinstance(f_in[k], h5py.Group)]
            root_vars = [k for k in f_in.keys()
                         if k not in group_names and k not in ('time', 'metadata')]

            if any(v.endswith('_at_bounce') for v in root_vars):
                f_out.create_group('bounce')

            for var in root_vars:
                if var == 'radius':
                    out_name = 'r'
                elif var == 'neutrino_energies':
                    out_name = 'Energy_groups'
                else:
                    out_name = var

                out_path = out_name
                if out_name == 'Energy_groups':
                    out_path = f"Spectra/{out_name}"
                elif out_name.endswith('_at_bounce'):
                    out_path = f"bounce/{out_name}"

                exponent = exp_dictionary.get(out_path, exp_dictionary.get(out_name, 0))
                scale = 10.0 ** exponent

                data = f_in[var][:]

                if out_name == 'volume':
                    # old outputs stored volume in code units; convert to cm^3
                    if np.max(data) > 1e25:
                        volume_conv = 1.0
                    else:
                        volume_conv = 1.0 / length_gf**3
                    data = data * volume_conv

                # Ghost slicing for spatial root variables (length n1); leave
                # energy grids and scalars (tbounce) alone
                is_spatial = (out_name in ('r', 'volume')) or out_name.endswith('_at_bounce')
                if is_spatial and nghosts > 0 and data.ndim > 0 and len(data) > 2 * nghosts:
                    data = data[nghosts:-nghosts, ...]

                data_scaled = (data / scale).astype(np.float32)
                f_out.create_dataset(out_path, data=data_scaled)

            # ==========================================
            # 3. TIME-DEPENDENT VARIABLES (/hydro and /M1)
            # ==========================================
            # Each dataset is already stacked as (nspace, nT); just strip
            # ghosts, scale, and downcast.
            for group_name in ('hydro', 'M1'):
                if group_name not in f_in:
                    continue

                for var in f_in[group_name].keys():
                    out_name = var
                    is_spectra_var = any(check_spectra in out_name for check_spectra in spectra_xg)

                    out_path = f"Spectra/{out_name}" if is_spectra_var else out_name

                    exponent = exp_dictionary.get(out_path, exp_dictionary.get(out_name, 0))
                    scale = 10.0 ** exponent

                    data = f_in[group_name][var][:]

                    # Spectra live on the energy grid and are NOT ghost-sliced
                    if not is_spectra_var and nghosts > 0 and data.ndim > 0:
                        data = data[nghosts:-nghosts, ...]

                    if data.shape[-1] != nT:
                        print(f"Warning: {group_name}/{var} has {data.shape[-1]} dumps, "
                              f"expected {nT}; writing as-is")

                    f_out.create_dataset(out_path, data=(data / scale).astype(np.float32))

        # Write dictionary to CSV
        df = pd.DataFrame.from_dict([exp_dictionary])
        df.to_csv(os.path.join(This_dir, 'xg_files_dictionary.csv'), header=True, index=False, mode='w')

    except Exception as e:
        print(f'{This_dir} big error: {e}')
        if os.path.exists(output_h5):
            os.remove(output_h5)
        return

    print(f'{This_dir} successful')
    return f'{This_dir} successful'


def Convert_burn_to_32bit(This_dir):
    """Compress burn.h5 to burn_files_32.h5: identical structure and dataset
    names, with every float dataset downcast to float32 (no scaling, no ghost
    stripping). Group/root attributes (e.g. /metadata) are copied verbatim.
    """
    print(f'{This_dir} burn.h5 in progress')

    input_h5 = os.path.join(This_dir, 'burn.h5')
    output_h5 = os.path.join(This_dir, 'burn_files_32.h5')

    try:
        with h5py.File(input_h5, 'r') as f_in, h5py.File(output_h5, 'w') as f_out:

            def copy_attrs(src, dst):
                for key, val in src.attrs.items():
                    dst.attrs[key] = val

            def copy_item(name, obj):
                if isinstance(obj, h5py.Group):
                    grp = f_out.create_group(name)
                    copy_attrs(obj, grp)
                else:
                    data = obj[:]
                    if np.issubdtype(data.dtype, np.floating):
                        data = data.astype(np.float32)
                    dset = f_out.create_dataset(name, data=data)
                    copy_attrs(obj, dset)

            copy_attrs(f_in, f_out)
            f_in.visititems(copy_item)

    except Exception as e:
        print(f'{This_dir} burn.h5 big error: {e}')
        if os.path.exists(output_h5):
            os.remove(output_h5)
        return

    print(f'{This_dir} burn.h5 successful')
    return f'{This_dir} burn.h5 successful'

