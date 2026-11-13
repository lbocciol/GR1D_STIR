#!/usr/bin/env python3
"""
HDF5 to Text Converter for GR1D

Converts xg.h5 and dat.h5 back to the legacy text output of output.F90 (.xg and .dat files),
i.e. what GR1D writes when compiled with HAVE_HDF5_OUTPUT=0.

xg.h5 layouts (both supported):
  layout 1 (older): one group /output_NNNNN per output, one dataset per variable inside
  layout 2 (/metadata attribute xg_layout = 2): one dataset per variable at the root, shape
           (n, nT) with n = n1 or number_groups, and /time of shape (nT); only the first
           len(time) outputs are complete
Grid variables use /radius as coordinate, spectra use /neutrino_energies. The bounce profiles
(*_at_bounce) and tbounce.dat are rebuilt from the root datasets written at bounce.

Usage: python hdf5_to_text.py <outdir | path/xg.h5 | path/dat.h5>
  with a folder both dat.h5 and xg.h5 in it are converted; text files go next to the h5 files
"""

import sys
import math
from pathlib import Path

import h5py
import numpy as np

MASS_GF = 5.02765209e-34  # GR1D_module mass_gf: grams -> code units (solar masses)


# ---------------------------------------------------------------------------------------------
# Fortran formatting

def fortran_e(value, width=18, digits=9, scale1p=True):
    """Fortran 1PEw.d (scale1p, d.ddd) or Ew.d (0.ddd) edit descriptor, correctly rounded"""
    if not math.isfinite(value):
        return ('NaN' if math.isnan(value) else ('-Infinity' if value < 0 else 'Infinity')).rjust(width)
    sig = digits + 1 if scale1p else digits
    mant, exp = f'{abs(value):.{sig - 1}e}'.split('e')  # correctly rounded significant digits
    exp = int(exp)
    if value == 0.0:
        exp = 0 if scale1p else -1
    if scale1p:
        m = mant
    else:
        m = '0.' + mant.replace('.', '')
        exp += 1
    e = f'E{exp:+03d}' if abs(exp) <= 99 else f'{exp:+04d}'
    return (('-' if value < 0 or (value == 0.0 and math.copysign(1, value) < 0) else '') + m + e).rjust(width)


def xg_value(value):
    # output_single / output_spectra write 0 for |value| < 1e-90
    return 0.0 if abs(value) < 1.0e-90 else value


def fortran_list_real(value):
    """ifort list-directed output of a real*8 (as in write(666,*) '"Time = ',time):
    fixed notation with 15 significant digits for 0.1 <= |x| < 1e15, else E with 16"""
    if value == 0.0:
        return '  0.000000000000000E+000'
    exp15 = int(f'{value:.14e}'.split('e')[1])  # exponent after rounding to 15 digits
    rounded = float(f'{value:.14e}')
    if 0.1 <= abs(rounded) and abs(value) < 1.0e15:
        decimals = max(15 - (exp15 + 1), 0)
        text = f'{rounded:.{decimals}f}' + ('.' if decimals == 0 else '')
        return text.rjust(19) + '     '
    mant, exp = f'{value:.15e}'.split('e')  # E notation: 16 significant digits
    return ('  ' if value > 0 else ' ') + f'{mant}E{int(exp):+04d}'


def header_line(time):
    return ' "Time = ' + fortran_list_real(time) + '\n'


BLANK_LINE = '  \n'  # write(666,*) " "


def write_block(out, time, coords, values):
    out.write(header_line(time))
    out.write(''.join(fortran_e(c) + fortran_e(xg_value(v)) + '\n' for c, v in zip(coords, values)))
    out.write(BLANK_LINE)
    out.write(BLANK_LINE)


def scalar_line(time, values):
    return fortran_e(time) + ''.join(fortran_e(v) for v in values) + '\n'


# ---------------------------------------------------------------------------------------------
# xg.h5 reading (both layouts)

def xg_layout(f):
    if 'metadata' in f and 'xg_layout' in f['metadata'].attrs:
        return int(np.ravel(f['metadata'].attrs['xg_layout'])[0])
    return 1 if any(k.startswith('output_') for k in f) else 2


def xg_series(f):
    """Returns (times, {variable: array (nT, n)}) for both layouts, outputs in time order"""
    if xg_layout(f) == 1:
        groups = [k for k in f if k.startswith('output_')]
        times = np.array([float(np.ravel(f[g]['time'][()])[0]) for g in groups])
        order = np.argsort(times, kind='stable')
        groups, times = [groups[i] for i in order], times[order]
        names = [k for k in f[groups[0]] if k != 'time']
        data = {n: np.array([f[g][n][()] for g in groups]) for n in names}
        return times, data
    times = f['time'][()]
    nT = times.size
    data = {}
    for name, obj in f.items():
        if isinstance(obj, h5py.Dataset) and obj.ndim == 2:
            data[name] = obj[:, :nT].T  # (n, nT) on file -> (nT, n)
    return times, data


def convert_xg_file(hdf5_file, output_dir):
    print(f'Converting {hdf5_file} ...')
    # whole file into memory: one sequential read instead of many small ones
    with h5py.File(hdf5_file, 'r', driver='core', backing_store=False) as f:
        print(f'  layout {xg_layout(f)}')
        radius = f['radius'][()]
        energies = f['neutrino_energies'][()] if 'neutrino_energies' in f else None
        times, data = xg_series(f)
        print(f'  {times.size} outputs, {len(data)} variables')

        for name, values in sorted(data.items()):
            if energies is not None and values.shape[1] == energies.size:
                coords = energies
            elif values.shape[1] == radius.size:
                coords = radius
            else:
                print(f'  warning: {name} has {values.shape[1]} points per output, skipped')
                continue
            with open(output_dir / f'{name}.xg', 'w') as out:
                for t, v in zip(times, values):
                    write_block(out, t, coords, v)

        # static profile
        if 'volume' in f:
            with open(output_dir / 'volume.xg', 'w') as out:
                write_block(out, 0.0, radius, f['volume'][()])

        # bounce profiles: mass_bary vs radius, the others vs enclosed mass in solar masses
        if 'tbounce' in f:
            t_bounce = float(np.ravel(f['tbounce'][()])[0])
            mass_code = f['mass_bary_at_bounce'][()] * MASS_GF if 'mass_bary_at_bounce' in f else None
            for name in sorted(k for k in f if k.endswith('_at_bounce')):
                coords = radius if name == 'mass_bary_at_bounce' else mass_code
                if coords is None:
                    print(f'  warning: no mass_bary_at_bounce, {name}.xg skipped')
                    continue
                with open(output_dir / f'{name}.xg', 'w') as out:
                    write_block(out, t_bounce, coords, f[name][()])
            if 'rho_at_bounce' not in f:
                print('  note: rho_at_bounce not stored in this file (older GR1D), rho_at_bounce.xg not written')
            if 'tbounce_dat' in f:
                with open(output_dir / 'tbounce.dat', 'w') as out:
                    out.write(''.join(fortran_e(v, 27, 18, scale1p=False) + '\n' for v in f['tbounce_dat'][()]))
            else:
                print('  note: tbounce_dat not stored in this file (older GR1D), tbounce.dat only from dat.h5')
    print(f'  done: {len(data)} .xg files')


# ---------------------------------------------------------------------------------------------
# dat.h5

def convert_dat_file(hdf5_file, output_dir):
    print(f'Converting {hdf5_file} ...')
    with h5py.File(hdf5_file, 'r', driver='core', backing_store=False) as f:
        s = f['scalars']
        time = s['time'][()]
        skip = {'time', 'tbounce', 'accretion_radii'}
        n_files = 0
        for name in sorted(k for k in s if k not in skip):
            data = s[name][()]
            # series written only after bounce are aligned with the end of time
            n = data.shape[-1]
            t = time[time.size - n:]
            with open(output_dir / f'{name}.dat', 'w') as out:
                if name in ('accretion_rates', 'accreted_mass') and 'accretion_radii' in s:
                    # header written by analysis.F90 when the run starts
                    radii = s['accretion_radii'][()][:, 0]
                    out.write('#Radii: ' + ''.join(fortran_e(r) for r in radii) + '\n')
                if data.ndim == 1:
                    out.write(''.join(scalar_line(ti, [vi]) for ti, vi in zip(t, data)))
                else:
                    out.write(''.join(scalar_line(ti, data[:, i]) for i, ti in enumerate(t)))
            n_files += 1

        # tbounce.dat: the full legacy version comes from xg.h5 (tbounce_dat); only t_bounce here
        if 'tbounce' in s and not (output_dir / 'tbounce.dat').exists():
            with open(output_dir / 'tbounce.dat', 'w') as out:
                out.write(fortran_e(float(np.ravel(s['tbounce'][()])[0]), 27, 18, scale1p=False) + '\n')
    print(f'  done: {n_files} .dat files')


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    path = Path(sys.argv[1]).resolve()
    if path.is_dir():
        files = [path / 'xg.h5', path / 'dat.h5']  # xg first: it has the full tbounce.dat
    else:
        files = [path]
    for file in files:
        if not file.is_file() or not h5py.is_hdf5(file):
            sys.exit(f'Error: {file} is not an HDF5 file')
        if file.name == 'xg.h5':
            convert_xg_file(file, file.parent)
        elif file.name == 'dat.h5':
            convert_dat_file(file, file.parent)
        else:
            sys.exit(f"Error: unrecognized file '{file.name}', expected xg.h5 or dat.h5")


if __name__ == '__main__':
    main()
