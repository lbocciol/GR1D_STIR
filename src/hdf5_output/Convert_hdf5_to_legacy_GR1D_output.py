#!/usr/bin/env python3
import os
import sys
import subprocess
import multiprocessing
import functools
from pathlib import Path
from Compactify_xgh5 import Convert_xg_to_32bit_stacked, Convert_burn_to_32bit

def run_conversions(folder, dir_dict):
    """Worker function to execute the conversion script on a specific folder"""
    script_path = dir_dict['script_path']

    names_to_convert = dir_dict['files_to_convert']

    # Check for the specific files the HDF5 converter expects
    for h5_name in names_to_convert:
        target_file = os.path.join(folder, h5_name)
        if os.path.exists(target_file):
            print(f"--- Processing: {h5_name} in {os.path.basename(folder)} ---")
            try:
                # Execute hdf5_to_text.py <path_to_h5>
                subprocess.run([sys.executable, script_path, target_file], check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error converting {target_file}: {e}")

    # 32-bit stacked xg output
    if 'xg_files_32.h5' in os.listdir(folder):
        print('xg_files_32.h5 already present in ' + folder)
    elif os.path.exists(os.path.join(folder, 'xg.h5')):
        Convert_xg_to_32bit_stacked(folder, nghosts=dir_dict['nghosts'])

    # 32-bit copy of the burn output (burn.h5 is HDF5-only, no text conversion)
    if os.path.exists(os.path.join(folder, 'burn.h5')):
        if 'burn_files_32.h5' in os.listdir(folder):
            print('burn_files_32.h5 already present in ' + folder)
        else:
            Convert_burn_to_32bit(folder)

# --- Configuration ---
main_dir_path = '/home/lbocciol/GR1D_burn'
out_dir       = '/home/lbocciol/GR1D_burn'

script_path = os.path.join(main_dir_path, "src", "hdf5_output", "hdf5_to_text.py")

# Get OMP_NUM_THREADS, or fallback to cpu_count - 1 if not set
env_threads = os.environ.get('OMP_NUM_THREADS')

if env_threads is not None:
    number_of_processes = int(env_threads)
else:
    number_of_processes = max(1, multiprocessing.cpu_count() - 1)

print(f"Using {number_of_processes} processes")

# skip some stuff
patterns_to_skip = ['backup', 'src', 'profiles']

# Cleaner folder calculation
# 1. Walk through out_dir to get all subdirectories
all_folders = [x[0] for x in os.walk(out_dir)]

# 2. Filter out folders matching skip patterns and ensure they contain HDF5 files
folders_to_convert = []
for f in all_folders:
    # Skip if any pattern is in the path
    if any(pattern in f for pattern in patterns_to_skip):
        continue

    # Check if the folder has at least one of the target HDF5 files
    if any(os.path.exists(os.path.join(f, name)) for name in ('xg.h5', 'dat.h5', 'burn.h5')):
        folders_to_convert.append(f)

dir_dict = {'main_dir_path'    : main_dir_path, \
            'script_path'      : script_path, \
            'folder_to_convert': folders_to_convert}
dir_dict['nghosts'] = 4

dir_dict['files_to_convert'] = ['xg.h5', 'dat.h5']
#dir_dict['files_to_convert'] = ['dat.h5']
#dir_dict['files_to_convert'] = []

if __name__ == '__main__':
    if not folders_to_convert:
        print("No valid folders containing xg.h5 or dat.h5 were found.")
        sys.exit(0)

    print(f"Starting conversion for {len(folders_to_convert)} folders...")
    pool = multiprocessing.Pool(number_of_processes)
    # Pass the dir_dict to the worker function via partial
    pool.map(functools.partial(run_conversions, dir_dict=dir_dict), folders_to_convert)
    pool.close()
    pool.join()
    print("Done.")

