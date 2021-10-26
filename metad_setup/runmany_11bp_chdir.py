import subprocess
import sys
import os
import shutil, glob
import pickle as pkl

idx = int(sys.argv[1])-1
temp_diff = int(sys.argv[2])
pace = 500

seq_dict = {'CCTATATATCC':[327, 315, 308, 307],
            'CGCATATATAT':[330, 312, 310, 312],
            'TTTTTTTTTTT':[338, 334, 325, 320],
            'TATAGCGATAT':[334, 329, 321, 303]
            }

seq_list = list(seq_dict.keys())
base_list = [0, 2, 4, 6]

inverse = False
if inverse: inv='_inv'
else: inv=''

# update plumed file with these values
sigma = 0.01
height = 0.6
bias_factor = 4
traj_save = 10000 
run_steps = 1000000000

## pair each base and temp 1-1
base_idx, seq_idx = idx%len(base_list), idx//len(base_list)
base, seq = base_list[base_idx], seq_list[seq_idx] 
temp = seq_dict[seq][base_idx] + temp_diff

# specify relevant file names:
main_file = '11bp_abasic_plumed.in'
plumed_file = f'11bp_abasic-{base}_plumed{inv}.dat'
conf_file = f'11bps/{seq}/conf_lammps_base-{base}_strand-1.in'

# account for control
if base == 0: 
    plumed_file = f'11bp_plumed{inv}.dat'
    conf_file = f'11bps/{seq}/conf_lammps.in'

# make an output directory
#dir_name = f'./11bp_runs{inv}/{seq}_base-{base}_temp-{temp}'
dir_name = f'./11bp_varyT{inv}_s-{sigma}_bf-{bias_factor}_ns-{run_steps:.0e}/'
subdir_name = f'{dir_name}/{seq}_base-{base}_temp-{temp}'

subprocess.call(['mkdir', dir_name])
subprocess.call(['mkdir', subdir_name])

# copy lammps and plumed files in to dir:
subprocess.call(['cp', main_file, plumed_file, subdir_name])

# change temp in main copy
with open(f'{subdir_name}/{main_file}', 'r') as file:
    lines = file.readlines()
    for n, line in enumerate(lines):
        if 'variable T equal' in line:
            lines[n] = f'variable T equal {temp}\n'
        elif 'fix pl all' in line:
            lines[n] = f'fix pl all plumed plumedfile {plumed_file} outfile p.log\n'
        elif 'dump mydump' in line:
            lines[n] = f'dump mydump all custom {traj_save} traj.xyz id type x y z\n'
        elif 'read_data' in line:
            lines[n] = f'read_data ../../{conf_file}\n'
        elif 'run ' in line:
            lines[n] = f'run {run_steps}\n'
    
    file.close()
with open(f'{subdir_name}/{main_file}', 'w') as file: 
    file.writelines(lines)
    file.close()
                                                            
# change vars in plumed copy
with open(f'{subdir_name}/{plumed_file}', 'r') as file:
    lines = file.readlines()
    for n, line in enumerate(lines):
        if 'FILE=' in line:
            if 'PRINT' in line:
                lines[n] = line.replace('FILE=COLVAR', f'FILE=COLVAR')
            else:
                lines[n] = f'FILE=HILLS\n'
        elif 'TEMP=' in line:
            lines[n] = f'TEMP={temp}\n'
        elif 'BIASFACTOR=' in line:
            lines[n] = f'BIASFACTOR={bias_factor}\n'
        elif 'HEIGHT=' in line:
            lines[n] = f'HEIGHT={height}\n'
        elif 'SIGMA=' in line:
            lines[n] = f'SIGMA={sigma}\n'
        elif 'PACE=' in line:
            lines[n] = f'PACE={pace}\n'
            
    file.close()
with open(f'{subdir_name}/{plumed_file}', 'w') as file: 
    file.writelines(lines)
    file.close()

os.chdir(subdir_name)
subprocess.call(['mpirun', '-np', '1', '/scratch/midway3/mikejones/Mike/lammps-7Aug19/src/lmp_mpi', '-in', f'{main_file}'])
