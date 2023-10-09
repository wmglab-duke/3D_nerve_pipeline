import os
import sys

sys.path.append(r'D:\ASCENT\fresh')
os.chdir(r'D:\ASCENT\fresh')
import pickle
import shutil

new_nsim = 1
dest_samples, dest_model, dest_sim = [3070, 3150, 3230], 0, 3000
source_sample, source_model, source_sim = 3000, 0, 3000
for dest_sample in dest_samples:
    with open(fr'samples\{dest_sample}\models\{dest_model}\sims\{dest_sim}\sim.obj', 'rb') as obj:
        sim = pickle.load(obj)

    destnsim = fr'samples\{dest_sample}\models\{dest_model}\sims\{dest_sim}\n_sims\{new_nsim}'
    sourcepath = fr'samples\{source_sample}\models\{source_model}\sims\{source_sim}\n_sims\0\data\outputs'

    if os.path.isdir(destnsim):
        shutil.rmtree(destnsim)
    os.makedirs(destnsim + r'\data\outputs')

    for file in os.listdir(sourcepath):
        if file.startswith('thresh'):
            if 'inner0' not in file:
                sys.exit('file is not inner 0')
            else:
                q = int(os.path.splitext(file.split('fiber')[-1])[0])
                l, k = sim.indices_fib_to_n(0, q)
            shutil.copyfile(
                sourcepath + f'/thresh_inner0_fiber{q}.dat',
                destnsim + fr'\data\outputs/thresh_inner{l}_fiber{k}.dat',
            )
    del sim
