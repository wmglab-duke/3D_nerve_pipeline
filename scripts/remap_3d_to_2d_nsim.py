import sys
import os
sys.path.append(r'D:\ASCENT\fresh')
os.chdir(r'D:\ASCENT\fresh')
import pickle
import shutil

new_nsim = 1
dest_samples, dest_model, dest_sim = [3070,3150,3230], 0, 3000
source_sample, source_model, source_sim = 3000, 0, 3000
for dest_sample in dest_samples:
    with open(r'samples\{}\models\{}\sims\{}\sim.obj'.format(dest_sample, dest_model, dest_sim), 'rb') as obj:
        sim = pickle.load(obj)
        
    destnsim = r'samples\{}\models\{}\sims\{}\n_sims\{}'.format(dest_sample, dest_model, dest_sim, new_nsim)
    sourcepath = r'samples\{}\models\{}\sims\{}\n_sims\0\data\outputs'.format(source_sample, source_model, source_sim)
    
    if os.path.isdir(destnsim):
        shutil.rmtree(destnsim)
    os.makedirs(destnsim+r'\data\outputs')
    
    for file in os.listdir(sourcepath):
        if file.startswith('thresh'):
            if 'inner0' not in file:
                sys.exit('file is not inner 0')
            else:
                q = int(os.path.splitext(file.split('fiber')[-1])[0])
                l,k = sim.indices_fib_to_n(0,q)
            shutil.copyfile(sourcepath+'/thresh_inner0_fiber{}.dat'.format(q),
                            destnsim+r'\data\outputs/thresh_inner{}_fiber{}.dat'.format(l,k))
    del sim
    