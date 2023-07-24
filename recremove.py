# go through each directory in samples
import os
import shutil

for dire in os.listdir("samples"):
    checkdir = os.path.join("samples", dire, "models", "0", "sims", "10")
    if os.path.isdir(checkdir):
        print("Removing " + checkdir)
        shutil.rmtree(checkdir)
        print("Removed " + checkdir)
