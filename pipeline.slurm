#!/bin/bash
#SBATCH -J v111_ascent
#SBATCH -p wmglab
#SBATCH -o output.txt
#SBATCH -e error.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=384:00:00
#SBATCH --mem=200G

python3 run pipeline -w 10 606 607 608 609 610 611 612 613 614 615 616 617 618 619 620 621 622 623 624 625 626 627 628 629 630 631 632 633 605 "$@"
