#check that number of lines are correct!
wc -l < "gnum.txt"


To run the script (past in the below code in the terminal containing the scripts), this is for starting calculations for many permutation matrices at once, just att numbers in the gnum.txt file for the size of the random gene list that you need a permutation matrix for:

# Count the number of lines in the file
G_FILE=gnum.txt
NUM_LINES=$(wc -l < "$G_FILE")
# Submit the SLURM job array
sbatch --array=1-$NUM_LINES "./slurm_gen_null_set.sh"
