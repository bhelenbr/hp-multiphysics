# Use the Bash Shell
#$ -S /bin/bash

# Put a unique job name here
# so that you can monitor it in the queue
#$ -N unsteady

# Notifying user by email at the beginning,end,abort,suspensions of the job run
#$ -M helenbrk@clarkson.edu
#$ -m eas

# Tell GE to run the job from the current working directory
#$ -cwd 

# Uncomment to pass all current environment variables
#$ -V
# Uncomment to pass a single environment variable
# #$ -v VAR_TO_PASS

# Redirecting standard output / error to files
# named "output" and "errors"
#$ -o output
#$ -e errors

# The max walltime for this job is 1 hour 31 minutes
# Set this longer than required for run!!!
#$ -l h_rt=240:00:00

local_dir=$(pwd | cut -c 23-)
echo "This job was started in ${local_dir}"
echo "The start time was $(date)"
echo "The job id is $JOB_ID"

# Executing the program
./convergence.command
