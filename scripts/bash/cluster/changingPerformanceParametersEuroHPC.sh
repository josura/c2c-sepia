# control if parameters are passed, if only three arguments are passed, the maximum number of processors is also passed
if [ $# -ne 3 ]; then
    echo $# arguments passed
    echo "Please pass the folder with the inputs and the folder where to save the outputs and the maximum number of processors"
    exit 1
fi
