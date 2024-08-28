# control if parameters are passed, if only three arguments are passed, the maximum number of processors is also passed
if [ $# -ne 4 ]; then
    echo $# arguments passed
    echo "Please pass the folder with the inputs and the folder where to save the outputs, the minimum number of processors and the maximum number of processors"
    exit 1
fi

inputsFolder=$1
outputFolder=$2
minProcessors=$3
maxProcessors=$4

#nodes list goes from 10000 to 100000
nodesList=$(seq 10000 10000 100000)

# iterate over all nodes in input folder
for nodes in ${nodesList[@]}; do
    # create the output folder
    outputFolderName="$outputFolder/${nodes}"
    if [ ! -d $outputFolderName ]; then
        mkdir -p $outputFolderName
    fi
    inputFolderName="$inputsFolder/${nodes}"
    # run the simulation, start from 1 processor to the maximum number of processors
    for i in $(seq $minProcessors $maxProcessors); do
        bash /home/josura/Projects/ccc/c2c-sepia/scripts/bash/cluster/performanceAnalysis/localChangingPerformanceParametersEuroHPC.sh $inputFolderName $outputFolderName $i
    done
done    

