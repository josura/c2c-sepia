# control if parameters are passed, if only three arguments are passed, the maximum number of processors is also passed
if [ $# -ne 3 ]; then
    echo $# arguments passed
    echo "Please pass the folder with the inputs and the folder where to save the outputs and the maximum number of processors"
    exit 1
fi

inputsFolder=$1
outputFolder=$2
maxProcessors=$3

#nodes list goes from 1000 to 100000
nodesList=$(seq 1000 1000 100000)

# iterate over all nodes in input folder
for nodes in ${nodesList[@]}; do
    # create the output folder
    outputFolderName="$outputFolder/nodes${nodes}"
    if [ ! -d $outputFolderName ]; then
        mkdir -p $outputFolderName
    fi
    # run the simulation, start from 1 processor to the maximum number of processors
    for i in $(seq 1 $maxProcessors); do
        bash submitPerformanceExperiments.sh $inputsFolder $outputFolderName $i
    done
done