ncores=10
# To Run : bash RunPar.sh
for i in $(seq 1 $ncores)
do    
    let j="$i"-1
    dir="Part_$j"
    if [ -d "${dir}" ];
    then
    	rm -rf $dir
    fi
    mkdir $dir
    cd $dir
    ln -s ../RUN ./
    #nohup ./RUN 14400 $i $ncores &# > output.log &
    nohup ./RUN 21100 $i $ncores &# > output.log &
    #nohup ./RUN 21200 $i $ncores &# > output.log &
    echo "PID of process n°" $i " = " $! 
    cd ../
done
# Wait for all background jobs to finish before printing the completion message
wait

# Once all jobs have finished
echo "All computations are done!"