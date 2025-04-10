# For accuracy look at line 36 of WriteCard.py
path=$(pwd)
ntot=1
nbin=30
Gene="Gen.txt"
part="ttx"
comm="generate p p > t t~ [QCD]"
name=$part"prod_NLO"
s=1300
	
echo "import model sm-no_b_mass" >> $Gene
echo $comm >> $Gene
echo "output "$name >> $Gene
echo "y" >> $Gene
echo "exit" >> $Gene
    
cd ~/mg5amcnlo
./bin/mg5_aMC $path/$Gene
cd $path
rm $Gene
rm dataTot/MGXsec_$part.txt
rm dataTot/dXsec_$name*
rm data/MGXsec_$part.txt
rm data/dXsec_$name*
    
cp analysis_HwU_pp_ttx_v2.f ~/mg5amcnlo/$name/FixedOrderAnalysis/.
cp setscales.f ~/mg5amcnlo/$name/SubProcesses/
#cp cuts.f ~/mg5amcnlo/$name/SubProcesses/.
#cp pdg2pdf_lhapdf6.f ~/mg5amcnlo/$name/Source/PDF/.
#cp BinothLHA.f ~/mg5amcnlo/$name/SubProcesses/. 

for j in $(seq 1 $ntot)
do
    ij=$(($j-1))
    lame="Launch"
    lame+=$name
    cmd="launch "$name" NLO"
    echo $cmd > $lame
    echo "set WZ 0" >> $lame
    echo "set WW 0" >> $lame
    echo "set WT 0" >> $lame
    echo "set WH 0" >> $lame
    echo "set MZ 91.1876" >> $lame
    #echo "set Mt 172.7" >> $lame
    echo "set Mt 173.2" >> $lame
    echo "set Gf 1.1663787e-05" >> $lame
    
    python WriteCard.py $name $j $ntot $nbin $s
    
    cd ~/mg5amcnlo
    mv $name/Cards/FO_analyse_card2.dat $name/Cards/FO_analyse_card.dat
    ./bin/mg5_aMC $path/$lame
    cd $name
    if [ $j -lt 10 ]
    then
	cp Events/run_0$j/MADatNLO.HwU $path/dataTot/dXsec_$name$j.txt
    else
	cp Events/run_$j/MADatNLO.HwU $path/dataTot/dXsec_$name$j.txt
    fi
    cd $path
    python Readhist.py $name $j		
    rm $lame
done

