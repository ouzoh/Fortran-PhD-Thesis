 #!/bin/bash
 ##Author    : Ogaga G Uzoh
 ##University: University College London
 ##Date      : 26-02-2013
 ##            Extract geometries
 ##Update    : 03-03-2013
 ##            Extract torsion, E(intra) and RMS gradient
 ##            also states the explicit atoms involved in torsion
 ##            Hence we know the bond associated with the flexible torsion
################################################################################
 logfile="./*.log"
################### NOTHING SHOULD BE CHANGED BELOW THIS LINE ##################

babel="/usr/bin/babel"
mkdir -p mol2_dir
cd mol2_dir
cp ../$logfile .

if [ ! -f $logfile ]; then
  echo "ERROR: res file $logfile does not exist. Exiting."
  echo "     : check file.in if $logfile present"
  exit
fi

echo "###########################"
echo "all input files are present"
echo "###########################"

#############
filename=$(basename "$logfile")
extension="${filename##*.}"
filename="${filename%.*}"

## Extract the optimised geometries #####################################
sed '1,/ Symbolic Z-matrix:/d;/Variables:/,$d' ../$filename.log > temp
sed '1d' temp > top_zmat

#extract the torsion been scanned from the *.log file
awk  '/       Variables:/,/ GradGrad/' $filename.log  > temp


awk  '/Scan/' temp > temp_d
##awk "NR==1" temp_d > temp
##awk  '/Scan/' temp > temp_d
##{ if (NR<11) print $2 }'
dof=$(awk '{print $1}' temp_d)
echo $dof
#extract the explict atoms associated with torsion
a1=$(awk "/$dof/ {print FNR}" top_zmat)
echo $a1
grep "$dof" top_zmat > temp_a
a2=$(awk  '{print $2}' temp_a)
a3=$(awk  '{print $4}' temp_a)
a4=$(awk  '{print $6}' temp_a)
echo "$dof atoms are:"
echo "$a1 $a2 $a3 $a4" 
rm -f top_a

nmol=`wc -l < top_zmat`
errtop=7
tvar=$[(${nmol}-1)+(${nmol}-2)+(${nmol}-3)]
errbottom=2
tval=$[$tvar+$errtop+$errbottom]
echo "variables:" > vari
echo "" >line
echo "constants:" > const

## extract zmatrix
awk '/-- Stationary point found./,/ GradGradGrad/' ../$filename.log > temp
echo "extract gzmat file start for"

##loop starts here
for j in {1..42}
  do
    echo "$j"
    sed -n "8,${tval}p" temp > temp1
    cp temp tempp
    sed '$d' temp1 > temp2
    sed '$d' temp2 > temp1
    awk '{print $2,"=",$3}' temp1 > bottom_zmat1
    grep $dof  bottom_zmat1 > tempdof

    sed "/$dof/d" bottom_zmat1 > bottom_zmat
    awk 'NR > v1'  v1="${tval}" tempp > temp 
    echo '%mem=1800MB'         >ttop_zmat
    echo '%nprocshared= 1'    >>ttop_zmat
    echo '%nosave'            >>ttop_zmat
    echo '%chk=obs'           >>ttop_zmat
    #echo '#MP2/6-31+G(d) POpt' >>ttop_zmat
    #echo '#PBE1PBE/6-31+G(d) density=current FormCheck=All' >> ttop_zmat
    echo '#MP2/6-31+G(d) density=current FormCheck=All' >> ttop_zmat
    #!constrained optimisation
    #!Extracts as gassian input
    #echo '#PBE1PBE/6-31+G(d) Opt=Z-matrix' >> ttop_zmat
    echo '#Pop=ChelpG nosymm' >>ttop_zmat
    echo ''                   >>ttop_zmat
    ddof=$(grep "$dof" temp1|awk '{printf("%d",$3)}')
    #echo "PBE0 Opt zmat${ddof}">>ttop_zmat
    echo "MP2 Opt zmat${ddof}">>ttop_zmat
    echo ''                   >>ttop_zmat
    echo '0 1'                >>ttop_zmat
    ################
    #!modify
    #include constant: for constrained optimisation
    cat ttop_zmat top_zmat vari bottom_zmat const tempdof line  > zmat${ddof}.com
    #exclude constant "standard"
    #cat ttop_zmat top_zmat vari bottom_zmat tempdof line  > zmat${ddof}.com

    rm -f ttop_zmat bottom_zmat*

    rm -f temp? tempdof 
    $babel -igzmat zmat${ddof}.com -omol2 zmat${ddof}.mol2
    #$babel -igzmat zmat${ddof}.com -oreport zmat${ddof}.report
    #$babel -igzmat zmat${ddof}.com -oxyz zmat${ddof}.xyz
    #rm -f zmat${ddof}.com
    #awk 'NR > 2 {print}' zmat${ddof}.xyz > zmat${ddof}.cxyz
    rm -f zmat${ddof}.com
    #cat ttop_zmat zmat${ddof}.cxyz bottom> zmat${ddof}S.com
    rm -f ttop_zmat zmat${ddof}.cxyz zmat${ddof}.xyz bottom 
  done

rm -f top_zmat line vari const 

echo "extract gzmat files done!"

####extract torsion been scanned and at each optimised conformation
##################ab initio energy, torsion and its rms derivatives

flag=EUMP2
if grep -Fxq "$flag" $filename.log
 #if calculation is a MP2 
then
    awk  '/ EUMP2:/,/Optimization completed./'  $filename.log  > temp
    grep ' EUMP2:\|Optimization completed.'   temp > temp_t
    grep -B1 "Optimization completed." temp_t > temp_1
    awk  '/ EUMP2:/'  temp_1 > temp_2
    awk  '{print $6}' temp_2 > temp_E
  #if calculation is not a MP2
else
    awk  '/ SCF Done:/,/Optimization completed./' $filename.log  > temp
    grep ' SCF Done:\|Optimization completed.' temp > temp_t
    grep -B1 "Optimization completed." temp_t > temp_1
    awk  '/SCF Done:/' temp_1 > temp_2
    awk  '{print $5}'  temp_2 > temp_E
fi    
    
    #extract rms gradient
    awk  '/Cartesian Forces:/,/Optimization completed./' $filename.log  > temp
    grep 'Cartesian Forces:\|Optimization completed.' temp > temp_t
    grep -B1 "Optimization completed." temp_t > temp_1
    awk  '/Cartesian Forces:/' temp_1 > temp_2
    awk  '{print $6}' temp_2 > temp_G
    
    #extract degree of freedom
    awk '/Optimization completed./,/ SCF Done:/' $filename.log  > temp
    awk "/\!      ${dof}/" temp > temp_d
    awk '{print $3}' temp_d > temp_d1
    
    #paste results
    paste temp_d1 temp_E          temp_G > data_E
    rm temp temp_*

#evaluate the intramolecular energy
min=$(sort -nk 2 data_E | awk 'NR==1 {print $2}')
echo 'the minimum is'" $min"
awk '{print '"-1*$min"'+$2}' data_E > temp
#converting Hatree/mol to kJ/mol
awk '{print $1*2625.49962}' temp > temp_eintra
#coverting Hatree/bohr to KJ/Angstrom
awk '{print $3*2625.49962/0.529177249}' data_E > temp_grad_conv
paste data_E temp_eintra temp_grad_conv > temp_out
echo "$nmol" > summary.r
echo "X(deg)  E(H)            Eintra(KJ/mol)  dE/dX(H/Bohr)      dE/dX(KJ/mol/A)" >> summary.r
awk '{printf "%-6.1f %15.9f %9.6f %19.9f %15.6f\n", $1,$2,$4,$3,$5}' temp_out >> summary.r
rm -f temp temp_* 
rm -f $logfile
cp summary.r ../
rm -f summary.r data_E zmat.mol2
