if [ $# -eq 2 ]
then

rm MAP_BP.jpg
rm CONTOUR_BP.jpg
echo "set term jpeg" >dessin.gnu
echo "set palette functions (sin(gray*2*pi)+1)/2, (sin(gray*2*pi + 2*pi/3)+1)/2,(sin(gray*2*pi + 4*pi/3)+1)/2" >>dessin.gnu
echo "unset key" >>dessin.gnu
echo "set output 'MAP_BP.jpg'" >>dessin.gnu
echo "set cbrange [0.0:7.5]" >>dessin.gnu
echo "set xrange [0:$1]" >>dessin.gnu
echo "set yrange [0:$2]" >>dessin.gnu
echo "set pm3d map" >>dessin.gnu
echo "splot 'TT.bin' binary array=$1x$2" >>dessin.gnu
echo "reset" >>dessin.gnu
echo "unset key" >>dessin.gnu
echo "set term jpeg" >>dessin.gnu
echo "set output 'CONTOUR_BP.jpg'" >>dessin.gnu
echo "set contour base" >>dessin.gnu
echo "set cntrparam level incremental 0.,0.3,7." >>dessin.gnu
echo "set xrange [0:$1]" >>dessin.gnu
echo "set yrange [0:$2]" >>dessin.gnu
echo "set sample 5,5" >>dessin.gnu
echo "set isosample 5,5" >>dessin.gnu
echo "unset surface" >>dessin.gnu
echo "set table 'cont.dat'" >>dessin.gnu
echo "splot 'TT.bin' binary array=$1x$2" >>dessin.gnu
echo "unset table" >>dessin.gnu
echo "plot 'cont.dat' with lines lt -1 lw 1.5" >>dessin.gnu

gnuplot dessin.gnu
rm cont.dat
#eog MAP.jpg &
#eog CONTOUR.jpg

#string=`echo $1 $2 | awk '{print "splot TT.bin binary array="$1"x"$2}'`
#set palette defined (0   "#086f00",  0.1 "#4fb847", 0.2 "#C5c471",0.5 "#cf3c08", 1 "#ffffff")
#set palette functions (sin(gray*2*pi)+1)/2, (sin(gray*2*pi + 2*pi/3)+1)/2,(sin(gray*2*pi + 4*pi/3)+1)/2

else

##########################
# warning messages here 
##########################

if [ $# -eq 1 ]
then
if [ $1 = -h ]
then
echo @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
echo please enter number of points along x and number of points along z
echo BP example 126 626
echo @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
exit 0
fi
fi

#########################
# missing arguments
#########################

echo we need two integers ...

fi


