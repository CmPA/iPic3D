#!/bin/bash
i=27

w=$(printf "%05d\n" $(($i *1000)))
echo $w

#convert  \( Nongyro$w.png -trim  Vex$w.png  -trim +append  -bordercolor white -border 10x10 -background White \) \(  Jez$w.png -trim VdZ$w.png -trim +append -bordercolor white -border 10x10 -background White \) -append combinata$w.png

#convert Jez$w.png -trim VedZ$w.png -trim VidZ$w.png -trim  +append combinata_Sotto.png
convert Jez$w.png -trim Vez$w.png -trim Viz$w.png -trim  +append combinata_Sotto.png

convert Nongyro$w.png -trim  Vex$w.png  -trim +append combinata_Sopra.png
