#!/bin/bash
for i in {0..47}
do
z=$(($i * 1000))
#echo $z

w=$(printf "%05d\n" $(($i *1000)))
echo $w
convert dBdtX$z.png -trim Epar$z.png  -trim  T*par$z.png -trim  T*perp2$z.png -trim +append  combo$w.png
done

convert combo*.png combo.gif