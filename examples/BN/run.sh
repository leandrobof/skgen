sed -i 's/Udd  Upp  Uss/0.2246  0.1510  0.2246/' N-N.skf
sed -i 's/Udd  Upp  Uss/0.4954  0.523305  0.4954/' B-B.skf

cp bulk dftb_in.hsd
../dftb+
dp_bands band.out BN_band
efermi=`grep "Fermi energy" detailed.out | awk '{print $5}' `
python plot2.py -1.5970
