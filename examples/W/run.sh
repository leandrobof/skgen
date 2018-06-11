sed -i 's/Udd  Upp  Uss/0.2  0.2  0.2/' W-W.skf
sed -i 's/-0.0473147/0.040/' W-W.skf

cp bulk dftb_in.hsd
dftb+
dp_bands band.out w_band
efermi=`grep "Fermi energy" detailed.out | awk '{print $5}' `
python plot2.py $efermi
