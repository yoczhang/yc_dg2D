#/bin/bash
sed -i "s/*/.*/g" ../sysmDerivate_log.txt
sed -i "s/\^/.^/g" ../sysmDerivate_log.txt
sed -i "s/PI/pi/g" ../sysmDerivate_log.txt
sed -i "s/x,y)/x,y) /g" ../sysmDerivate_log.txt

#sed -i "s/Phi = 0/Phi = 0.*ones(size(X))/g" DxPhi.txt
#sed -i "s/Phi = 0/Phi = 0.*ones(size(X))/g" DyPhi.txt
#sed -i "s/Phi = 6/Phi = 6.*ones(size(X))/g" DyPhi.txt
