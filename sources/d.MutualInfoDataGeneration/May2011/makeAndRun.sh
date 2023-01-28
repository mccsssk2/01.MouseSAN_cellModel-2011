# make directorries and run the san binary
cc -o san Oct29MouseSANPNEW.c -lm
mkdir /net/atrium/space/sanjay1/SKcprogs_2004to2010/mouse/Oct2010_SAN/Peripheral_ALL/Run1
mkdir /net/atrium/space/sanjay1/SKcprogs_2004to2010/mouse/Oct2010_SAN/Peripheral_ALL/Run2
mkdir /net/atrium/space/sanjay1/SKcprogs_2004to2010/mouse/Oct2010_SAN/Peripheral_ALL/Run3
mkdir /net/atrium/space/sanjay1/SKcprogs_2004to2010/mouse/Oct2010_SAN/Peripheral_ALL/Run4
mkdir /net/atrium/space/sanjay1/SKcprogs_2004to2010/mouse/Oct2010_SAN/Peripheral_ALL/Run5
mkdir /net/atrium/space/sanjay1/SKcprogs_2004to2010/mouse/Oct2010_SAN/Peripheral_ALL/Run6
#
cp san Peripheral_ALL/Run1/ 
cp san Peripheral_ALL/Run2/ 
cp san Peripheral_ALL/Run3/ 
cp san Peripheral_ALL/Run4/ 
cp san Peripheral_ALL/Run5/ 
cp san Peripheral_ALL/Run6/ 