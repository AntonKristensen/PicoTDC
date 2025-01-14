cd /home/anton/Applications/TOPAS/OpenTOPAS-build
cmake ../OpenTOPAS -DCMAKE_INSTALL_PREFIX=../OpenTOPAS-install -DTOPAS_EXTENSIONS_DIR=/home/anton/Desktop/gitrepos/CERN/Simulations/topasextensions
sudo make -j4
sudo mv topas ../OpenTOPAS-install/bin/topas
