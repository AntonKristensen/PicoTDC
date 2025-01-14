# need to run as root
echo "Starting proceedure...."

cd openocd
src/openocd -s tcl -f interface/cmsis-dap.cfg -f target/rp2040.cfg -c "init;halt;" &




# Force CPU to 120MHz (40MHz XTAL)
#(gdb) set *((unsigned int)0x40028000) = 0x00000001 # Sets refdiv=1
#(gdb) set *((unsigned int)0x40028004) = 0x00000004 # Powers down PLL
#(gdb) set *((unsigned int)0x40028008) = 0x00000024 # Sets feedback divisor to 36
#(gdb) set *((unsigned int)0x4002800C) = 0x00062000 # Sets the two divisors to 6 and 2

# Force USB to 48MHz (40MHz XTAL)
#(gdb) set *((unsigned int)0x4002c000) = 0x00000001 # Sets refdiv=1
#(gdb) set *((unsigned int)0x4002c004) = 0x00000004 # Powers down PLL
#(gdb) set *((unsigned int)0x4002c008) = 0x00000024 # Sets feedback divisor to 36
#(gdb) set *((unsigned int)0x4002c00C) = 0x00065000 # Sets the two divisors to 6 and 5



gdb-multiarch -ex "target extended-remote localhost:3333" \
	-ex "set *((unsigned int)0x40028000) = 0x00000001" \
	-ex "set *((unsigned int)0x40028004) = 0x00000004" \
	-ex "set *((unsigned int)0x40028008) = 0x00000024" \
	-ex "set *((unsigned int)0x4002800C) = 0x00062000" \
	-ex "x 0x40028000" \
	-ex "x 0x40028004" \
	-ex "x 0x40028008" \
	-ex "x 0x4002800c" \
	-ex "set *((unsigned int)0x4002c000) = 0x00000001" \
	-ex "set *((unsigned int)0x4002c004) = 0x00000004" \
	-ex "set *((unsigned int)0x4002c008) = 0x00000069" \
	-ex "set *((unsigned int)0x4002c00C) = 0x00065000" \
	-ex "x 0x4002c000" \
	-ex "x 0x4002c004" \
	-ex "x 0x4002c008" \
	-ex "x 0x4002c00c" \
	-ex "cont" &

sleep 1
killall gdb-multiarch
killall openocd

echo "Done"
