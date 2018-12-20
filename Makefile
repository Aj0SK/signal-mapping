installminimap:
	git clone https://github.com/lh3/minimap2 && cd minimap2 && make

testminimap:
	./minimap2/minimap2 -a test/MT-human.fa test/MT-orang.fa > testout/test.sam
