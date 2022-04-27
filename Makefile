# subdirectories
all:
	make -C ./src
	make -C ./tests

clean:
	make -C ./src clean

test:
	make -C ./tests test
