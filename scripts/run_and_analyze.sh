#!/usr/bin/bash

# assume you are running from the scripts directory!

declare ROOT_DIR="/home/seang/Dev/Git/UT_Austin/ConnectivityAlg"
declare SCRIPTS_DIR="$ROOT_DIR/scripts/"
declare TESTS_DIR="$ROOT_DIR/tests"
declare BUILD_DIR="$ROOT_DIR/build"
declare TESTS_OUT_DIR="$BUILD_DIR/tests"
declare DATA_DIR="$ROOT_DIR/data"
declare DATA_OUT_DIR="$DATA_DIR/outputs"
declare DATA_FILE="$DATA_OUT_DIR/con_decay.bin"

cd $TESTS_DIR

if make; then 
	cd $TESTS_OUT_DIR
	./connectionTest
	cd $SCRIPTS_DIR
	source venv/bin/activate
	# ./analyze_connectivity.py $DATA_FILE
	# deactivate
fi

