#!/bin/bash

if [[ $OSTYPE == 'linux-gnu'* ]] || [[ $OSTYPE == 'darwin'* ]];
then  

	echo $OSTYPE	
	
	echo "Preparing a new conda enviroment called PCN"
	conda create -n PCN python=3.10
	conda install -n PCN pip
	pip install configparser

	echo "Preparing config file"
	cd pcn
	cd tools
	python create_config_file.py
	
	echo "Installing libraries"
	conda activate PCN
	pip install scipy
	pip install regex
	pip install numpy
	pip install fuzzy-c-means
	pip install git+https://github.com/palash1992/GEM.git
	pip install cdlib
	conda install -c conda-forge -c schrodinger pymol-bundle
	pip install pytz
	pip install python-dateutil
	pip install node2vec==0.4.4
	pip install pyqt5
	pip install leidenalg
	pip install numpy==1.21.6
fi
