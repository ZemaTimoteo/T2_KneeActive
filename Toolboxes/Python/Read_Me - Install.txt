Install Instructions

Install Anaconda
	- Read and download from the following folder
	https://docs.anaconda.com/free/anaconda/install/index.html

Create python 3.7
	-On command line type:
	conda create --name qMRI  python=3.7

Activate environment
	-On command line type:
	conda activate qMRI

Create libraries file
	req.txt file

Install Libraries
	conda install --file req.txt
	OR on pycharm command line (change directory to the folder where 'req.txt' is saved:
	pip install -r req.txt