#!/bin/bash


# load modules
module purge
module load Anaconda3
module list


# create conda environment
conda create -y -n vContact2 python=3
source activate vContact2


# install vContact2
# see https://bitbucket.org/MAVERICLab/vcontact2/src/master/
conda install -y -c bioconda vcontact2 mcl blast diamond


# deactivate conda environment
conda deactivate


# get ClusterONE
wget https://paccanarolab.org/static_content/clusterone/cluster_one-1.0.jar


# update on 12 Jul 2022
conda update -y -c bioconda vcontact2
