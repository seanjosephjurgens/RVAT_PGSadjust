# Script to install the SAIGE docker for analyses on UKB RAP
#### See: https://saigegit.github.io/SAIGE-doc/docs/UK_Biobank_WES_analysis.html 

docker pull wzhou88/saige:1.0.9
docker save -o saige_1.0.9.tar.gz wzhou88/saige:1.0.9

# make it readable for other users in your project
chmod a+r saige_1.0.9.tar.gz

# store it in your image folder on DNAnexus
dx mkdir docker_images
dx cd docker_images/
dx upload saige_1.0.9.tar.gz --destination path/to/docker/
