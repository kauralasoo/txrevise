### Build the Docker container
#sudo docker build -t andreasvija/txrevise .

### Push to DockerHub
#docker push andreasvija/txrevise

### Build a local copy of the Singularity container
cd scripts
singularity build txrevise.img docker://andreasvija/txrevise:latest
