### Build the Docker container
sudo docker build -t kauralasoo/txrevise .

### Push to DockerHub
docker push kauralasoo/txrevise

### Build a local copy of the Singularity container
singularity build txrevise.img docker://kauralasoo/txrevise:latest
