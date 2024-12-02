
cd /Users/geertvangeest/Documents/repositories/cancer-variants-training
docker run --rm -v $PWD:/config/project -p 8443:8443 geertvangeest/cancer-variants-vscode:latest
