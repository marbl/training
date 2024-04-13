docker build --platform linux/amd64  -t assembly-102 .
docker save --output assembly-102.tar assembly-102:latest
gzip assembly-102.tar
