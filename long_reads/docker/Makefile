TAG = weisburd/long-reads

all: build push

build:
	cd ../.. && docker build -t ${TAG}:latest -f long_reads/docker/Dockerfile .

push:
	docker push ${TAG}:latest
	docker push ${TAG} | tee /dev/tty | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt
	cat sha256.txt && [ ! -z "`cat sha256.txt`" ] && sed -i.bak "s/long-reads@sha256:[^\"]*/"long-reads@`cat sha256.txt`"/"  ../*.py && rm ../*.bak
