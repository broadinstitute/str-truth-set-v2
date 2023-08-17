TAG = weisburd/hprc-pipeline

all: build push

build:
	cd ../.. && docker build -t ${TAG}:latest -f hprc_pipeline/docker/Dockerfile .

push:
	docker push ${TAG}:latest
	docker push ${TAG} | tee /dev/tty | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt
	cat sha256.txt && [ ! -z "`cat sha256.txt`" ] && sed -i.bak "s/@sha256:[^\"]*/"@`cat sha256.txt`"/"  ../*.py && rm ../*.bak
