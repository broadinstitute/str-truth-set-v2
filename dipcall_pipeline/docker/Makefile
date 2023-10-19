TAG = weisburd/create-hprc-vcfs

all: build push

build:
	cd ../.. && docker build -t ${TAG}:latest -f create_hprc_vcfs/docker/Dockerfile .

push:
	docker push ${TAG}:latest
	docker push ${TAG} | tee /dev/tty | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt
	cat sha256.txt && [ ! -z "`cat sha256.txt`" ] && sed -i.bak "s/create-hprc-vcfs@sha256:[^\"]*/"create-hprc-vcfs@`cat sha256.txt`"/"  ../*.py && rm ../*.bak
