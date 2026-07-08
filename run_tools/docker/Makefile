TAG = weisburd/run-tools

all: build push

build:
	cd ../.. && docker build -t ${TAG}:latest -f run_tools/docker/Dockerfile .

push:
	docker push ${TAG}:latest
	docker push ${TAG} | tee /dev/tty | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt
	cat sha256.txt && [ ! -z "`cat sha256.txt`" ] && sed -i.bak "s/run-tools@sha256:[^\"]*/"run-tools@`cat sha256.txt`"/"  ../*.py && rm ../*.bak
