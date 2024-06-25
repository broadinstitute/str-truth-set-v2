TAG = weisburd/bwa

all: build push

build:
	docker build -t ${TAG}:latest .

push:
	docker push ${TAG}:latest
	docker push ${TAG} | tee /dev/tty | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt
	cat sha256.txt && [ ! -z "`cat sha256.txt`" ] && sed -i.bak "s/bwa@sha256:[^\"]*/"bwa@`cat sha256.txt`"/"  ../*.py && rm ../*.bak
