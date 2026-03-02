IMAGE   ?= needletail-server
TAG     ?= latest
PORT    ?= 8002

.PHONY: docker docker-run docker-push clean

## Build the Docker image
docker:
	docker build -t $(IMAGE):$(TAG) .

## Run the container (interactive, auto-remove)
docker-run: docker
	docker run --rm -it -p $(PORT):8002 $(IMAGE):$(TAG)

## Push to a registry (set IMAGE=registry/repo first)
docker-push: docker
	docker push $(IMAGE):$(TAG)

## Cargo release build (no Docker)
build:
	cargo build --release -p needletail-server

## Run all core tests
test:
	cargo test -p needletail-core

## Clean build artifacts
clean:
	cargo clean
