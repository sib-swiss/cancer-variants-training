#!/usr/bin/env bash

cd "$HOME"/Documents/cancer-variants-training/Docker/vscode
podman build --format docker -t thefericcio/cancer-variants-vscode:latest .
podman push thefericcio/cancer-variants-vscode:latest