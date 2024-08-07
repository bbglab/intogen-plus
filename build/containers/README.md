# Containers

This folder contains all the containers needed by the intogen pipelines specifically.

Every container will have an independent folder here or a link to another folder in the repository. Every of those folders need to include a `Dockerfile` with the specs to build the Docker image.

## CI/CD

The Github workflow `containers.yaml` will run the following jobs:
- Run some quality checks to make sure that we follow the best practices
- Determine which containers need to be built and/or pushed:
    - When the push corresponds with a tag, all the containers will be built and pushed.
    - When this is a regular push into a branch, only the containers that were modified will be built to make sure that they are not broken by the changes.

You can run the workflow locally with the following command:

```shell
act -W '.github/workflows/containers.yaml'
```

Note that you will need to install the tool [act](https://github.com/nektos/act).

## Local Development

You can run the following commands (note that you can skip the `-C build/containers` argument if you run the commands from the `build/containers` folder):

### Running quality checks

We run two types of checks:
- [shellcheck](https://github.com/koalaman/shellcheck) to make sure that bash scripts follow the best practices
- [hadolint](https://github.com/hadolint/hadolint) to make sure that Dockerfiles follow the best practices

You will need to install them in your computer before running the following commands:

- To run both checks:

    ```shell
    make -C build/containers checks
    ```

- Or individually with:

    ```shell
    make -C build/containers shellcheck
    ```

    ```shell
    make -C build/containers hadolint
    ```

### Building all the containers

```shell
make -C build/containers build-all
```

### Building only the containers that were modified

```shell
make -C build/containers build-changed
```
