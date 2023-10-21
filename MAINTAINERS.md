## Docker images

The docker images for `dnmtools` are hosted in GitHub Container registry.  The
process of building and pushing the image to the registry is handled by the
workflow specified in
[docker-build.yml](https://github.com/smithlabcode/dnmtools/blob/master/.github/workflows/docker-build.yml).
The build instruction is in
[Dockerfile](https://github.com/smithlabcode/dnmtools/blob/master/Dockerfile).
You can see the published images
[here](https://github.com/smithlabcode/dnmtools/pkgs/container/dnmtools).

The workflow is triggered either manually or automatically by a tag event of
type `v*.*.*`, which is intended for new releases. Currently, publishing the
images can happen only to commits tagged by a version number. This is intended
to associate every docker image with a version number. This means that there is
no option to push the image for the latest commit if it is not tagged by
a version number.

### Automatic build and publish in a tag event

In a tag event of type `v*.*.*`, such as new release or retagging of versoin
number, this work flow is triggered to build and publish the image for the
tagged version number. The published image is tagged with SHA hash and the
version number.  It is also taged with `latest` if the version number is the
latest.

### Manual build (and publish)

Manual trigger is intedned to test the image build processes as well as publish
an image for an existing version.  In
[Actions](https://github.com/smithlabcode/dnmtools/actions), go to `Docker image
build` under `All workflows` and click `Run workflow` and choose from the
following options:

1. `Build latest commit`: for testing for the latest commit
2. `Build existing version`: for testing a particular version
3. `Build + push existing version`: for publishing a particular version  

For options 2 and 3, specify the version number in the form `v*.*.*`. If not
specified, the workflow will assume the latest verion.

### Use scenarios 

**Before a new release**: It is a good idea to test image building before a new
release. Manually trigger the workflow with opiton 1. If it builds with no
issues, make a new release and the image will automatically be built and
published. 

**Publish an existing version**: It is possible to publish a docker image for an
existing version by option 3 in the manual trigger. First, test build using
option 2, and then publish using option 3.  The published image is tagged with
SHA hash and the version number.  It is also taged with `latest` if the version
number is the latest. If option 3 is deployed with a version number for which
a docker image already exists, it will simply rebuild and update the existing
image.

**Deleting an image**: If you have owner access to `smithlabcode`, you can
delete an image by going
[here](https://github.com/smithlabcode/dnmtools/pkgs/container/dnmtools/versions)
and manually delete a version.



## Installation
The image can be pulled by one of the following commands.

```bash
docker pull ghcr.io/smithlabcode/dnmtools:latest
docker pull ghcr.io/smithlabcode/dnmtools:[7-DIGIT SHA]
docker pull ghcr.io/smithlabcode/dnmtools:v[VERSION NUMBER] #(e.g. v1.4.2)
```

