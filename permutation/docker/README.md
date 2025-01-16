## Docker image for running convergence permutation tests in Terra

Copyright (c) 2024-Present, [Ryan L. Collins](mailto:Ryan_Collins@dfci.harvard.edu) and the Dana-Farber Cancer Institute.  
Distributed under terms of the [GNU GPL v2.0 License](/LICENSE) (see `LICENSE`).  

---  

This subdirectory contains the Docker build information for a stand-alone Docker image used for permutation testing.  

This image is available on Dockerhub here:
`vanallenlab/gsc-perm:latest`  

---  

Build instructions (not necessary for most users):  
```
export TAG=my_build_tag
docker build \
	--platform linux/x86_64 \
	--progress plain \
	--tag vanallenlab/gsc-perm:$TAG \
	.
docker push vanallenlab/gsc-perm:$TAG
```