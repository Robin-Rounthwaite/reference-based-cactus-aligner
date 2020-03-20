build:
	# Build custom docker image
	echo $(USER)
	docker build -f Dockerfile -t $(USER)-cactus-minimap-aligner .

debug:
	# Run our image with the local version of the app and shell into it (remove --user to have sudo power?)
	docker run -it --rm \
		-v `pwd`:/app \
		--entrypoint /bin/bash \
		$(USER)-cactus-minimap-aligner
		# --user=`id -u`:`id -g` (goes above --entrypoint, with a backslash at end)

run:
	# Run the image in a container locally
	docker run -it --rm \
		$(USER)-cactus-minimap-aligner
		# $(USER)-cactus_minimap_aligner -c 3 foobar

push:
	# Push our image to dockerhub for running in k8s
	#NOTE: do we want the rrounthrl/ubuntu image? or ubuntu:18:04? Why rrounthrl?
	docker tag $(USER)-cactus-minimap-aligner $(DOCKERHUB_ACCOUNT)/cactus-minimap-aligner
	#NOTE: do we want the rrounthrl/ubuntu image? or ubuntu:18:04? Why rrounthrl?
	docker push $(DOCKERHUB_ACCOUNT)/cactus-minimap-aligner

run-job:
	# Run a kubernetes job with our image, prefix with USERNAME and timestamp
	TS=`date +"%Y%m%d-%H%M%S"` envsubst < job.yml | kubectl create -f -

delete-my-jobs:
	# Delete jobs prefixed with USERNAME
	kubectl get jobs -o custom-columns=:.metadata.name \
		| grep '^$(USER)*' | xargs kubectl delete jobs

update-secrets:
	# Update secrets from our AWS file so we can access S3 in k8s
	kubectl delete secrets/shared-s3-credentials
	kubectl create secret generic shared-s3-credentials --from-file=credentials=../cgl-shared-s3-credentials
