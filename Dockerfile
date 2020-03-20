#NOTE: to configure aws, run your command like docker run -e AWS_DEFAULT_REGION='[your region] -e AWS_ACCESS_KEY_ID='[your access ID] -e AWS_SECRET_ACCESS_KEY='[your access key] rrounthw/minimap_toil
FROM python
# FROM ubuntu:18.04 Note: for some reason, toil install doesn't work on ubuntu base.

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
  wget \
  curl \
  time \
  python3 \
  python3-pip \
  python3-setuptools \
  && rm -rf /var/lib/apt/lists/*

#install aws CLI
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
RUN unzip awscliv2.zip
RUN ./aws/install
# role_arn=arn:aws:iam::781907127277:role/developer
# mfa_serial=arn:aws:iam::652235167018:mfa/rrounthw@ucsc.edu

RUN curl -L https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 | tar -jxvf -
RUN mv /minimap2-2.17_x64-linux/minimap2 /usr/bin/ 

RUN pip3 install --upgrade pip

WORKDIR /app
ADD requirements.txt /requirements.txt
RUN pip3 install --no-cache-dir -r /requirements.txt

ADD minimap2_cactus_aligner_toil.py /app

WORKDIR /app
# ENTRYPOINT time python3 minimap2_cactus_aligner_toil.py --assemblies_dir assembly_contigs --ref_file hg38.fa --output_file alignments_10k_context_20_mapq_cutoff_2_remap_thresh_100.sam file:my-job-store && \
#   aws s3 cp alignments_10k_context_20_mapq_cutoff_2_remap_thresh_100.sam s3://vg-k8s/users/rrounthw/large_1_alignments_10k_context_20_mapq_cutoff_2_remap_thresh_100.sam --profile vg-developer



# # # small run
# ADD hg38_chr21.fa /app
# ADD small_chr21/assemblies /app/assembly_contigs
# COPY small_chr21/assemblies /app/assembly_contigs
# ENTRYPOINT time python3 minimap2_cactus_aligner_toil.py --assemblies_dir assembly_contigs --ref_file hg38_chr21.fa --output_file alignments_10k_context_20_mapq_cutoff_2_remap_thresh_100.sam file:my-job-store
# ENTRYPOINT time python3 minimap2_cactus_aligner_toil.py --assemblies_dir assembly_contigs --ref_file hg38_chr21.fa --output_file alignments_10k_context_20_mapq_cutoff_2_remap_thresh_100.sam file:my-job-store && \
  # aws s3 cp alignments_10k_context_20_mapq_cutoff_2_remap_thresh_100.sam s3://vg-k8s/users/rrounthw/small_alignments_10k_context_20_mapq_cutoff_2_remap_thresh_100.sam
# ENTRYPOINT aws s3 cp minimap2_cactus_aligner_toil.py s3://vg-k8s/users/rrounthw/minimap2_cactus_aligner_toil.py

# ENTRYPOINT time python3 minimap2_cactus_aligner_toil.py file:my-job-store









###################
# For running the command without a shell (which doesn't work with time, I think because it's an environment variable)
# ENTRYPOINT ["time", "python3", "minimap2_cactus_aligner_toil.py", "file:my-job-store"]
