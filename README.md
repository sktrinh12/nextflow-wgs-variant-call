# Germline Variant Calling Pipeline (GATK4 + Nextflow)

## Overview

This pipeline performs **germline short variant discovery (SNPs and INDELs)** from **paired-end human whole-genome sequencing (WGS) data** following the **GATK4 Best Practices** [workflow](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-)

The workflow includes:

 → FASTQ - Read QC
 → BWA-MEM - Alignment
 → MarkDuplicates - Duplicate marking
 → BQSR - Base quality score recalibration
 → HaplotypeCaller (-ERC GVCF) - Variant calling
 → GenomicsDBImport
 → GenotypeGVCFs - Variant selection (SNPs / INDELs)
 → VQSR (or CNN / hard filters)

```bash
FASTQ
 └─ BWA-MEM
     └─ MarkDuplicates
         └─ BQSR
             └─ HaplotypeCaller (per-sample, sharded)
               └─ Merge per-sample GVCFs
                   └─ GenomicsDBImport (interval-aware, cohort)
                       └─ GenotypeGVCFs
                           └─ GatherVcfs (if needed)
                               └─ FILTERING (CNN or VQSR)
                                   └─ Funcotator
                                       └─ FINAL VCF

### FILTERING
cohort.raw.vcf.gz
 └─ VariantRecalibrator (SNP)
     └─ ApplyVQSR (SNP)
         └─ VariantRecalibrator (INDEL)
             └─ ApplyVQSR (INDEL)
```

---

## Preliminary Steps

### 1. Install requirements

Tools needed:

- **Nextflow** (≥ 22.x)
- **Docker** (for container execution)
- **kubectl** (if using k8s executor)
- **aws cli**

Check:

```bash
nextflow -version
docker --version
aws --version
kubectl version
```

### 2. Download example input data

The original bash script uses paired-end reads from the 1000 Genomes Project.

Run the following commands before starting nextflow since only need to do once

```bash
wget -P reads ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_1.filt.fastq.gz
wget -P reads ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_2.filt.fastq.gz
```

### 3. Download reference and known-sites files

```bash
# Reference genome
wget -P assets https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

# Known sites for BQSR
wget -P assets https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget -P assets https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx
```

### 4. Use custom docker image to index ref file

Also convert the already downloaded `Homo_sapiens_assembly38 vcf file to vcf.gz` to `tbi`

```bash
docker run --rm \
  -v "$(pwd)/reference:/ref" \
  -w /ref \
  sktrinh12/bionfox-tools:1.1 \
  bwa index hg38.fa
  samtools faidx hg38.fa
  bgzip -c Homo_sapiens_assembly38.dbsnp138.vcf > Homo_sapiens_assembly38.dbsnp138.vcf.gz
  gatk CreateSequenceDictionary -R hg38_chr22.fa
  gatk IndexFeatureFile -I Homo_sapiens_assembly38.dbsnp138.vcf.gz
```

If testing locally, may want to subset the reference human genome for testing, so `bwa_mem` doesn't need so much memory (OOM error)

Extract chr20, chr21, chr22 for slightly more realistic testing

```bash
docker run --rm \
  -v "$(pwd)/reference:/ref" \
  -w /ref \
  sktrinh12/bionfox-tools:1.1 \
  samtools faidx hg38.fa chr20 chr21 chr22 > hg38_chr20-22.fa
  bwa index hg38_chr20-22.fa
  bgzip -c Homo_sapiens_assembly38.dbsnp138.vcf > Homo_sapiens_assembly38.dbsnp138.vcf.gz
  gatk CreateSequenceDictionary -R hg38_chr22.fa
  gatk IndexFeatureFile -I Homo_sapiens_assembly38.dbsnp138.vcf.gz
```

In addition, for testing purposes, subset the number of lines of the fastq file

```bash
docker run --rm \
  -v "$(pwd)/reads:/reads" \
  -w /reads \
  sktrinh12/bionfox-tools:1.1 \
  bash -c "seqkit head -n 50000 SRR062634_1.filt.fastq.gz > test_1.fastq.gz &&  seqkit head -n 50000 SRR062634_2.filt.fastq.gz > test_2.fastq.gz"
```

### Download the reference files

```bash
# Define the base URL
BASE_FTP="ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle"

# HapMap 3.3
curl -O "${BASE_FTP}/hg38/hapmap_3.3.hg38.vcf.gz"
curl -O "${BASE_FTP}/hg38/hapmap_3.3.hg38.vcf.gz.tbi"

# 1000G Omni 2.5
curl -O "${BASE_FTP}/hg38/1000G_omni2.5.hg38.vcf.gz"
curl -O "${BASE_FTP}/hg38/1000G_omni2.5.hg38.vcf.gz.tbi"

# 1000G Phase 1 High Confidence SNPs
curl -O "${BASE_FTP}/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
curl -O "${BASE_FTP}/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi"

# funcotator bundle download
curl -O "${BASE_FTP}/funcotator/funcotator_dataSources.v1.7.20200521g.tar.gz"
```

### Upload to S3

Use the `--exclude` flag for `"*chr22"` only if a test run was done

```bash
aws s3 sync reference s3://${BUCKET_NAME}/nextflow/reference/ \
    --exclude "*chr22*" \
    --exclude ".nextflow*" \
    --exclude "funcotator*.gz" \
    --dryrun \
    --profile ${AWS_PROFILE_NAME}
```

### Directory Layout

Recommended project structure:

```bash
.
├── reference
│  ├── hg38.fa
│  ├── Homo_sapiens_assembly38.dbsnp138.vcf
│  └── Homo_sapiens_assembly38.dbsnp138.vcf.idx
├── reads
│    ├── SRR062634_1.fastq.gz
│    └── SRR062634_2.fastq.gz
├── main.nf
├── modules
│  └── local
│     └── scripts.nf
├── nextflow.config
└── results
```

### Running the Pipeline

Local execution (Docker)

```bash
nextflow run main.nf -profile local -resume
```

### AWS Batch execution

Make sure:
* An AWS Batch queue exists

```bash
nextflow run main.nf -profile awsbatch -resume
```

* An S3 bucket is configured as the Nextflow work directory

### K8s execution

```bash
kubectl exec -it nextflow-launcher -- sh
```

### Output

Results are written to the `results/` directory:

```bash
/workspace/results/
├── variants/
│   └── cohort.final.annotated.vcf.gz
├── vcf_stats/
│   ├── cohort.stats
│   └── plots/
│       ├── depth.png
│       ├── substitutions.png
│       ├── indels.png
│       └── ...
```

---

## Kubernetes setup

Install the EFS CSI driver and configure EFS so pods can use it as shared storage for Nextflow scratch/workspace files. Not needed if using `awsbatch`, only for `k8s` executor.

The EKS cluster must already have:

* An EFS file system created
* Mount targets in each AZ
* An IAM role configured for IRSA (`AmazonEKS_EFS_CSI_DriverRole`)

Install the EFS CSI driver using IRSA (do not let Helm create the service account).

```bash
kubectl create serviceaccount efs-csi-controller-sa \
  -n kube-system \
  --dry-run=client -o yaml | kubectl apply -f -
```

Annotate the service account with the IAM role used by the EFS CSI controller.

```bash
kubectl annotate serviceaccount efs-csi-controller-sa \
  -n kube-system \
  eks.amazonaws.com/role-arn=arn:aws:iam::<ACCOUNT_ID>:role/AmazonEKS_EFS_CSI_DriverRole
```

Install the EFS CSI driver via Helm, reusing the existing service account.

```bash
helm upgrade --install aws-efs-csi-driver \
  --namespace kube-system \
  aws-efs-csi-driver/aws-efs-csi-driver \
  --set controller.serviceAccount.create=false \
  --set controller.serviceAccount.name=efs-csi-controller-sa
```

Create a StorageClass backed by EFS using access points.

```bash
kubectl apply -f - <<EOF
apiVersion: storage.k8s.io/v1
kind: StorageClass
metadata:
  name: efs-sc
provisioner: efs.csi.aws.com
parameters:
  provisioningMode: efs-ap
  fileSystemId: <FROM_TERRAFORM_OUTPUT>
  directoryPerms: "700"
EOF
```

Create a temporary pod with the EFS volume mounted, used only for staging data.

```bash
kubectl apply -f - <<EOF
apiVersion: v1
kind: Pod
metadata:
  name: efs-uploader
spec:
  restartPolicy: Never
  containers:
  - name: uploader
    image: amazonlinux:2
    command: ["sleep", "infinity"]
    volumeMounts:
    - name: workspace
      mountPath: /workspace
  volumes:
  - name: workspace
    persistentVolumeClaim:
      claimName: nextflow-pvc
EOF
```

Wait for the pod to be ready.

```bash
kubectl wait pod/efs-uploader --for=condition=Ready
```

Copy local read files into EFS by streaming them through `kubectl exec`.

```bash
for file in reads/SRR*; do
  kubectl exec -i efs-uploader -- sh -c \
    "cat > /workspace/reads/$(basename "$file")" < "$file"
done
```

Install utilities in the pod and extract reference data into EFS.

```bash
kubectl exec efs-uploader -- yum install -y tar gzip

tar \
  --exclude='reference/.*' \
  --exclude='reference/hg38_chr22*' \
  --exclude='reference/funcotator*.tar.gz' \
  -czf - reference | \
kubectl exec -i efs-uploader -- tar \
  --no-same-owner \
  --no-same-permissions \
  -xzf - -C /workspace
```

#### Misc `kubectl` commands to understand and observe the pods whilst executing tasks

```bash
# delete pods that have completed or in not pending or running state
kubectl get pods -l nextflow.io/app=nextflow,nextflow.io/processName=BWA_MEM \
  --field-selector=status.phase!=Running,status.phase!=Pending \
  -o name | xargs -r kubectl delete

# get pods that are running
kubectl get pods --field-selector=status.phase=Running

# get pods that are running in table format
kubectl get pods -o json | jq -r '
.items[] |
select(
  .status.containerStatuses[0].state.terminated?.reason != "Error"
) |
"\(.metadata.name)\t\(.status.phase)"'

# view pods for particular nextflow process
kubectl get pods -l nextflow.io/app=nextflow,nextflow.io/processName=CREATE_INTERVALS

# confirm that nextflow applies the proper cpu/mem requests for the pods
k get po -l 'nextflow.io/app=nextflow,nextflow.io/processName=HAPLOTYPE_CALLER' | awk 'NR==2{print $1}' | xargs -I {} kubectl describe po {} | grep -A 10 "Requests:"
```

This approach above was an initial attempt to stage data into EFS from a local machine. It works, but is slow, manual, and not ideal for large datasets or reproducibility.

The next section replaces this approach with a cleaner and more scalable pattern using S3 as the source of truth and letting the pipeline pull data directly from S3 instead.

Optional cleanup if the driver is no longer needed:

```bash
k delete -f pvc.yaml
k delete -f nf-launcher.yaml
helm uninstall aws-efs-csi-driver --namespace kube-system
```

---

### AWS Batch setup for Nextflow

This section documents the minimal AWS-side setup used to run a Nextflow pipeline with AWS Batch as the executor.

#### 1. Identify VPC and usable subnets

Get the default VPC ID:

```bash
VPC_ID=$(aws ec2 describe-vpcs \
  --filters Name=isDefault,Values=true \
  --query 'Vpcs[0].VpcId' \
  --output text \
  --profile ${AWS_PROFILE_NAME})
```

List all subnets in that VPC:

```bash
aws ec2 describe-subnets \
  --filters Name=vpc-id,Values=$VPC_ID \
  --query 'Subnets[].SubnetId' \
  --output text \
  --profile ${AWS_PROFILE_NAME}
```

Pick a subnet that has outbound internet access (required for pulling container images, accessing S3, etc.). Check that the subnet is associated with a route table that has a default route to an Internet Gateway (public) or has a NAT (private)
- can select multiple, since the subnet flag accepts a list

```bash
aws ec2 describe-route-tables \
  --filters Name=association.subnet-id,Values=<SUBNET_ID> \
  --query 'RouteTables[].Routes[?DestinationCidrBlock==`0.0.0.0/0`]' \
  --profile ${AWS_PROFILE_NAME}
```

Use one such subnet in all Batch compute environments.

---

#### 2. Create a security group for Batch compute nodes

Create a security group in the selected VPC:

```bash
aws ec2 create-security-group \
  --group-name nextflow-batch-sg \
  --description "SG for Nextflow AWS Batch compute nodes" \
  --vpc-id $VPC_ID \
  --profile ${AWS_PROFILE_NAME}
```

Allow all outbound traffic (use the output from the command above for `group-id`):

```bash
aws ec2 authorize-security-group-egress \
  --group-id <SECURITY_GROUP> \
  --protocol -1 \
  --cidr 0.0.0.0/0 \
  --profile ${AWS_PROFILE_NAME}
```

---

#### 3. Create IAM roles

Batch service role (used by AWS Batch itself):

```bash
aws iam create-role \
  --role-name AWSBatchServiceRole \
  --assume-role-policy-document '{
    "Version": "2012-10-17",
    "Statement": [{
      "Effect": "Allow",
      "Principal": { "Service": "batch.amazonaws.com" },
      "Action": "sts:AssumeRole"
    }]
  }' \
  --profile ${AWS_PROFILE_NAME}

aws iam attach-role-policy \
  --role-name AWSBatchServiceRole \
  --policy-arn arn:aws:iam::aws:policy/service-role/AWSBatchServiceRole \
  --profile ${AWS_PROFILE_NAME}
```

EC2 instance role for Batch compute nodes:

```bash
aws iam create-role \
  --role-name NextflowBatchInstanceRole \
  --assume-role-policy-document '{
    "Version": "2012-10-17",
    "Statement": [{
      "Effect": "Allow",
      "Principal": { "Service": "ec2.amazonaws.com" },
      "Action": "sts:AssumeRole"
    }]
  }' \
  --profile ${AWS_PROFILE_NAME}
```

Attach required policies (minimum for ECS, ECR, and S3 access):

```bash
aws iam attach-role-policy \
  --role-name NextflowBatchInstanceRole \
  --policy-arn arn:aws:iam::aws:policy/service-role/AmazonEC2ContainerServiceforEC2Role \
  --profile ${AWS_PROFILE_NAME}

aws iam attach-role-policy \
  --role-name NextflowBatchInstanceRole \
  --policy-arn arn:aws:iam::aws:policy/AmazonEC2ContainerRegistryReadOnly \
  --profile ${AWS_PROFILE_NAME}

aws iam attach-role-policy \
  --role-name NextflowBatchInstanceRole \
  --policy-arn arn:aws:iam::aws:policy/AmazonS3FullAccess \
  --profile ${AWS_PROFILE_NAME}
```

Create an instance profile and add the role (required by Batch):

```bash
aws iam create-instance-profile \
  --instance-profile-name NextflowBatchInstanceProfile \
  --profile ${AWS_PROFILE_NAME}

aws iam add-role-to-instance-profile \
  --instance-profile-name NextflowBatchInstanceProfile \
  --role-name NextflowBatchInstanceRole \
  --profile ${AWS_PROFILE_NAME}
```

---

#### 4. Create an EC2 launch template

Batch EC2 instances need a larger root volume for Nextflow work directories:

```bash
aws ec2 create-launch-template \
  --launch-template-name nf-ebs-disk \
  --launch-template-data '{
    "BlockDeviceMappings": [
      {
        "DeviceName": "/dev/xvda",
        "Ebs": {
          "VolumeSize": 200,
          "VolumeType": "gp3",
          "DeleteOnTermination": true,
          "Throughput": 250,
          "Iops": 3000
        }
      }
    ]
  }' \
  --profile ${AWS_PROFILE_NAME}
```

---

#### 5. Create AWS Batch compute environments

On-demand compute environment:

```bash
aws batch create-compute-environment \
  --compute-environment-name nf-ondemand-ce \
  --type MANAGED \
  --state ENABLED \
  --service-role arn:aws:iam::<ACCOUNT_ID>:role/AWSBatchServiceRole \
  --compute-resources '{
    "type": "EC2",
    "minvCpus": 0,
    "maxvCpus": 50,
    "desiredvCpus": 0,
    "instanceTypes": ["c6i.large", "c6i.xlarge", "c6i.2xlarge"],
    "subnets": ["<SUBNET_ID>"],
    "securityGroupIds": ["<SECURITY_GROUP>"],
    "instanceRole": "arn:aws:iam::<ACCOUNT_ID>:instance-profile/NextflowBatchInstanceProfile",
    "launchTemplate": {
      "launchTemplateName": "nf-ebs-disk"
    }
  }' \
  --profile ${AWS_PROFILE_NAME}
```

Spot compute environment:

```bash
aws batch create-compute-environment \
  --compute-environment-name nf-spot-ce \
  --type MANAGED \
  --state ENABLED \
  --service-role arn:aws:iam::<ACCOUNT_ID>:role/AWSBatchServiceRole \
  --compute-resources '{
    "type": "SPOT",
    "minvCpus": 0,
    "maxvCpus": 100,
    "desiredvCpus": 0,
    "instanceTypes": ["c6i.large", "c6i.xlarge", "c6i.2xlarge"],
    "subnets": ["<SUBNET_ID>"],
    "securityGroupIds": ["<SECURITY_GROUP>"],
    "instanceRole": "arn:aws:iam::<ACCOUNT_ID>:instance-profile/NextflowBatchInstanceProfile",
    "allocationStrategy": "SPOT_CAPACITY_OPTIMIZED",
    "bidPercentage": 80,
    "launchTemplate": {
      "launchTemplateName": "nf-ebs-disk"
    }
  }' \
  --profile ${AWS_PROFILE_NAME}
```

---

#### 6. Create a job queue

Create a job queue pointing at the Spot environment:

```bash
aws batch create-job-queue \
  --job-queue-name nf-spot-queue \
  --state ENABLED \
  --priority 1 \
  --compute-environment-order order=1,computeEnvironment=nf-spot-ce \
  --profile ${AWS_PROFILE_NAME}
```

This job queue is what Nextflow references in `nextflow.config`.

---

#### 7. Useful inspection and cleanup commands

List compute environments:

```bash
aws batch describe-compute-environments --profile ${AWS_PROFILE_NAME}
```

List jobs in a queue with readable timestamps:

```bash
aws batch list-jobs --job-queue nf-spot-queue --profile ${AWS_PROFILE_NAME} | \
  jq '
    .jobSummaryList[] |
    (.createdAt |= (. / 1000 | todate)) |
    (.startedAt |= (if . then (. / 1000 | todate) else . end))
  '
```

Disable and delete resources when cleaning up:

```bash
aws batch update-job-queue --job-queue nf-spot-queue --state DISABLED --profile ${AWS_PROFILE_NAME}
aws batch delete-job-queue --job-queue nf-spot-queue --profile ${AWS_PROFILE_NAME}

aws batch update-compute-environment --compute-environment nf-spot-ce --state DISABLED --profile ${AWS_PROFILE_NAME}
aws batch delete-compute-environment --compute-environment nf-spot-ce --profile ${AWS_PROFILE_NAME}
```
