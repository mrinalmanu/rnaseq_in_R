TAGS=$(ls fastqs/SRX*.fastq.gz | xargs -n 1 basename | sed 's/.fastq.gz//')

for TAG in $TAGS; do
  OUTDIR="fastqc/$TAG"; mkdir -p "$OUTDIR" 
  fastqc -o "$OUTDIR" "fastqs/$TAG.fastq.gz" |& tee "$OUTDIR/$TAG.fastqc.log"
done
------------------------------------------
TAGS=$(ls fastqs/SRX31956*.fastq.gz | xargs -n 1 basename | sed 's/_[1,2].fastq.gz//' | uniq)
for TAG in $TAGS; do
  HISAT_IDX=/mnt/reference/Gencode_mouse/release_M20/GRCm38.primary_assembly

  # aligning to the genome reference

  OUTDIR="hisat2/$TAG"; mkdir -p "$OUTDIR"
  date
  hisat2 -p 8 -x ${HISAT_IDX} \
  -1 "fastqs/$TAG*_1.fastq.gz" -2 "fastqs/$TAG*_2.fastq.gz" \
  2> "$OUTDIR/$TAG.hisat2.log" \
  | samtools view -bS - > "$OUTDIR/$TAG.raw.bam"
  date
done
--------------------------------------------
TAGS=$(ls fastqs/SRX31956*.fastq.gz | xargs -n 1 basename | sed 's/_[1,2].fastq.gz//' | uniq)
for TAG in $TAGS; do
  OUTDIR="hisat2/$TAG"
  date
  samtools sort -@ 8 -O bam "$OUTDIR/$TAG.raw.bam" > "$OUTDIR/$TAG.bam" && \
  samtools index "$OUTDIR/$TAG.bam" && \
  rm -v "$OUTDIR/$TAG.raw.bam"
  date
done
--------------------------------------------
TAGS=$(ls fastqs/SRX31956*.fastq.gz | xargs -n 1 basename | sed 's/_[1,2].fastq.gz//' | uniq)
for TAG in $TAGS; do
  OUTDIR="hisat2/$TAG"
  bamCoverage -b "$OUTDIR/$TAG.bam" -o "$OUTDIR/$TAG.cov.bw" |& tee "$OUTDIR/$TAG.bamcov.log"

done
---------------------------------------------
TAGS=$(ls fastqs/SRX31956*.fastq.gz | xargs -n 1 basename | sed 's/_[1,2].fastq.gz//' | uniq)
for TAG in $TAGS; do
  OUTDIR="hisat2/$TAG"
  REFGENE_MODEL=/mnt/reference/Gencode_mouse/release_M20/mm10_Gencode_VM18.bed
  infer_experiment.py -i "$OUTDIR/$TAG.bam" \
  -r $REFGENE_MODEL | tee "$OUTDIR/$TAG.infer_experiment.txt"

  REFGENE_MODEL=/mnt/reference/Gencode_mouse/release_M20/mm10_Gencode_VM18.bed
  read_distribution.py -i "$OUTDIR/$TAG.bam" \
  -r $REFGENE_MODEL | tee "$OUTDIR/$TAG.read_distribution.txt"
done
---------------------------------------------
TAGS=$(ls fastqs/SRX31956*.fastq.gz | xargs -n 1 basename | sed 's/_[1,2].fastq.gz//' | uniq)
for TAG in $TAGS; do
  GTF=/mnt/reference/Gencode_mouse/release_M20/gencode.vM20.annotation.gtf

  OUTDIR="featureCounts/$TAG"; mkdir -p "$OUTDIR"
  date
  featureCounts -a "$GTF" -s 0 -p -o "$OUTDIR/$TAG.fc.txt" \
  "hisat2/$TAG/$TAG.bam" |& tee "$OUTDIR/$TAG.fc.log"
  date

  head "$OUTDIR/$TAG.fc.txt"
  wc -l "$OUTDIR/$TAG.fc.txt"
done
---------------------------------------------
mkdir kallisto
for TAG in $TAGS; do
  KALLISTO_IDX=/mnt/reference/Gencode_mouse/release_M20/gencode.vM20.transcripts.kalliso.idx

  OUTDIR="kallisto/$TAG"; mkdir -p "$OUTDIR"
  date
  # —single -l and -s option should be set for each dataset separately, 200+-50 is most common for single end
  kallisto quant -i $KALLISTO_IDX -t 8 \
  —single -l 200 -s 50 \
  —plaintext \
  -o $OUTDIR \
  fastqs/$TAG.fastq.gz |& tee $OUTDIR/$TAG.kallisto.log
  date
done
---------------------------------------------
multiqc . -x -f 

OUTDIR="mmquant"; mkdir -p "$OUTDIR"
GTF=/mnt/reference/Gencode_mouse/release_M20/gencode.vM20.annotation.gtf
date
mmquant -a "$GTF" -s U -o "$OUTDIR/mmq.txt" \
-r hisat2/*/*.bam |& tee "$OUTDIR/mmq.log"
date