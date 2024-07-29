(usage)=

# Usage

Assuming that you've followed the {ref}`installations steps <installation>`, you're now ready to use this package.

Start with entering the container's command line by executing:

```docker
docker exec -it $(docker ps -aqf "ancestor=sungshic/gsap") /bin/bash
```

Or, by using the following, if the image was rebuilt locally:

```docker
docker exec -it $(docker ps -aqf "ancestor=gsap") /bin/bash
```

In the container's command line, execute the following to run the GSAP pipeline against the example data:

```bash
python src/gsap -F ./data_input/sample/SP2_S6_L001_R1_001.fastq.gz -R ./data_input/sample/SP2_S6_L001_R2_001.fastq.gz -N testgenome -A ./data_input/refseq/AL009126_v11/AL009126.fasta -B ./data_input/refseq/AL009126_v11/AL009126.gb -C "hello" -o "B. subtilis" -m "DNA" -O ./data_output/out.gb -T ./data_output/
```

Now, it would take some time for the pipeline to fully perform the sequence assembly, analysis, and annotation.
After a long wait, somewhere between 3 to 4 hours for the example shown here on a powerful laptop machine, the resulting Genbank file will appear under the designated path. In this example:

```bash
ls -lh data_output/out.gb
```
