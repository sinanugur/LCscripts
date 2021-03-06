directories, files, = glob_wildcards("data/raw/{samplegroup}/{sample}.fastq.gz")

humangenome="/home/sium/data/humangenome/hg38"

# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

rule all:
	input:
		expand("analyses/bowtie_mappings_genome_multi/{samplegroup}/{sample}.truncated.collapsed.fasta.gz.hg38.local.multi.sorted.bam",zip,samplegroup=directories,sample=files),
		"classification/miRJanus_file_statistics.csv",
		expand("data/stats/{samplegroup}/{sample}.truncated.FastQC",zip,samplegroup=directories,sample=files)
		#expand("analyses/bam_to_bed/{samplegroup}/{sample}.truncated.collapsed.fasta.gz.hg38.local.multi.sorted.bed.gz",zip,samplegroup=directories,sample=files)

rule adapter_trimming:
	input:
		"data/raw/{samplegroup}/{sample}.fastq.gz"

	output:
		"data/trimmed/{samplegroup}/{sample}.truncated.fastq.gz"
	
	threads: 15

	shell:
		# "Adapter_removal.sh {input} {output}"
		"""
		out={output}
		basename=${{out/.truncated.fastq.gz/}}
		AdapterRemoval --file1 {input} --basename $basename --trimns --trimqualities --gzip --threads {threads} --minlength 17 --minquality 30;
		mv "$basename".truncated.gz "$basename".truncated.fastq.gz;
		mv "$basename".discarded.gz "$basename".discarded.fastq.gz;
		"""


rule collapsing_reads:
	input:
		"data/trimmed/{samplegroup}/{sample}.truncated.fastq.gz"
	
	output:
		"data/collapsed/{samplegroup}/{sample}.truncated.collapsed.fasta.gz"

	shell:
		"""
		#collapse_reads.sh {input} {output} 
		zcat {input} | fastx_collapser | gzip > {output}
				
		"""

rule bowtie2_mapping:
	input:
		"data/collapsed/{samplegroup}/{sample}.truncated.collapsed.fasta.gz"

	output:
		"analyses/bowtie_mappings_genome_multi/{samplegroup}/{sample}.truncated.collapsed.fasta.gz.hg38.local.multi.sorted.bam"

	threads: 15

	shell:
		#"bowtie_map_local_all.sh /home/sium/data/humangenome/hg38 {input} {output}"
		"""
		bowtie2 --sensitive-local -k 10 -f -p {threads} -x {humangenome} -U <(zcat {input}) | samtools view -bS - | samtools sort - -o {output}
		"""

rule file_stats:
	input:
		"data/raw/{samplegroup}/{sample}.fastq.gz",
		"data/trimmed/{samplegroup}/{sample}.truncated.fastq.gz",
		"data/collapsed/{samplegroup}/{sample}.truncated.collapsed.fasta.gz"

	output:
		"data/stats/{samplegroup}/{sample}.sample.stats.txt"

	shell:
		"""
		samplename=$(sample_name.sh {input[2]})
		type=$(gawk -v samplename=$samplename '{{if(match($0, " "samplename" ")) print $3;}}' /home/sium/projects/miRJanus/classification/all.samples.latest.txt)
		
		md5=$(md5sum {input[0]} | awk '{{print $1}}')
		seq_count=$(zcat {input[0]} | echo $((`wc -l`/4)))
		trimmed_seq_count=$(zcat {input[1]} | echo $((`wc -l`/4)))
		collapsed_reads=$(zcat {input[2]} | grep ">" -c)
		spike_in_counts=$(zcat {input[2]} | awk '{{if(/^TCACCGGGTGTAAATCAGCTTG$/) print line; else line=$0;}}' | awk '/>/{{match($0,/-([0-9]+)$/,m); print m[1]}}')
	
		filename=$(basename {input[0]})
		gc=$(seqtk fqchk {input[0]} | awk '/ALL/{{print $4+$5}}')
		if [ -L {input[0]} ]; then physical=$(readlink {input[0]}); else physical={input[0]}; fi
		
		filesize=$(du -h --block-size=M $physical | cut -f1)

		echo -e "md5sum\tPhysical Location\tFilename\tCase-Control\tSample ID\t#Reads\t#Reads After Trimming\t#Collapsed Reads\t#Spike In Counts\t#File Size(MB)\tGC" > {output} 
		echo -e $md5"\t"$physical"\t"$filename"\t"$type"\t"$samplename"\t"$seq_count"\t"$trimmed_seq_count"\t"$collapsed_reads"\t"$spike_in_counts"\t"$filesize"\t"$gc >> {output}

		"""

rule fastqc_stats:
	input:
		"data/trimmed/{samplegroup}/{sample}.truncated.fastq.gz"

	output:
		directory("data/stats/{samplegroup}/{sample}.truncated.FastQC")

	threads: 15

	shell:
		"""
		mkdir {output};
		fastqc --nogroup --extract -o {output} -t {threads} {input};
		"""

rule combine_file_stats:
	input:
		expand("data/stats/{samplegroup}/{sample}.sample.stats.txt",zip,samplegroup=directories,sample=files)

	output:
		"classification/miRJanus_file_statistics.csv"
	
	run:
		shell("""echo "md5sum\tphysical\tfilename\tcase_control\tsample\ttotal_reads\ttotal_reads_after_trimming\tcollapsed_reads\tspike_in_counts\tfile_size(MB)\tgc" > {output};""")	
		for i in input:
			shell("""cat {i} | grep -v "#" >> {output};""")


rule bam_to_bed:
	input:
		"analyses/bowtie_mappings_genome_multi/{samplegroup}/{sample}.truncated.collapsed.fasta.gz.hg38.local.multi.sorted.bam"
	output:
		"analyses/bam_to_bed/{samplegroup}/{sample}.truncated.collapsed.fasta.gz.hg38.local.multi.sorted.bed.gz"
	run:
		shell("bedtools bamtobed -i {input} | awk 'match($0,/[0-9]+-([0-9]+)/,m) {{for(i=1;i <= m[1];i++) print}}' | gzip > {output}")

