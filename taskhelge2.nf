nextflow.enable.dsl = 2
params.storedir = "${baseDir}/cache"
params.outdir = "${baseDir}/outf"
params.catdir = "${baseDir}"
//params.url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=M21012&rettype=fasta&retmode=text"
params.url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id="
params.accession = "M21012"
params.outfile =  params.accession.concat(".fasta")
params.appendix = "&rettype=fasta&retmode=text"
params.git =  "https://gitlab.com/dabrowskiw/cq-examples.git"

process downloadRefs {
    storeDir params.outdir
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    input:
        val url
    output:
        path params.outfile
    script:     
    """
    wget "${url}" -O  ${params.outfile}
    """
}//downloads the reference sequence as specified by accession and 


process catAcc {   //concatenates reference file (accession) to samples 
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
  
    script:          
  """
    cat ${params.catdir}/outf/seq*.fasta >  ${params.catdir}/outf/collect.fasta 
    cat ${params.catdir}/outf/${params.outfile} ${params.catdir}/outf/collect.fasta > "${params.catdir}/outf/allfasta.fasta"
  """
  }

process runMafft {
    storeDir "${params.outdir}"
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    container "https://depot.galaxyproject.org/singularity/mafft%3A7.505--hec16e2b_0" 
input: 
   path fastafile
output: 
   path "${fastafile}.m"
script:    
"""
mafft --retree 2 --maxiterate 1000  ${fastafile} > "${fastafile}.m"
""" 
}   //singularity run https://depot.galaxyproject.org/singularity/mafft%3A7.505--hec16e2b_0

process runTrimal {
    storeDir "${params.outdir}"
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    container "https://depot.galaxyproject.org/singularity/trimal%3A1.4.1--hc9558a2_4"
input: 
   path fastamfile
output: 
   path "${fastamfile}.t" 
   path "${fastamfile}.html"

script: 
"""
trimal -in ${fastamfile} -out  "${fastamfile}.t"  -automated1 -htmlout "${fastamfile}.html"
"""
}

workflow {
    url = params.url + params.accession + params.appendix
    reference = downloadRefs(url)
    catAcc()
    fastainfile = Channel.fromPath("${params.catdir}/outf/allfasta.fasta")
    fastamfile = runMafft(fastainfile)
    runTrimal(fastamfile)    
}

/*
N E X T F L O W  ~  version 22.04.5
Launching `taskhelge2.nf` [suspicious_koch] DSL2 - revision: cd5a61b76e
executor >  local (3)
[skipped  ] process > downloadRefs  [100%] 1 of 1, stored: 1 ✔
[6e/ae2dea] process > catAcc        [100%] 1 of 1 ✔
[2e/e8d644] process > runMafft (1)  [100%] 1 of 1 ✔
[8f/279e14] process > runTrimal (1) [100%] 1 of 1 ✔
[skipping] Stored process > downloadRefs
*/
