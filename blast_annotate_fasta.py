#!/usr/bin/env python

import os, sys
from Bio import SeqIO
import subprocess
from progress.bar import FillingSquaresBar
from datetime import datetime
from random import randint
import argparse

def main():
    options=argparse.ArgumentParser(sys.argv[0],
        description='Script for annotating multiple sequences in fasta format using a list of reference genomes in GenBank format (downloaded from NCBI or in NCBI-compliant format)',
        prefix_chars='-',
        add_help=True,
        epilog='Written by Chrispin Chaguza, Yale School of Public Health, Yale University, 2021')
    options.add_argument('-s','--sequences',action='store',required=True,nargs=1,
        metavar='input_sequences',dest='input_sequences',help='Input (multi-) fasta file containing nucleotide sequences to be annotated')
    options.add_argument('-r','--references',action='store',required=True,nargs=1,
        metavar='input_references',dest='input_references',help='Input file containing locations to the reference genomes to be used for annotation (one per line)')
    options.add_argument('-o','--output',action='store',required=True,nargs=1,
        metavar='output_annotations_file',dest='output_annotations_file',help='Output file containing the annotated nucleotide sequences')

    options=options.parse_args()

    seq_kmers=options.input_sequences[0:][0]
    gb_files=options.input_references[0:][0]
    output_rep_file=options.output_annotations_file[0:][0]

    rand_str_val="tmp."+str(randint(0,10000))+"."+str(randint(0,10000))

    tmp_blast_output = "tmp."+rand_str_val+".bl.txt"
    tmp_genome_names = "tmp."+rand_str_val+".gb.txt"
    tmp_rscript_file = "tmp."+rand_str_val+".R"

    with open(tmp_genome_names,"w") as tmp_file:
        max_prog_bar = 1 if len([loopCount for loopCount in open(gb_files,"r")])==0 else len([loopCount for loopCount in open(gb_files,"r")])
        status = FillingSquaresBar('['+str(datetime.now().strftime("%H:%M:%S"))+'] step 1/6', max=max_prog_bar)
        for j in [str(r).strip() for r in open(gb_files,"r")]:
            in_seq_data = SeqIO.read(j,"genbank")
            out_seq_data = open(str(in_seq_data.id)+".fasta","w")
            SeqIO.write(in_seq_data,out_seq_data,"fasta")
            tmp_file.write(str(in_seq_data.id)+".fasta\n")
            status.next()
    status.next()
    status.finish()

    cmd = []
    cmd.append("cat `cat")
    cmd.append(tmp_genome_names)
    cmd.append("|")
    cmd.append("tr \"\n\" \" \"`")
    cmd.append(">")
    cmd.append(rand_str_val+"db")
    cmd = " ".join(cmd)
    retval = subprocess.call(cmd,shell=True,stdout=open(os.devnull,'w'),stderr=subprocess.STDOUT)

    cmd = []
    cmd.append("makeblastdb")
    cmd.append("-in")
    cmd.append(rand_str_val+"db")
    cmd.append("-dbtype")
    cmd.append("nucl")
    cmd=" ".join(cmd)
    retval = subprocess.call(cmd,shell=True,stdout=open(os.devnull,'w'),stderr=subprocess.STDOUT)

    cmd = []
    cmd.append("rm -rf")
    cmd.append(tmp_blast_output)
    cmd = " ".join(cmd)
    retval = subprocess.call(cmd,shell=True,stdout=open(os.devnull,'w'),stderr=subprocess.STDOUT)
    
    max_prog_bar = 1 if len([loopCount for loopCount in SeqIO.parse(seq_kmers,"fasta")])==0 else len([loopCount for loopCount in SeqIO.parse(seq_kmers,"fasta")])
    status = FillingSquaresBar('['+str(datetime.now().strftime("%H:%M:%S"))+'] step 2/6', max=max_prog_bar)
    for r in SeqIO.parse(seq_kmers,"fasta"):
        with open(rand_str_val+"lp","w") as fhandle:
            fhandle.write(">"+str(r.id)+"\n"+str(r.seq))

        cmd = []
        cmd.append("blastn")
        cmd.append("-db")
        cmd.append(rand_str_val+"db")
        cmd.append("-query")
        cmd.append(rand_str_val+"lp")
        cmd.append("-task")
        cmd.append("blastn-short")
        cmd.append("-outfmt")
        cmd.append("6")
        cmd.append("-word_size")
        cmd.append("31")
        cmd.append("-num_threads")
        cmd.append("4")
        cmd.append("-max_target_seqs")
        cmd.append("2000")
        cmd.append("-evalue")
        cmd.append("0.0001")
        cmd.append("-perc_identity")
        cmd.append("90")
        cmd.append("-ungapped")
        cmd.append("|")
        cmd.append("sed")
        cmd.append("\"s/Query_1/"+str(r.id)+"/g\"")
        cmd.append(">>")
        cmd.append(tmp_blast_output)
        cmd = " ".join(cmd)
        retval = subprocess.call(cmd,shell=True,stdout=open(os.devnull,'w'),stderr=subprocess.STDOUT)
        status.next()
    status.next()
    status.finish()    

    annotfiles = []

    with open(gb_files,"r") as gb_fhandle:
        max_prog_bar = 1 if len([loopCount for loopCount in open(gb_files,"r")])==0 else len([loopCount for loopCount in open(gb_files,"r")])
        status = FillingSquaresBar('['+str(datetime.now().strftime("%H:%M:%S"))+'] step 3/6', max=max_prog_bar)
        for gb_name in gb_fhandle:
            annotfiles.append(str(gb_name).strip())
            status.next()
        status.next()
        status.finish()

    db={}
    
    max_prog_bar = 1 if len([loopCount for loopCount in annotfiles])==0 else len([loopCount for loopCount in annotfiles])
    status = FillingSquaresBar('['+str(datetime.now().strftime("%H:%M:%S"))+'] step 4/6', max=max_prog_bar)
    for h in annotfiles:
        db[str(h)] = SeqIO.read(h,"genbank")
        status.next()
    status.next()
    status.finish()

    blast_outfile = tmp_blast_output
    blast = []

    max_prog_bar = 1 if len([loopCount for loopCount in open(blast_outfile,"r")])==0 else len([loopCount for loopCount in open(blast_outfile,"r")])
    status = FillingSquaresBar('['+str(datetime.now().strftime("%H:%M:%S"))+'] step 5/6', max=max_prog_bar)
    with open(blast_outfile,"r") as blast_res:
        for i in blast_res:
            blast.append(str(i).strip().split("\t"))
            status.next()
    status.next()
    status.finish()

    fheader = []
    fheader.append("Variant")
    fheader.append("Genome")
    fheader.append("MatchLength")
    fheader.append("GenomeStart")
    fheader.append("GenomeEnd")
    fheader.append("GeneStart")
    fheader.append("GeneEnd")
    fheader.append("GeneID")
    fheader.append("LocusTag")
    fheader.append("Product")
    fheader.append("Evalue")
    fheader = "\t".join(fheader)

    with open(output_rep_file,"w") as blast_report:
        blast_report.write(fheader+"\n")
        max_prog_bar = 1 if len([loopCount for loopCount in blast])==0 else len([loopCount for loopCount in blast])
        status = FillingSquaresBar('['+str(datetime.now().strftime("%H:%M:%S"))+'] step 6/6', max=max_prog_bar)
        for j in blast:
            for i in annotfiles:
                annot = db[str(i)]
                if annot.id == j[1]:
                    check = True
                    for l,k in enumerate(annot.features):
                        if k.type == "CDS" and k.location.start.position!=0:
                            start = k.location.start.position
                            end = k.location.end.position
                            product = ""
                            gene = ""
                            if 'product' in k.qualifiers:
                                product = k.qualifiers['product'][0]
                            if 'gene' in k.qualifiers:
                                gene = k.qualifiers['gene'][0]
                            if 'locus_tag' in k.qualifiers:
                                locus=k.qualifiers['locus_tag'][0]
                            if ((int(j[8])<=start) and (int(j[9])>=start)) or ((int(j[8])>=start) and (int(j[8])<=end)):
                                outp = []
                                outp.append(str(j[0]))
                                outp.append(str(j[1]))
                                outp.append(str(j[3]))
                                outp.append(str(j[8]))
                                outp.append(str(j[9]))
                                outp.append(str(start))
                                outp.append(str(end))
                                outp.append(str(gene))
                                outp.append(str(locus))
                                outp.append(str(product))
                                outp.append(str(j[10]))
                                outp = "\t".join(outp)
                                blast_report.write(outp+"\n")
                                check=False
                            else:
                                pass
                        else:
                            pass
                    if check:
                        outp = []
                        outp.append(str(j[0]))
                        outp.append(str(j[1]))
                        outp.append(str(j[3]))
                        outp.append(str(j[8]))
                        outp.append(str(j[9]))
                        outp.append("N/A")
                        outp.append("N/A")
                        outp.append("")
                        outp.append("Intergenic")
                        outp.append("")
                        outp.append(str(j[10]))
                        outp = "\t".join(outp)
                        blast_report.write(outp+"\n")
                else:
                    pass
            status.next()
        status.next()
        status.finish()

    with open(tmp_rscript_file,"w") as rhandle:
        rhandle.write("#!/usr/bin/env Rscript\n")
        rhandle.write("library(dplyr)\n")
        rhandle.write("library(magrittr)\n")
        rhandle.write("library(tidyr)\n")
        rhandle.write("MM<-dplyr::as_tibble(read.csv(\""+str(output_rep_file)+"\",header=T,sep=\"\\t\",comment.char=\"?\"))\n")
        rhandle.write("MM %>% dplyr::group_by(Variant) %>% arrange(Product,Variant) %>% dplyr::mutate(Annot.group = data.table::rleid(Product)) %>%\n")
        rhandle.write("dplyr::mutate(Product = tidyr::replace_na(Product, \"\")) %>% dplyr::mutate(GeneID = tidyr::replace_na(GeneID, \"\")) %>%\n")
        rhandle.write("dplyr::arrange(Variant,Annot.group) %>% dplyr::mutate(MostCommonProductOtherGenomes=names(rev(sort(table( Product )))[1]) ) %>%\n")
        rhandle.write("dplyr::arrange(desc(GeneID),MostCommonProductOtherGenomes) %>% arrange(Variant,desc(Product),Evalue) %>%\n")
        rhandle.write("dplyr::select(-Annot.group) %>% dplyr::mutate(MostCommonProductOtherGenomes=gsub(\"(^[[:alpha:]])\", \"\\\\U\\\\1\", MostCommonProductOtherGenomes, perl=TRUE)) %>%\n")
        rhandle.write("dplyr::mutate(Product=gsub(\"(^[[:alpha:]])\", \"\\\\U\\\\1\", Product, perl=TRUE)) %>% dplyr::filter(row_number()==1) %>%\n")
        rhandle.write("dplyr::mutate(MostCommonProductOtherGenomes=ifelse(MostCommonProductOtherGenomes==\"\",\"Intergenic region\",MostCommonProductOtherGenomes)) %>%\n")
        rhandle.write("dplyr::mutate(Product=ifelse(Product==\"\",\"Intergenic region\",Product)) %>% write.table(file=\"Summary."+str(output_rep_file)+"\",sep=\"\t\",row.names=FALSE)")

    rcmd = []
    rcmd.append("Rscript")
    rcmd.append(tmp_rscript_file)
    rcmd = " ".join(rcmd)
    retval = subprocess.call(rcmd,shell=True,stdout=open(os.devnull,'w'),stderr=subprocess.STDOUT)

    cmd = []
    cmd.append("rm -rf")
    cmd.append(tmp_blast_output)
    cmd.append(tmp_genome_names)
    cmd.append(tmp_rscript_file)
    cmd.append(rand_str_val+"*")
    cmd = " ".join(cmd)
    retval = subprocess.call(cmd,shell=True,stdout=open(os.devnull,'w'),stderr=subprocess.STDOUT)

    variantList = ["Variant"]
    for r in SeqIO.parse(seq_kmers,"fasta"):
        variantList.append(str(r.id))

    variantData=[]

    if os.path.exists("Summary."+str(output_rep_file)):
        with open("Summary."+str(output_rep_file),'r') as fhandle:
            for i in fhandle:
                variantData.append(str(i).strip().replace('"','').split('\t'))

        check=False
        with  open("Summary.final."+str(output_rep_file),'w') as fhandle:
            for k in variantList:
                check=False
                for l in variantData:
                    if l[0]==k:
                        check=True
                        break
                if check:
                    fhandle.write('\t'.join(l)+"\n")
                else:
                    fhandle.write(str(k)+str("\t--"*(len(l)-1))+"\n")

    else:
        with  open("Summary.final."+str(output_rep_file),'w') as fhandle:
            headers=["Variant","Genome","MatchLength","GenomeStart","GenomeEnd","GeneStart","GeneEnd","GeneID","LocusTag","Product","Evalue","MostCommonProductOtherGenomes"]
            fhandle.write("\t".join(headers)+"\n")

            for i in SeqIO.parse(seq_kmers,"fasta"):
                fhandle.write(str(i.id)+str("\t--"*(len(headers)-1))+"\n")

if __name__=="__main__":
    main()
