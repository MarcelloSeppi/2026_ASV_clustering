import qiime2
import pandas as pd

import numpy as np
import pandas as pd

import traceback
import os
from os.path import join

from typing import Literal

from qiime2.plugins import rescript
from qiime2.plugins import feature_table
from qiime2.plugins import feature_classifier

from qiime2.plugins import cutadapt
from qiime2.plugins import dada2
from qiime2.plugins import demux
















########################################################
########################################################
################ CLASSIFIER PREPARATION ################






#-------------------------------------------------------------------
## Import gtdb representative sequences (reps) and taxonomy

# import reference reps database with corresponding taxonomy
# extract ASVs (seq,tax)
def import_parse_reps_dataset(seq_path,tax_path,force_download=False):
    #from qiime2.plugins import rescript

    if os.path.exists(seq_path) and os.path.exists(tax_path) and not force_download:
        seq = qiime2.Artifact.load(seq_path)
        tax = qiime2.Artifact.load(tax_path)

    else:
        res = rescript.methods.get_gtdb_data(
            version='202.0',
            db_type='SpeciesReps',
            domain='Bacteria'       # any of ["Bacteria","Archea","Both"]
        )
        seq = res.gtdb_sequences
        tax = res.gtdb_taxonomy
        
        # save
        os.makedirs(os.path.dirname(seq_path),exist_ok=True)
        os.makedirs(os.path.dirname(tax_path),exist_ok=True)
        seq.save(seq_path)      
        tax.save(tax_path)
    
    return seq, tax


## example usage

# paths
#seq_path=join(save_folder,'0_ASVs_seq.qza')
#tax_path=join(save_folder,'0_ASVs_tax.qza')

# run
#seq,tax = import_parse_reps_dataset(seq_path,tax_path,force_download=False)

# reimport
#seq = qiime2.Artifact.load(seq_path)
#tax = qiime2.Artifact.load(tax_path)

#-------------------------------------------------------------------
#-------------------------------------------------------------------
## Upsteram selection of ASVs with UHGG-compatible taxonomies

# match GTDB references with UHGG taxonomy (requires UHGG_genomes-all_metadata.tar.gz file)
def filter_gtdb_reps_with_uhgg_taxonomy(seq,tax,seq_path,tax_path,matching_level="genus",UHGG_genomes_all_metadata_path="./metadata/UHGG_genomes-all_metadata.tar.gz"):
    from qiime2.plugins import feature_table
    
    import numpy as np
    import pandas as pd

    def genus_level_cast(slt):
        # slt := "species_level_taxonomy"
        return slt.str.split(";").str[:-1].str.join(";")

    # import gtdb reps with corresponding taxonomy 
    gtdb_202_bac_tax = tax.view(pd.DataFrame)
    GTDB_spec  = gtdb_202_bac_tax["Taxon"]                    # species level taxonomy (of each reps)
    GTDB_genus = genus_level_cast(gtdb_202_bac_tax["Taxon"])  # genus level taxonomy   (of each reps)

    # import UHGG taxonomy
    UHGGmd     = pd.read_csv(UHGG_genomes_all_metadata_path,sep="\t",index_col=0)
    UHGG_spec  = UHGGmd["Lineage"].drop_duplicates()                   # unique UHGG species level taxonomies
    UHGG_genus = genus_level_cast(UHGGmd["Lineage"]).drop_duplicates() # unique UHGG genus level taxonomies

    # filter reps according to UHGG compatibility, at species or genus level. (keep if reps taxonomy is found in UHGG)
    isin_spec  = GTDB_spec .isin(UHGG_spec )
    isin_genus = GTDB_genus.isin(UHGG_genus)
    filtered_GTDB_byspec  = GTDB_spec.loc[isin_spec ]
    filtered_GTDB_bygenus = GTDB_spec.loc[isin_genus] # keep species-level depth also for genusl-level selection

    #return GTDB_genus,UHGG_genus
    
    # show compatibility metrics
    print({"all_reps":len(gtdb_202_bac_tax),"spec_compatible_reps":len(filtered_GTDB_byspec),"genus_compatible_reps":len(filtered_GTDB_bygenus)})

    # make qiime metadata object for filtering
    assert matching_level in ["genus","species"]
    if matching_level=="genus":
        filtered_GTDB = filtered_GTDB_bygenus
    if matching_level=="species":
        filtered_GTDB = filtered_GTDB_byspec
    
    #display(filtered_GTDB_bygenus)
    #display(filtered_GTDB_byspec)
    metadata = qiime2.Metadata(pd.DataFrame(data=None,index=filtered_GTDB.index))

    # filter sequences
    filtered_seq = feature_table.methods.filter_seqs(
        data=seq,
        metadata=metadata
    ).filtered_data
    filtered_tax = qiime2.Artifact.import_data('FeatureData[Taxonomy]',filtered_GTDB)
    
    filtered_seq.save(seq_path)
    filtered_tax.save(tax_path)

    return filtered_seq,filtered_tax


## example usage

# paths
#seq_path=join(save_folder,'1_ASVs_seq-uhgg_fitered.qza')
#tax_path=join(save_folder,'1_ASVs_tax-uhgg_fitered.qza')

# run
#seq,tax = filter_gtdb_reps_with_uhgg_taxonomy(seq,tax,seq_path,tax_path,matching_level="species")

# reimport
#seq = qiime2.Artifact.load(seq_path)
#tax = qiime2.Artifact.load(tax_path)

#-------------------------------------------------------------------



#-------------------------------------------------------------------
## Culling low-quality reps

def cull_sequences(seq,tax,seq_path,tax_path):
    #from qiime2.plugins import rescript

    culled = rescript.methods.cull_seqs(sequences=seq)
    culled_seq = culled.clean_sequences
    culled_seq.save(seq_path)
    tax       .save(tax_path)
    return culled_seq,tax

## example usage

# paths
#seq_path=join(save_folder,'2_ASVs_seq-full_seq_culled.qza')
#tax_path=join(save_folder,'2_ASVs_tax-full_seq_culled.qza')

# run
#seq,tax = cull_sequences(seq,tax,seq_path,tax_path)

# reimport
#seq = qiime2.Artifact.load(seq_path)
#tax = qiime2.Artifact.load(tax_path)


#-------------------------------------------------------------------
#-------------------------------------------------------------------
## Dereplicate full-length reps

def dereplicate(seq,tax,seq_path,tax_path):
    #from qiime2.plugins import rescript

    derep = rescript.methods.dereplicate(sequences=seq,taxa=tax,mode='lca')
    drep_seq = derep.dereplicated_sequences
    drep_tax = derep.dereplicated_taxa
    drep_seq.save(seq_path)
    drep_tax.save(tax_path)
    return drep_seq,drep_tax

## example usage

# paths
#seq_path=join(save_folder,'3_ASVs_seq-dereplicated.qza')
#tax_path=join(save_folder,'3_ASVs_tax-dereplicated.qza')

# run
#seq,tax = dereplicate(seq,tax,seq_path,tax_path)

# reimport
#seq = qiime2.Artifact.load(seq_path)
#tax = qiime2.Artifact.load(tax_path)

#-------------------------------------------------------------------
#-------------------------------------------------------------------
## Extract amplicon region from reps

def extract_amplicon_region(seq,tax,seq_path,tax_path,f_primer,r_primer):
    #from qiime2.plugins import feature_classifier

    amplicon = feature_classifier.methods.extract_reads(
                    sequences=seq,
                    f_primer=f_primer,
                    r_primer=r_primer,
                    read_orientation='forward',
                    n_jobs=3
                )
    amp_seq = amplicon.reads
    amp_seq.save(seq_path)
    tax    .save(tax_path)
    return amp_seq,tax

# example usage

# paths
#seq_path=join(save_folder,'4_ASVs_seq-extracted_amplicons.qza')
#tax_path=join(save_folder,'4_ASVs_tax-extracted_amplicons.qza')

# run
#f_primer='GTGCCAGCMGCCGCGGTAA'
#r_primer='GGACTACHVGGGTWTCTAAT'
#seq,tax = extract_amplicon_region(seq,tax,seq_path,tax_path,f_primer,r_primer)

# reimport
#seq = qiime2.Artifact.load(seq_path)
#tax = qiime2.Artifact.load(tax_path)

#-------------------------------------------------------------------
#-------------------------------------------------------------------
## Fit classifier over reps

def fit_classifier(seq,tax,classifier_path):
    #from qiime2.plugins import feature_classifier

    classifier = feature_classifier.methods.fit_classifier_naive_bayes(
                    reference_reads=seq,
                    reference_taxonomy=tax
                )
    classifier.classifier.save(classifier_path)
    return classifier.classifier

## example usage

# paths
#classifier_path=join(save_folder,'6_ASVs_classifier.qza')

# run
#classifier = fit_classifier(seq,tax,classifier_path)

# reimport
#classifier = qiime2.Artifact.load(classifier_path)

#-------------------------------------------------------------------
#-------------------------------------------------------------------





############################ FULL CLASSIFIER PREPARATION #########################

#-------------------------------------------------------------------
# full classification pipeline
def prepare_GTDB_classifier(raw_refs_folder, save_folder,UHGG_matching_level,f_primer,r_primer,UHGG_metadata,start_from=0):

    #### **Import/download gtdb representative sequences (reps) and taxonomy**
    
    # step index
    I = 0

    # add I-step intermediates paths to list
    seq_path = [join(save_folder,'0.1_ASVs_seq-uhgg_fitered.qza')]
    tax_path = [join(save_folder,'0.1_ASVs_tax-uhgg_fitered.qza')]

    # run step if after starting position
    if start_from <= I:
        print("\n0) Import/download gtdb representative sequences (reps) and taxonomy")

        # make folder in which classifier's intermediate and final results will be stored
        os.makedirs(save_folder,exist_ok=True)

        # paths
        seq_path_gtdb=join(raw_refs_folder,'0_ASVs_seq.qza')
        tax_path_gtdb=join(raw_refs_folder,'0_ASVs_tax.qza')
        
        # run
        seq,tax = import_parse_reps_dataset(seq_path_gtdb,tax_path_gtdb,force_download=False)


        ##### **Upsteram selection of ASVs with UHGG-compatible taxonomies**
        print("0.1) Upsteram selection of ASVs with UHGG-compatible taxonomies")

        if UHGG_matching_level is not None:

            # run
            seq,tax = filter_gtdb_reps_with_uhgg_taxonomy(seq,tax,seq_path[I],tax_path[I],matching_level=UHGG_matching_level,UHGG_genomes_all_metadata_path=UHGG_metadata)
    


    #### **Culling low-quality reps**
    
    # step index
    I = 1 

    # add I-step intermediates paths to list
    seq_path.append(join(save_folder,'1_ASVs_seq-full_seq_culled.qza'))
    tax_path.append(join(save_folder,'1_ASVs_tax-full_seq_culled.qza'))

    # run step if after starting position
    if start_from <= I:
        print("\n1) Culling low-quality reps")

        # load intermediates if starting step
        if start_from == I: 
            seq = qiime2.Artifact.load(seq_path[I-1])
            tax = qiime2.Artifact.load(tax_path[I-1])
    
        # run
        seq,tax = cull_sequences(seq,tax,seq_path[I],tax_path[I])


    #### **Dereplicate full-length reps**

    # step index
    I = 2

    # add I-step intermediates paths to list
    seq_path.append(join(save_folder,'2_ASVs_seq-dereplicated.qza'))
    tax_path.append(join(save_folder,'2_ASVs_tax-dereplicated.qza'))

    # run step if after starting position
    if start_from <= I:
        print("\n2) Dereplicate full-length reps")

        # load intermediates if starting step
        if start_from == I: 
            seq = qiime2.Artifact.load(seq_path[I-1])
            tax = qiime2.Artifact.load(tax_path[I-1])
    
        # run
        seq,tax = dereplicate(seq,tax,seq_path[I],tax_path[I])



    #### **Extract amplicon region from reps**

    # step index
    I = 3 

    # add I-step intermediates paths to list
    seq_path.append(join(save_folder,'3_ASVs_seq-extracted_amplicons.qza'))
    tax_path.append(join(save_folder,'3_ASVs_tax-extracted_amplicons.qza'))

    # run step if after starting position
    if start_from <= I:
        print("\n3) Extract amplicon region from reps")

        # load intermediates if starting step
        if start_from == I: 
            seq = qiime2.Artifact.load(seq_path[I-1])
            tax = qiime2.Artifact.load(tax_path[I-1])
    
        # run
        seq,tax = extract_amplicon_region(seq,tax,seq_path[I],tax_path[I],f_primer,r_primer)

  

    #### **Dereplicate again amplicon-regon reps**

    # step index
    I = 4 

    # add I-step intermediates paths to list
    seq_path.append(join(save_folder,'4_ASVs_seq-amplicons_dereplicated.qza'))
    tax_path.append(join(save_folder,'4_ASVs_tax-amplicons_dereplicated.qza'))

    # run step if after starting position
    if start_from <= I:
        print("\n4) Dereplicate again amplicon-regon reps")

        # load intermediates if starting step
        if start_from == I: 
            seq = qiime2.Artifact.load(seq_path[I-1])
            tax = qiime2.Artifact.load(tax_path[I-1])
    
        # run
        seq,tax = dereplicate(seq,tax,seq_path[I],tax_path[I])



    #### **Fit classifier over reps**

    # step index
    I = 5 

    # add I-step intermediates paths to list
    classifier_path=join(save_folder,'5_ASVs_classifier.qza')

    # run step if after starting position
    if start_from <= I:
        print("\n5) Fit classifier over reps")

        # load intermediates if starting step
        if start_from == I: 
            seq = qiime2.Artifact.load(seq_path[I-1])
            tax = qiime2.Artifact.load(tax_path[I-1])
    
        # run
        classifier = fit_classifier(seq,tax,classifier_path)



    #### END
    return classifier

## example usage

# import primers
#primers_table = pd.read_csv("../primers_table.tsv",sep="\t",index_col=0)
#UHGG_genomes_all_metadata_path = "./metadata/UHGG_genomes-all_metadata.tar.gz"

# classifier selection
#raw_refs_folder     = "GTDB_raw_refs"

#for UHGG_matching_level in ["species","genus",None]:
#    for primers_set in ["515f_806r_A","515f_806r_B"]:
#        save_folder = f"GTDB_202_bac/UHGG_2.0.1_matching_{UHGG_matching_level}/primers_set_{primers_set}/"
#        print(f">>>> {save_folder} <<<<")
#        # run
#        if UHGG_matching_level is None:
#            classifier = prepare_GTDB_classifier(raw_refs_folder,
#                                             save_folder,
#                                             UHGG_matching_level,
#                                             primers_table["f_primer"].loc[primers_set],
#                                             primers_table["r_primer"].loc[primers_set],
#                                             UHGG_genomes_all_metadata_path
#                                             )
#-------------------------------------------------------------------

























########################################################
########################################################
################## DATA PREPROCESSING ##################





#************* Utility functions ***************

# clear qiime folder (deletes all files and directoris in it)
# preserves manifest.tsv file
def clear_preprocessing_folder(path):
    import os
    import shutil

    try:
        # Check if the path exists
        if os.path.exists(path) and os.path.isdir(path):
            # Iterate over each item in the directory
            for item in os.listdir(path):
                if item != "manifest.tsv":
                    item_path = os.path.join(path, item)
                    if os.path.isdir(item_path):
                        # If the item is a folder, delete it and all its contents
                        shutil.rmtree(item_path)
                    else:
                        # If it's a file, delete it
                        os.remove(item_path)
            print(f"Folder '{path}' has been cleared.")
        else:
            print(f"Path '{path}' does not exist or is not a directory.")
    except Exception as e:
        print(f"An error occurred: {e}")



# given a folder containing runs and a metadata path
# finds all runs and divides them by primer
# if batch_size!=None further divides in batches
# creates directory tree, main_data_folder/primers_name/batch/manifest.tsv
def compile_manifest_files(raw_runs_folder,preprocessing_dir,metadata_file_path,additional_groupby_fields=["instrument_model"],batch_size=None,full_preprocessing_path_cleaning=False,preprocessing_type:Literal["by_group","by_run"]="by_group"):
    import os
    
    assert preprocessing_type in ["by_group","by_run"]
    
    # assemble "manifest file" with run_id, absolute_path and orientation/direction

    # all runs in folder
    run_acc = pd.Index([ f.split("/")[-1].split(".")[0]    for f in os.listdir(raw_runs_folder) if f.endswith(".fastq.gz")])
    
    # absolute paths to run files
    rel_paths = [join(raw_runs_folder,f"{ra}.fastq.gz") for ra in run_acc] # standard relative path within the project
    abs_paths = [os.path.abspath(rp) for rp in rel_paths]
    
    # orientation of the reads (defaul values)
    orientation = ["forward"]*len(run_acc)

    # "sampleid", "absolute-filepath", and "direction" are expected qiime keywords (additional fields are ignored by qiime)
    # "pcr_primers" and "pcr_primers_name" are additional fields taken from metadata but necessary for data processing
    # the dada2 denoising is performed over groups homogeneous with respect to primers and additional fields
    meta = pd.read_parquet(metadata_file_path).set_index("run_accession").loc[run_acc]
    manifest = pd.DataFrame({"sampleid"          : run_acc,
                             "absolute-filepath" : abs_paths,
                             "direction"         : orientation,
                             "relative-filepath" : rel_paths,
                             "pcr_primers_name"  : meta["pcr_primers_name"],
                             "f_primer"          : meta["pcr_primers"].str.split(";").str[0].str.split(":").str[-1].str.strip(),
                             "r_primer"          : meta["pcr_primers"].str.split(";").str[1].str.split(":").str[-1].str.strip(),
                             }).set_index("sampleid")
    for field in additional_groupby_fields:
        manifest[field] = meta[field]
    

    if preprocessing_type=="by_run":
        # complete target directory
        preprocessing_dir = join(preprocessing_dir,"by_run")
        os.makedirs(preprocessing_dir,exist_ok=True)

        # clear target directory
        if full_preprocessing_path_cleaning:
            clear_preprocessing_folder(preprocessing_dir)
        

        manifest_file_paths = pd.Series(index=pd.Index([],name=manifest.index.name))
        for run in manifest.index:
            # prepare folder
            folder_name = join(preprocessing_dir,run)
            os.makedirs(folder_name,exist_ok=True)

            manifest_file_path = join(folder_name,"manifest.tsv")
            manifest.loc[[run]].to_csv(manifest_file_path,sep="\t")
            manifest_file_paths.loc[run] = manifest_file_path
        
    elif preprocessing_type=="by_group":
        # complete target directory
        preprocessing_dir = join(preprocessing_dir,"by_group")
        os.makedirs(preprocessing_dir,exist_ok=True)
        
        # clear target directory
        if full_preprocessing_path_cleaning:
            clear_preprocessing_folder(preprocessing_dir)
        
        
        groupby_fields = np.concatenate([["pcr_primers_name"]+additional_groupby_fields])

        # select group with homegeneous primers
        def make_manifest_file_path(manifest_group):
        
            # extract group index (by construction, same for all rows)
            group_index = manifest_group.iloc[0].loc[groupby_fields]

            # manage batching (fix maximum batch length)
            if batch_size is not None:
                L = len(manifest_group)
                N_batch = int(np.ceil(L/float(batch_size)))
                batches = np.array_split(manifest_group.index,N_batch) 
            else:
                batches = [manifest_group.index]

            manifest_file_paths = pd.Series(index=pd.Index([],name="batch"))
            for i,batch in enumerate(batches):
                # prepare folder
                folder_name = join(preprocessing_dir,*group_index,f"batch_{i}")
                os.makedirs(folder_name,exist_ok=True)

                manifest_file_path = join(folder_name,"manifest.tsv")
                manifest.loc[batch].to_csv(manifest_file_path,sep="\t")

                manifest_file_paths.loc[i] = manifest_file_path

            return manifest_file_paths

        # for each homogeneous group creates the manifest file and the corresponding filepath
        manifest_file_paths = manifest.groupby(by=groupby_fields.tolist(),group_keys=True).apply(make_manifest_file_path)
                
    return manifest_file_paths



# updates the absolute path to raw run writte in manifest files
# (necessary when proccessing on new machine environment)
# (safer than recomputing the manifest files: no clear folder option)
def update_runs_absolute_paths(manifest_file_paths):
    for mfp in manifest_file_paths:
        manifest = pd.read_csv(mfp,sep="\t")
        manifest["absolute-filepath"] = manifest["relative-filepath"].apply(lambda rel_path:os.path.abspath(rel_path))
        manifest.to_csv(mfp,sep="\t",index=False)
    return None



#*************** Preprocessing Pipeline **************

#--------------------------------------------------------------
## 0) import data into qiime format
def qiime2_import_raw_fastq(manifest_file_path):
    seq_raw = qiime2.Artifact.import_data('SampleData[SequencesWithQuality]',
                                  manifest_file_path,
                                  view_type='SingleEndFastqManifestPhred33V2')
    return seq_raw
#---------------------------------------------------------------
#---------------------------------------------------------------
## 1) remove primers
def remove_primers(base_qiime_dir, seq_raw, cores=3):

    # select primers
    f_primer, r_primer = pd.read_csv(join(base_qiime_dir,"manifest.tsv"),sep="\t").iloc[0].loc[["f_primer","r_primer"]]
    
    # remove primers from seq
    remove_primers_r = cutadapt.methods.trim_single(
        demultiplexed_sequences=seq_raw,
        front=[f_primer],adapter=[r_primer],
        cores=cores
    )
    seq_noprimers = remove_primers_r.trimmed_sequences
    return seq_noprimers
#---------------------------------------------------------------
#---------------------------------------------------------------
## 2) quality trimming

# choose truncation length
# from visualize quality trends across reads length
# rule of thumb --> quality falls ~below 30 --> truncate
def suggest_trunc_from_visualization(demux_summary_path, q_threshold=30):
    from zipfile import ZipFile
    
    # demux_summary: qiime2.sdk.visualization.Visualization
    with ZipFile(demux_summary_path, 'r') as z:
        quality_filepath = [name for name in z.namelist() if name.endswith("/data/forward-seven-number-summaries.tsv")][0]

        # find per-base quality file (inspect names if needed)
        with z.open(quality_filepath) as f:
            quality_quantiles = pd.read_csv(f, sep='\t', index_col=0)

    medians = quality_quantiles.T["50%"].to_numpy()

    # find first position where median < threshold
    below = np.where(medians < q_threshold)[0]

    if below.size == 0:
        return len(medians)  # quality never drops below threshold -> keep full length
    else:
        # trunc_len should be position index (1-based read length) where quality still >= threshold
        trunc_len = int(below[0])  # zero-based -> trunc_len equals first index below threshold
        if trunc_len == 0:
            # very low quality at start — better to set trunc_len small or re-evaluate
            return 0
        else:
            return trunc_len

# compute visualization object and return truncation length suggestion
def  compute_truncation_length_from_quality_scores(seq,qiime_dir,q_threshold=30):
    
    # compute visualization object
    demux_summary = demux.visualizers.summarize(seq)
    demux_summary_path = join(qiime_dir,"2_quality_summary_visualization.qzv")
    demux_summary.visualization.save(demux_summary_path)
    
    # compute truncation length
    truncation_length = suggest_trunc_from_visualization(demux_summary_path,q_threshold)

    return truncation_length
#---------------------------------------------------------------
#---------------------------------------------------------------
## 3) denoise and ASVs ectraction
def denoise_ASVs_extraction(seq,truncation_length,cores=3):

    if truncation_length==0:
        raise ValueError(f"{truncation_length=}, therefore there is no sequence left to be processed")

    # denoise
    denoise_r = dada2.methods.denoise_single(
        demultiplexed_seqs  =seq,
        trunc_len=truncation_length, # choose from quality plot
        trim_left=0,                 # if you need to remove extra bases
        n_threads=cores
    )

    # extract results
    asv_tab         = denoise_r.table                   
    asv_seq         = denoise_r.representative_sequences
    denoising_stats = denoise_r.denoising_stats         

    return asv_tab, asv_seq, denoising_stats
#---------------------------------------------------------------







def create_patched_r_script(output_path):
    """
    Creates a patched version of run_dada.R that exposes OMEGA_A as a CLI argument.
    The original script is not modified.
    
    Parameters
    ----------
    output_path : str
        Path where the patched R script will be saved.
    """
    import shutil
    import os

    r_script_path = shutil.which('run_dada.R')

    # safety check: make sure we're not reading an already-patched file
    with open(r_script_path, 'r') as f:
        content = f.read()

    if r_script_path is None:
        raise FileNotFoundError("run_dada.R not found in PATH")

    if 'omega_a' in content:
        raise RuntimeError(
            f"The original R script at {r_script_path} already contains 'omega_a' — "
            "it may have been patched on disk previously. "
            "Please restore the original file before running this function."
        )


    patches = [
        # 1) Add omega_a to optparse option list
        (
            '''  make_option(c("--band_size"), action="store", default='NULL', type='character',
              help="When set, banded Needleman-Wunsch alignments are performed.")
)''',
            '''  make_option(c("--band_size"), action="store", default='NULL', type='character',
              help="When set, banded Needleman-Wunsch alignments are performed."),
  make_option(c("--omega_a"), action="store", default='NULL', type='character',
              help="Threshold for calling new ASVs. Lower = more conservative. Default 1e-40.")
)'''
        ),
        # 2) Add OMEGA_A parsing after BAND_SIZE parsing
        (
            "BAND_SIZE <- if(opt$band_size=='NULL') NULL else as.integer(opt$band_size)",
            """BAND_SIZE <- if(opt$band_size=='NULL') NULL else as.integer(opt$band_size)
OMEGA_A <- if(opt$omega_a=='NULL') 1e-40 else as.numeric(opt$omega_a)"""
        ),
        # 3) Single-end first pass
        (
            "dds[[j]] <- dada(drp, err=err, multithread=multithread,HOMOPOLYMER_GAP_PENALTY=HOMOPOLYMER_GAP_PENALTY,\n                     BAND_SIZE=BAND_SIZE, verbose=FALSE)",
            "dds[[j]] <- dada(drp, err=err, multithread=multithread,HOMOPOLYMER_GAP_PENALTY=HOMOPOLYMER_GAP_PENALTY,\n                     BAND_SIZE=BAND_SIZE, OMEGA_A=OMEGA_A, verbose=FALSE)"
        ),
        # 4) Single-end pseudo-pool second pass
        (
            "priors=pseudo_priors, HOMOPOLYMER_GAP_PENALTY=HOMOPOLYMER_GAP_PENALTY,\n                       BAND_SIZE=BAND_SIZE, verbose=FALSE)",
            "priors=pseudo_priors, HOMOPOLYMER_GAP_PENALTY=HOMOPOLYMER_GAP_PENALTY,\n                       BAND_SIZE=BAND_SIZE, OMEGA_A=OMEGA_A, verbose=FALSE)"
        ),
        # 5) Paired-end first pass (F and R)
        (
            "ddsF[[j]] <- dada(drpF, err=err, multithread=multithread, verbose=FALSE)\n    drpR <- derepFastq(filtsR[[j]])\n    ddsR[[j]] <- dada(drpR, err=errR, multithread=multithread, verbose=FALSE)",
            "ddsF[[j]] <- dada(drpF, err=err, multithread=multithread, OMEGA_A=OMEGA_A, verbose=FALSE)\n    drpR <- derepFastq(filtsR[[j]])\n    ddsR[[j]] <- dada(drpR, err=errR, multithread=multithread, OMEGA_A=OMEGA_A, verbose=FALSE)"
        ),
        # 6) Paired-end pseudo-pool second pass (F and R)
        (
            "ddsF[[j]] <- dada(drpF, err=err, priors=pseudo_priorsF,\n                        multithread=multithread, verbose=FALSE)\n      drpR <- derepFastq(filtsR[[j]])\n      ddsR[[j]] <- dada(drpR, err=errR, priors=pseudo_priorsR,\n                        multithread=multithread, verbose=FALSE)",
            "ddsF[[j]] <- dada(drpF, err=err, priors=pseudo_priorsF,\n                        multithread=multithread, OMEGA_A=OMEGA_A, verbose=FALSE)\n      drpR <- derepFastq(filtsR[[j]])\n      ddsR[[j]] <- dada(drpR, err=errR, priors=pseudo_priorsR,\n                        multithread=multithread, OMEGA_A=OMEGA_A, verbose=FALSE)"
        ),
    ]

    for old, new in patches:
        if old not in content:
            raise ValueError(f"Could not find the following string to patch:\n{old}")
        content = content.replace(old, new)

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as f:
        f.write(content)

    print(f"Patched R script saved to: {output_path}")

    # verify all OMEGA_A occurrences
    count = content.count('OMEGA_A')
    print(f"OMEGA_A appears {count} times in patched script (expected 9)")

    # set as executable
    import stat
    os.chmod(output_path, os.stat(output_path).st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)



def apply_omega_a_patch(patched_r_script_path):
    """
    Monkey-patches _denoise_single() to use the patched R script and
    accept an omega_a parameter.

    Parameters
    ----------
    patched_r_script_path : str
        Path to the patched run_dada.R script created by create_patched_r_script().
    """
    import os
    import tempfile
    import subprocess
    import q2_dada2._denoise as _denoise

    if not os.path.exists(patched_r_script_path):
        raise FileNotFoundError(f"Patched R script not found at: {patched_r_script_path}")

    # check if patch has already been applied
    if getattr(_denoise._denoise_single, '_is_omega_a_patched', False):
        print("Patch already applied, skipping.")
        return


    # keep a reference to the original functions we still need
    _check_inputs  = _denoise._check_inputs
    _denoise_helper = _denoise._denoise_helper
    run_commands   = _denoise.run_commands

    def _patched_denoise_single(demultiplexed_seqs, trunc_len, trim_left, max_ee, trunc_q,
                                max_len, pooling_method, chimera_method,
                                min_fold_parent_over_abundance, allow_one_off,
                                n_threads, n_reads_learn, hashed_feature_ids,
                                homopolymer_gap_penalty, band_size, retain_all_samples,
                                omega_a=1e-40):  # <-- new parameter

        _valid_keys = set(_denoise._valid_inputs.keys())
        _check_inputs(**{k: v for k, v in locals().items() if k in _valid_keys})

        if trunc_len != 0 and trim_left >= trunc_len:
            raise ValueError("trim_left (%r) must be smaller than trunc_len (%r)"
                             % (trim_left, trunc_len))
        if max_len != 0 and max_len < trunc_len:
            raise ValueError("trunc_len (%r) must be no bigger than max_len (%r)"
                             % (trunc_len, max_len))
        max_len = 'Inf' if max_len == 0 else max_len

        with tempfile.TemporaryDirectory() as temp_dir_name:
            biom_fp = os.path.join(temp_dir_name, 'output.tsv.biom')
            track_fp = os.path.join(temp_dir_name, 'track.tsv')

            cmd = [patched_r_script_path,          # <-- patched script
                   '--input_directory', str(demultiplexed_seqs),
                   '--output_path', biom_fp,
                   '--output_track', track_fp,
                   '--filtered_directory', temp_dir_name,
                   '--truncation_length', str(trunc_len),
                   '--trim_left', str(trim_left),
                   '--max_expected_errors', str(max_ee),
                   '--truncation_quality_score', str(trunc_q),
                   '--max_length', str(max_len),
                   '--pooling_method', str(pooling_method),
                   '--chimera_method', str(chimera_method),
                   '--min_parental_fold', str(min_fold_parent_over_abundance),
                   '--allow_one_off', str(allow_one_off),
                   '--num_threads', str(n_threads),
                   '--learn_min_reads', str(n_reads_learn),
                   '--homopolymer_gap_penalty', str(homopolymer_gap_penalty),
                   '--band_size', str(band_size),
                   '--omega_a', str(omega_a)]        # <-- new argument

            try:
                run_commands([cmd])
            except subprocess.CalledProcessError as e:
                if e.returncode == 2:
                    raise ValueError(
                        "No reads passed the filter. trunc_len (%r) may be longer"
                        " than read lengths, or other arguments (such as max_ee"
                        " or trunc_q) may be preventing reads from passing the"
                        " filter." % trunc_len)
                else:
                    raise Exception("An error was encountered while running DADA2"
                                    " in R (return code %d), please inspect stdout"
                                    " and stderr to learn more." % e.returncode)

            return _denoise_helper(biom_fp, track_fp, hashed_feature_ids,
                                   retain_all_samples)

    # mark the patched function so we can detect it
    _patched_denoise_single._is_omega_a_patched = True

    _denoise._denoise_single = _patched_denoise_single
    print("Patch applied: _denoise_single now accepts omega_a and uses patched R script")



PATCHED_R_SCRIPT_PATH = './patched_dada2/run_dada_patched.R'

def denoise_ASVs_extraction2(seq, truncation_length, create_patch=False, omega_a=1e-40, cores=3):
    import q2_dada2._denoise as _denoise
    from q2_dada2._denoise import SingleLanePerSampleSingleEndFastqDirFmt
    
    # --- patch setup ---
    if create_patch:
        create_patched_r_script(PATCHED_R_SCRIPT_PATH)
    apply_omega_a_patch(PATCHED_R_SCRIPT_PATH)

    # --- validate ---
    if truncation_length == 0:
        raise ValueError(f"{truncation_length=}, therefore there is no sequence left to be processed")

    # --- unwrap QIIME2 artifact ---
    seq_unwrapped = seq.view(SingleLanePerSampleSingleEndFastqDirFmt)

    # --- denoise ---
    table, rep_sequences, metadata = _denoise._denoise_single(
        demultiplexed_seqs             = seq_unwrapped,
        trunc_len                      = truncation_length,
        trim_left                      = 0,
        max_ee                         = 2.0,
        trunc_q                        = 2,
        max_len                        = 0,
        pooling_method                 = 'independent',
        chimera_method                 = 'consensus',
        min_fold_parent_over_abundance = 1.0,
        allow_one_off                  = False,
        n_threads                      = cores,
        n_reads_learn                  = 1000000,
        hashed_feature_ids             = True,
        homopolymer_gap_penalty        = 'NULL',
        band_size                      = '16',
        retain_all_samples             = True,
        omega_a                        = omega_a
    )

    return table, rep_sequences, metadata






################### full_preprocessing_pipeline #######################
#---------------------------------------------------------------
def full_preprocessing_pipeline(manifest_file_path,cores=3,q_threshold=30,from_step=0,to_step=3,clear_path=False,overwrite=True,skip_failed=False):
    
    # extract pipeline dir from manifest file path
    qiime_dir = os.path.dirname(manifest_file_path)
    
    ## PREAMBLE
    print(f"\n\n>>>> processing {manifest_file_path} <<<")
    print(f"#samples = {len(pd.read_csv(manifest_file_path))}")
    Exit1 = False
    Exit2 = False

    # check if there is on disk a computed result for the last pipeline step
    computed_last_step = os.path.exists(join(qiime_dir,'3_asv_seq.qza')) and os.path.exists(join(qiime_dir,'3_asv_tab.qza')) and os.path.exists(join(qiime_dir,'3_dada2_stats.qza'))
    if computed_last_step:
        print("alredy exsists...")
        if not overwrite:
            print("do not recompute...")
            Exit1=True
        else:
            print("recompute...")
            
    
    # check if there is a "FAILED" fail indicating a failed run
    is_FAILED = os.path.exists(join(qiime_dir,"FAILED"))
    if is_FAILED:
        print("flagget as FAILED...")
        if skip_failed:
            print("do not recompute...")
            Exit2=True
        else:
            print("recompute...")
        
    
    ## COMPUTATION
    if Exit1 or Exit2:
        return None
    
    # compute preprecessing for assigned group
    else:
        try:  
            
            # clear path from all files except manifest
            if clear_path:
                clear_preprocessing_folder(qiime_dir)
  

            ## step 0 import in qiime
            I = 0
            if from_step<=I and to_step>=I:
                print("\n\n-----------------------------------")
                print("\n0) import in qiime")
                # run
                seq = qiime2_import_raw_fastq(manifest_file_path)
                # save
                seq.save(join(qiime_dir,'0_seq_raw.qza'))


            ## step 1 remove primers
            I = 1
            if from_step<=I and to_step>=I:
                print("\n\n-----------------------------------")
                print("\n1) remove primers")
                # if first step, load results from previous steps
                if from_step==I:
                    seq = qiime2.Artifact.load(join(qiime_dir,'0_seq_raw.qza'))
                # run
                remove_primers(qiime_dir,seq,cores=cores)
                # save
                seq.save(join(qiime_dir,"1_seq_noprimers.qza") )


            ## step 2 select truncation length
            I = 2
            if from_step<=I and to_step>=I:
                print("\n\n-----------------------------------")
                print("\n2) compute truncation length")
                # if first step, load results from previous steps
                if from_step==I:
                    seq = qiime2.Artifact.load(join(qiime_dir,'1_seq_noprimers.qza'))
                # run
                truncation_length = compute_truncation_length_from_quality_scores(seq,qiime_dir,q_threshold)
                # save
                with open(join(qiime_dir,"2_truncation_length.txt"),mode="w") as f:
                    f.write(f"{truncation_length=}")

            
            ## step 3 denoise and ASVs extraction
            I = 3
            if from_step<=I and to_step>=I:
                print("\n\n-----------------------------------")
                print("\n3) denoise and ASVs ectraction")
                # if first step, load results from previous steps
                if from_step==I:
                    seq = qiime2.Artifact.load(join(qiime_dir,'1_seq_noprimers.qza'))
                    with open(join(qiime_dir,"2_truncation_length.txt"),mode="r") as f:
                        truncation_length=int(f.readline().split("=")[-1])
                # run
                asv_tab, asv_seq, denoising_stats = denoise_ASVs_extraction(seq,truncation_length,cores=cores)
                # save
                asv_tab        .save(join(qiime_dir,'3_asv_tab.qza'    ))  # feature table (ASV counts per sample)
                asv_seq        .save(join(qiime_dir,'3_asv_seq.qza'    ))  # ASV sequences
                denoising_stats.save(join(qiime_dir,'3_dada2_stats.qza'))  # stats


            # remove "FAILED" file if present (only when previously failed run is completed)
            FAILED = join(qiime_dir,"FAILED")
            if os.path.exists(FAILED):
                os.remove(FAILED)
            
        except Exception as e:
            print("\n\n######RUN FAILED... continue#######\n\n")
            dir = os.path.dirname(manifest_file_path)
            with open(join(dir,"FAILED"),mode="w") as file:
                error_message = traceback.format_exc()
                file.write(f"{e}")
                file.write(error_message)
        return None



























