#include "lpo.h"
#include "msa_format.h"
#include "align_score.h"
#include <pthread.h>


int nseq=0;
LPOSequence_T *seq=NULL,*seqRef=NULL,*seqUnco=NULL;
ResidueScoreMatrix_T score_matrix; /* DEFAULT GAP PENALTIES*/
FILE* seq_ifile=NULL;
char* comment=NULL;
int ibundle=ALL_BUNDLES;
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
int max_input_seqs = 4;
int* sizeArray;
int** seqArray;
int use_aggressive_fusion=0;
int do_progressive=1;
char* pair_score_file=NULL;
int do_global=1;
int do_preserve_sequence_order=1;

static LPOSequence_T *read_partial_order_file (char *po_filename, char *subset_filename, int remove_listed_seqs, int keep_all_links, int do_switch_case, ResidueScoreMatrix_T *mat);


void buildAndAnalysePOMSA (int n_input_seqs, LPOSequence_T* lpo_out, LPOSequence_T **input_seqs, ResidueScoreMatrix_T score_matrix, FILE* seq_ifile);

void freeMem(LPOSequence_T **input_seqs, int nseq, LPOSequence_T *seq, int n_input_seqs);

typedef struct s_ThreadsParams {
	int nb;
	LPOSequence_T* lpo_out;
	LPOSequence_T** in_seqs;
} *ThreadsParams;

void* doThreadJob(void* params) {
	ThreadsParams p = (ThreadsParams) params;
	int tnb = p->nb;
	LPOSequence_T** input_seqs = p->in_seqs;
	LPOSequence_T* lpo_out = p->lpo_out;
	int n_input_seqs;
	int i;
	for (i = 0; i < sizeArray[tnb]; i++) {
		n_input_seqs = 0;
		input_seqs[n_input_seqs++] = &(seqRef[seqArray[tnb][i]]);
		input_seqs[n_input_seqs++] = &(seqUnco[seqArray[tnb][i]]);
		input_seqs[n_input_seqs++] = &(seq[seqArray[tnb][i]]);
		initialize_seqs_as_lpo(1,&(seqRef[seqArray[tnb][i]]),&score_matrix);
		initialize_seqs_as_lpo(1,&(seqUnco[seqArray[tnb][i]]),&score_matrix);
		initialize_seqs_as_lpo(1,&(seq[seqArray[tnb][i]]),&score_matrix);
		buildAndAnalysePOMSA (n_input_seqs, lpo_out, input_seqs, score_matrix, seq_ifile);
		n_input_seqs++;
	}
	
	for (i = 0; i < n_input_seqs; i++) {
		free_lpo_sequence(p->in_seqs[i],0);
	}
	FREE(p->in_seqs);
	free(lpo_out);
	free(params);
}

int main(int argc,char *argv[])
{
  int i,j,nframe_seq=0,use_reverse_complement=0;
  int do_switch_case=dont_switch_case,do_analyze_bundles=0;
  int nseq_in_list=0,n_input_seqs=0;
  char score_file[256],seq_file[256],po_list_entry_filename[256],*al_name="test align";
  LPOSequence_T *lpo_out=NULL,*frame_seq=NULL,*dna_lpo=NULL,*lpo_in=NULL;
  FILE *errfile=stderr,*logfile=NULL,*lpo_file_out=NULL,*po_list_file=NULL,*seqCorrected_ifile=NULL,*seqUncorrected_ifile=NULL,*seqReference_ifile=NULL;
  char *print_matrix_letters=NULL,*fasta_out=NULL,*po_out=NULL,*matrix_filename="./src/poa-graph/blosum80.mat",
    *seq_filename=NULL, *unco_seq_filename=NULL,*ref_seq_filename=NULL,*frame_dna_filename=NULL,*po_filename=NULL,*po2_filename=NULL,
    *po_list_filename=NULL, *hbmin=NULL,*numeric_data=NULL,*numeric_data_name="Nmiscall",
    *dna_to_aa=NULL,*aafreq_file=NULL,*termval_file=NULL,
    *bold_seq_name=NULL,*subset_file=NULL,*subset2_file=NULL,*rm_subset_file=NULL,
    *rm_subset2_file=NULL;
  float bundling_threshold=0.9;
  int exit_code=0,count_sequence_errors=0,please_print_snps=0,
    report_consensus_seqs=0,report_major_allele=0;
  int show_allele_evidence=0,please_collapse_lines=0,keep_all_links=0;
  int remove_listed_seqs=0,remove_listed_seqs2=0,please_report_similarity;
  char *reference_seq_name="CONSENS%d",*clustal_out=NULL;
  char* nbT = NULL;
  int nThreads;

  black_flag_init(argv[0],PROGRAM_VERSION);
  
  if (argc<2) {
    fprintf(stderr,"\nUsage: %s [OPTIONS] MATRIXFILE\n"
"Align a set of sequences or alignments using the scores in MATRIXFILE.\n"
"Example: %s -read_fasta multidom.seq -clustal m.aln blosum80.mat\n\n"
"INPUT:\n"
"  -read_fasta FILE       Read in FASTA sequence file.\n"
"  -read_msa FILE         Read in MSA alignment file.\n"
"  -read_msa2 FILE        Read in second MSA file. \n" 
"  -subset FILE           Filter MSA to include list of seqs in file.\n"
"  -subset2 FILE          Filter second MSA to include list of seqs in file.\n"
"  -remove FILE           Filter MSA to exclude list of seqs in file.\n"
"  -remove2 FILE          Filter second MSA to exclude list of seqs in file.\n"
"  -read_msa_list FILE    Read an MSA from each filename listed in file.\n"
"  -tolower               Force FASTA/MSA sequences to lowercase\n"
"                           (nucleotides in our matrix files)\n"
"  -toupper               Force FASTA/MSA sequences to UPPERCASE\n"
"                           (amino acids in our matrix files)\n"
"\nALIGNMENT:\n"
"  -do_global             Do global alignment.\n"
"  -do_progressive        Perform progressive alignment using a guide tree\n"
"                           built by neighbor joining from a set of\n"
"                           sequence-sequence similarity scores.\n"
"  -read_pairscores FILE  Read tab-delimited file of similarity scores.\n"
"                           (If not provided, scores are constructed\n"
"                           using pairwise sequence alignment.)\n"
"  -fuse_all              Fuse identical letters on align rings.\n"
"  -threads				  Number of threads to use\n"
"\nANALYSIS:\n"
"  -hb                    Perform heaviest bundling to generate consensi.\n"
"  -hbmin VALUE           Include in heaviest bundle sequences with\n"
"                           percent ID (as a fraction) >= value.\n"
"\nOUTPUT:\n"
"  -pir FILE              Write out MSA in PIR format.\n"
"  -clustal FILE          Write out MSA in CLUSTAL format.\n"
"  -po FILE               Write out MSA in PO format.\n"
"  -preserve_seqorder     Write out MSA with sequences in their input order.\n"
"  -printmatrix LETTERS   Print score matrix to stdout.\n"
"  -best                  Restrict MSA output to heaviest bundles (PIR only).\n"
"  -v                     Run in verbose mode (e.g. output gap penalties).\n\n"
"  NOTE:  One of the -read_fasta, -read_msa, or -read_msa_list arguments\n" 
"         must be used, since a sequence or alignment file is required.\n\n"
"For more information, see http://www.bioinformatics.ucla.edu/poa.\n\n"
	    ,argv[0],argv[0]);
    exit(-1);
  }
	fasta_out = "default_output_msa.fasta";
  FOR_ARGS(i,argc) { /* READ ALL THE ARGUMENTS */
    ARGGET("-pir",fasta_out); /* SAVE FASTA-PIR FORMAT ALIGNMENT FILE */
    ARGGET("-corrected_reads_fasta",seq_filename); /* READ FASTA FILE FOR ALIGNMENT */
    ARGGET("-uncorrected_reads_fasta",unco_seq_filename); /* READ FASTA FILE FOR ALIGNMENT */
    ARGGET("-reference_reads_fasta",ref_seq_filename); /* READ FASTA FILE FOR ALIGNMENT */
    ARGGET("-threads", nbT);
  }
  nThreads = atoi(nbT);
  
  do_switch_case = 1;

  /** CHECK FOR CONFLICTING FLAGS **/
  
  if (po_list_filename && (po_filename || po2_filename)) {
    WARN_MSG(USERR,(ERRTXT, "Error: The -read_po_list and -read_po flags cannot be used at the same time.\nExiting."), "$Revision: 1.2.2.9 $");
    exit_code = 1;
    goto free_memory_and_exit;
  }
  
  if (((subset_file || rm_subset_file) && !po_filename) || ((subset2_file || rm_subset2_file) && !po2_filename)) {
    WARN_MSG(USERR,(ERRTXT, "Error: Each -subset/-remove flag must have a corresponding -read_po flag.\nExiting."),"$Revision: 1.2.2.9 $");
    exit_code = 1;
    goto free_memory_and_exit;
  }
  
  if ((subset_file && rm_subset_file) || (subset2_file && rm_subset2_file)) {
    WARN_MSG(USERR,(ERRTXT, "Error: The -subset and -remove flags cannot be used at the same time.\nExiting."),"$Revision: 1.2.2.9 $");
    exit_code = 1;
    goto free_memory_and_exit;
  }
  
  if (rm_subset_file) {
    subset_file = rm_subset_file;
    remove_listed_seqs = 1;
  }
  if (rm_subset2_file) {
    subset2_file = rm_subset2_file;
    remove_listed_seqs2 = 1;
  }
  
  if (hbmin)
    bundling_threshold=atof(hbmin);  

  if (!matrix_filename ||
      read_score_matrix(matrix_filename,&score_matrix)<=0){/* READ MATRIX */
    WARN_MSG(USERR,(ERRTXT,"Error reading matrix file %s.\nExiting",
		    matrix_filename ? matrix_filename: "because none specified"),"$Revision: 1.2.2.9 $");
    exit_code=1; /* SIGNAL ERROR CONDITION */
    goto free_memory_and_exit;
  }
  
  if (print_matrix_letters) /* USER WANTS US TO PRINT A MATRIX */
    print_score_matrix(stdout,&score_matrix,print_matrix_letters
		       /*"ARNDCQEGHILKMFPSTWYV"*/);
  
  
  /** READ INPUT FILES **/
  
  n_input_seqs = 0;
  
  if (seq_filename && unco_seq_filename && ref_seq_filename) {
    seqCorrected_ifile = fopen (seq_filename, "r");
    seqUncorrected_ifile = fopen (unco_seq_filename, "r");
    seqReference_ifile = fopen (ref_seq_filename, "r");
    if (seqCorrected_ifile == NULL || seqUncorrected_ifile == NULL || seqReference_ifile == NULL) {
      WARN_MSG(USERR,(ERRTXT,"Couldn't open sequence file %s.\nExiting",
		      seq_filename),"$Revision: 1.2.2.9 $");
      exit_code=1; /* SIGNAL ERROR CONDITION */
      goto free_memory_and_exit;
    }
    nseq = read_fasta (seqCorrected_ifile, &seq, do_switch_case, &comment);
    nseq = read_fasta (seqUncorrected_ifile, &seqUnco, do_switch_case, &comment);
    nseq = read_fasta (seqReference_ifile, &seqRef, do_switch_case, &comment);
    fclose (seqCorrected_ifile);
    fclose (seqUncorrected_ifile);
    fclose (seqReference_ifile);
    if (nseq == 0) {
      WARN_MSG(USERR,(ERRTXT,"Error reading sequence file %s.\nExiting",
		      seq_filename),"$Revision: 1.2.2.9 $");
      exit_code=1; /* SIGNAL ERROR CONDITION */
      goto free_memory_and_exit;
    }
    seq_ifile=fopen(fasta_out,"w");
    nThreads = nseq > nThreads ? nThreads : nseq;
    pthread_t threads[nThreads];
    
    seqArray = malloc(nThreads * sizeof(int*));
    sizeArray = calloc(nThreads, sizeof(int));
    int cpt = 0;
	seqArray[0] = malloc((nseq / nThreads + nseq % nThreads) * sizeof(int));
	int j;
	for (j = 0; j < (nseq / nThreads + nseq % nThreads); j++) {
		seqArray[0][j] = cpt;
		cpt++;
		sizeArray[0]++;
	}
    for (i = 1; i < nThreads; i++) {
		seqArray[i] = malloc(nseq / nThreads * sizeof(int));
		for (j = 0 ; j < nseq / nThreads; j++) {
			seqArray[i][j] = cpt;
			sizeArray[i]++;
			cpt++;
		}
	}
	
	for (i = 0; i < nThreads; i++) {
		ThreadsParams params = malloc(sizeof(*params));
		params->nb = i;
		params->lpo_out = NULL;
		params->in_seqs = NULL;
		CALLOC (params->in_seqs, max_input_seqs, LPOSequence_T *);
		pthread_create(&threads[i], NULL, doThreadJob, params);
	}
	
	for (i = 0; i < nThreads; i++) {
		pthread_join(threads[i], NULL);
	}
	
	pthread_mutex_destroy(&mutex);
    fclose(seq_ifile);

  }

 free_memory_and_exit: /* FREE ALL DYNAMICALLY ALLOCATED DATA!!!! */
 
  if (nseq>0) FREE (seq);
  exit (exit_code);
}

void freeMem(LPOSequence_T **input_seqs, int nseq, LPOSequence_T *seq, int n_input_seqs){
	int i,j;
	for (i=0;i<n_input_seqs;i++) {
    free_lpo_sequence(input_seqs[i],(j==nseq));
  }
  FREE (input_seqs);
  if (nseq>0) FREE (seq);
}

void buildAndAnalysePOMSA (int n_input_seqs, LPOSequence_T* lpo_out, LPOSequence_T **input_seqs, ResidueScoreMatrix_T score_matrix, FILE* seq_ifile) {
	lpo_out = buildup_progressive_lpo(n_input_seqs, input_seqs, &score_matrix, use_aggressive_fusion, do_progressive, pair_score_file, matrix_scoring_function, do_global, do_preserve_sequence_order);
	  
	if (comment) { /* SAVE THE COMMENT LINE AS TITLE OF OUR LPO */
		FREE(lpo_out->title);
		lpo_out->title=strdup(comment);
	}
	
	pthread_mutex_lock(&mutex);
    write_lpo_bundle_as_fasta(seq_ifile,lpo_out,score_matrix.nsymbol,
		score_matrix.symbol,ibundle);
	pthread_mutex_unlock(&mutex);
}



static LPOSequence_T *read_partial_order_file (char *po_filename, char *subset_filename, int remove_listed_seqs, int keep_all_links, int do_switch_case, ResidueScoreMatrix_T *mat)
{
  LPOSequence_T *lpo_in;
  FILE *po_file=NULL, *subset_file=NULL;
  
  if (!po_filename)
    return NULL;
  
  po_file = fopen (po_filename, "r");
  if (!po_file) {
    WARN_MSG (USERR, (ERRTXT,"Couldn't open MSA file %s.\nExiting.",po_filename), "$Revision: 1.2.2.9 $");
    return NULL;
  }
  
  if (subset_filename) {
    subset_file = fopen (subset_filename, "r");
    if (!subset_file) {
      WARN_MSG (USERR, (ERRTXT,"Couldn't open subset file %s.\nExiting.",subset_filename), "$Revision: 1.2.2.9 $");
      return NULL;
    }
  }
  
  if (subset_file) {
    lpo_in = read_msa_select (po_file, UNKNOWN_MSA, subset_file, keep_all_links, remove_listed_seqs, do_switch_case, mat);
    fclose (subset_file);
    fclose (po_file);
    if (lpo_in==NULL || lpo_in->nsource_seq == 0) {
      WARN_MSG (USERR, (ERRTXT,"MSA file %s, filtered with subset file %s, couldn't be read or contains no sequences.\nExiting.", po_filename, subset_filename), "$Revision: 1.2.2.9 $");
      return NULL;
    }
  }
  else {
    lpo_in = read_msa (po_file, UNKNOWN_MSA, do_switch_case, mat);
    fclose (po_file);
    if (lpo_in==NULL || lpo_in->nsource_seq == 0) {
      WARN_MSG (USERR, (ERRTXT,"MSA file %s couldn't be read or contains no sequences.\nExiting.", po_filename), "$Revision: 1.2.2.9 $");
      return NULL;
    }
  }

  return lpo_in;
}
