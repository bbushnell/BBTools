package test.prok;

import java.io.File;
import map.ObjectMap;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import prok.NcrnaMitoPartitions;
import prok.RrnaResourceIO;
import structures.ByteBuilder;

/** Hand-derived bin loads and exact padded bases test whole-accession partitioning.
 * @author Ganyu, Raiden */
public final class NcrnaMitoPartitionsTest {
	public static void main(String[] args){
		check(args.length==1 && args[0].startsWith("outdir="),"Fresh fixture directory required");final String d=args[0].substring(7);check(new File(d).mkdir(),"Fixture path exists");
		final String header="accession\ttaxid\torganelle_class\torganism\tstatus\tbases\ttopology\trrna_features\textracted\tunsupported\tgbff\tsequence_source\tsequence_path\tcontig_location\n";
		final String a=record("A",8,false,2),b=record("B",6,true,1),c=record("C",4,false,0),z=record("D",10,true,1);
		write(d+"/records.tsv",header+a+b+c+z);write(d+"/reversed.tsv",header+z+c+b+a);
		write(d+"/loci.tsv","locus_id\taccession\ttaxid\torganelle_class\trrna_index\traw_location\tproduct\tstatus\tpartial\tpseudo\tpseudogene\tstrand\tmin_coordinate\tmax_coordinate\tsegments_in_transcript_order\tcrosses_origin\textracted_length\tgbff\n"
			+locus("A_1","A",4)+locus("A_2","A",3)+locus("B_1","B",2)+locus("D_1","D",5));
		final String genomes=">A tid_1\nACGTACGT\n>B tid_1\nAACCGG\n>C tid_1\nTTAA\n>D tid_1\nACGTACGTAA\n";
		write(d+"/genomes.fa",genomes);write(d+"/candidates.fa",">A_1\nACGT\n>A_2\nTTT\n>B_1\nCC\n>D_1\nGGGGG\n");
		run(d,"records.tsv","genomes.fa","ok",2,36);run(d,"reversed.tsv","genomes.fa","permuted",2,36);
		check(read(d+"/ok.membership.tsv").equals(read(d+"/permuted.membership.tsv")),"Metadata order must not affect partition membership");
		check(read(d+"/ok.membership.tsv").contains("genomic\t0\tD\tD\t10\t14\n") && read(d+"/ok.membership.tsv").contains("genomic\t0\tC\tC\t4\t4\n"),"Longest-first bins must be D+C and B+A, eighteen bases each");
		check(read(d+"/ok.membership.tsv").contains("candidate\t0\tA\tA_1\t4\t4\n") && read(d+"/ok.membership.tsv").contains("candidate\t0\tA\tA_2\t3\t3\n"),"All annotation sequences from one accession stay together");
		final ObjectMap<String,String> actual=new ObjectMap<String,String>(String.class,String.class);
		for(int part=0;part<2;part++){RrnaResourceIO.stream(d+"/ok.genomic."+part+".fa",r->{check(actual.put(RrnaResourceIO.first(r.id),new String(r.bases,java.nio.charset.StandardCharsets.US_ASCII))==null,"Duplicate genome in partition union");});}
		check(actual.size()==4 && "ACGTACGT".equals(actual.get("A")) && "AACCGGAACC".equals(actual.get("B")) && "TTAA".equals(actual.get("C")) && "ACGTACGTAAACGT".equals(actual.get("D")),"Linear genomes remain exact and each circular prefix is appended once");
		final ObjectMap<String,String> candidates=new ObjectMap<String,String>(String.class,String.class);
		for(int part=0;part<2;part++){RrnaResourceIO.stream(d+"/ok.candidate."+part+".fa",r->{check(candidates.put(RrnaResourceIO.first(r.id),new String(r.bases,java.nio.charset.StandardCharsets.US_ASCII))==null,"Duplicate candidate in partition union");});}
		check(candidates.size()==4 && "ACGT".equals(candidates.get("A_1")) && "TTT".equals(candidates.get("A_2")) && "CC".equals(candidates.get("B_1")) && "GGGGG".equals(candidates.get("D_1")),"Every complete candidate must remain byte-exact");
		final String partition=read(d+"/ok.partitions.tsv");check(partition.contains("genomic\t0\t2\t18\t") && partition.contains("genomic\t1\t2\t18\t") && partition.contains("candidate\t0\t2\t7\t") && partition.contains("candidate\t1\t2\t7\t"),"Native counts/base loads must match independent hand sums");
		write(d+"/missing.fa",genomes.replace(">C tid_1\nTTAA\n",""));bad(()->run(d,"records.tsv","missing.fa","missing",2,36));
		write(d+"/duplicate.fa",genomes+">A\nACGTACGT\n");bad(()->run(d,"records.tsv","duplicate.fa","duplicate",2,36));
		write(d+"/short.fa",genomes.replace("AACCGG","AACCG"));bad(()->run(d,"records.tsv","short.fa","short",2,36));
		bad(()->run(d,"records.tsv","genomes.fa","total",2,35));bad(()->run(d,"records.tsv","genomes.fa","empty",4,36));
		subset(d,genomes);
		System.out.println("MITO_PARTITIONS_TEST_PASS whole_accessions=true padded_bases=true candidate_bytes=true balanced_bins=true permutation=true negatives=5");
	}
	/** Four genomes have independently known usable counts2/1/0/1 and padded lengths8/10/4/14. */
	static void subset(String d,String genomes){
		final String original=read(d+"/loci.tsv");
		final String labels=original.replace(locus("A_1","A",4),locus("A_1","A",4).replace("\tRNA\t","\t5S\t"))
			.replace(locus("A_2","A",3),locus("A_2","A",3).replace("\tRNA\t","\tLSU\t"));
		write(d+"/labels.tsv",labels);
		subsetRun(d,"labels.tsv","genomes.fa","count",1,true,32);
		check(read(d+"/count.population.tsv").contains("genomic\t3\t20\t28\t"),"Metadata-only selection must count B+C+D: three genomes, twenty original and twenty-eight expanded bases");
		check(read(d+"/count.eligibility.tsv").contains("A\t2\t2\tfalse\t8\t0\t8\n") && read(d+"/count.eligibility.tsv").contains("C\t0\t0\ttrue\t4\t0\t4\n"),"Label-free rule excludes two usable5S+LSU annotations and retains the no-rRNA genome");
		check(!new File(d+"/count.genomic.0.fa").exists() && !new File(d+"/count.partitions.tsv").exists(),"Metadata census must not emit sequence payloads or claim emitted partitions");
		subsetRun(d,"labels.tsv","genomes.fa","selected",1,false,32);
		for(String suffix:new String[]{"selection.tsv","membership.tsv","population.tsv","eligibility.tsv"}){check(read(d+"/count."+suffix).equals(read(d+"/selected."+suffix)),"Census and emission must share exact selection and membership: "+suffix);}
		final ObjectMap<String,String> selected=new ObjectMap<String,String>(String.class,String.class);
		for(int i=0;i<3;i++){RrnaResourceIO.stream(d+"/selected.genomic."+i+".fa",r->{check(selected.put(RrnaResourceIO.first(r.id),new String(r.bases,java.nio.charset.StandardCharsets.US_ASCII))==null,"Subset must not duplicate an accession");});}
		check(selected.size()==3 && selected.get("A")==null && "AACCGGAACC".equals(selected.get("B")) && "TTAA".equals(selected.get("C")) && "ACGTACGTAAACGT".equals(selected.get("D")),"Only eligible complete genomes with exact circular padding may be emitted");
		check(!new File(d+"/selected.genomic.3.fa").exists() && !new File(d+"/selected.candidate.0.fa").exists(),"Cap partitions at available genomes and never rewrite accepted candidate inputs");
		for(String flag:new String[]{"partial","pseudo"}){
			final String row=locus("A_2","A",3),changed=row.replace("\tEXTRACTED\tfalse\tfalse\t",flag.equals("partial")?"\tEXTRACTED\ttrue\tfalse\t":"\tEXTRACTED\tfalse\ttrue\t");
			write(d+'/'+flag+".tsv",original.replace(row,changed));subsetRun(d,flag+".tsv","genomes.fa",flag,1,true,2);
			check(read(d+'/'+flag+".population.tsv").contains("genomic\t4\t28\t36\t") && read(d+'/'+flag+".eligibility.tsv").contains("A\t2\t1\ttrue\t"),"A partial or pseudo annotation must not count as a second usable rRNA");
		}
		subsetRun(d,"labels.tsv","genomes.fa","zero",0,true,32);
		check(read(d+"/zero.population.tsv").contains("genomic\t1\t4\t4\t"),"Threshold zero selects only the independently known unannotated genome");
		write(d+"/missing_excluded.fa",genomes.replace(">A tid_1\nACGTACGT\n",""));bad(()->subsetRun(d,"labels.tsv","missing_excluded.fa","missing_excluded",1,false,2));
		System.out.println("MITO_GENOMIC_SUBSET_TEST_PASS metadata_only=true usable_flags=true label_free=true excluded_audited=true exact_padding=true capped_partitions=true missing_excluded_rejected=true");
	}
	static void subsetRun(String d,String loci,String genomes,String out,int max,boolean metadata,int parts){
		final java.util.ArrayList<String> args=new java.util.ArrayList<String>(java.util.Arrays.asList("records="+d+"/records.tsv","loci="+d+'/'+loci,"outprefix="+d+'/'+out,"genomiconly=t","metadataonly="+metadata,"genomicmaxusable="+max,"partitions="+parts,"window=5","expectedrecords=4","expectedloci=4","expectedsearchbases=36","expectedcandidatebases=14","t=1","fastawrap=1023"));
		if(!metadata){args.add("genomes="+d+'/'+genomes);}NcrnaMitoPartitions.main(args.toArray(new String[args.size()]));
	}
	static void run(String d,String records,String genomes,String out,int parts,int bases){NcrnaMitoPartitions.main(new String[]{"records="+d+'/'+records,"loci="+d+"/loci.tsv","genomes="+d+'/'+genomes,"candidates="+d+"/candidates.fa","outprefix="+d+'/'+out,"partitions="+parts,"window=5","expectedrecords=4","expectedloci=4","expectedsearchbases="+bases,"expectedcandidatebases=14","t=1","fastawrap=1023"});}
	static String record(String id,int length,boolean circular,int loci){return id+"\t1\tmitochondrion\tfixture\tFIXTURE\t"+length+'\t'+(circular?"circular":"linear_or_unspecified")+'\t'+loci+'\t'+loci+"\t0\tfixture.gbff\tGBFF_ORIGIN\tfixture.gbff\tNA\n";}
	static String locus(String id,String accession,int length){return id+'\t'+accession+"\t1\tmitochondrion\t1\t1.."+length+"\tRNA\tEXTRACTED\tfalse\tfalse\tNA\t+\t1\t"+length+"\t1.."+length+":+\tfalse\t"+length+"\tfixture.gbff\n";}
	static void bad(Runnable action){boolean caught=false;try{action.run();}catch(IllegalArgumentException e){caught=true;}check(caught,"Invalid corpus/partition must fail loudly");}
	static void write(String path,String text){final ByteStreamWriter out=RrnaResourceIO.writer(path);try{out.print(text);}finally{RrnaResourceIO.finish(out);}}
	static String read(String path){final ByteFile in=RrnaResourceIO.open(path);final ByteBuilder b=new ByteBuilder();try{for(byte[] line;(line=in.nextLine())!=null;){b.append(line).nl();}}finally{check(!in.close(),"Fixture read failed");}return b.toString();}
	static void check(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
}
