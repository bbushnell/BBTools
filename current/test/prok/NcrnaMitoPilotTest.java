package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import prok.NcrnaMitoPilot;
import prok.NcrnaMitoPilot.MappedHit;
import structures.ByteBuilder;

/** Independent coordinate oracles and whole-CLI selection/conservation fixtures.
 * @author Ganyu */
public final class NcrnaMitoPilotTest {
	public static void main(String[] args){
		require(args.length==1 && args[0].startsWith("outdir="),"Require outdir=NEW_DIRECTORY");
		final String root=args[0].substring(7);require(new File(root).mkdir(),"Fixture directory must be new");
		final byte[] a="ACGTACGTAA".getBytes(StandardCharsets.US_ASCII);
		require(new String(NcrnaMitoPilot.expand(a,5,true),StandardCharsets.US_ASCII).equals("ACGTACGTAAACGT"),"Circular search appends the declared prefix exactly");
		require(NcrnaMitoPilot.expand(a,100,true).length==19 && NcrnaMitoPilot.expand(a,100,false).length==10,"At most length-1 prefix bases, and no extension of linear genomes");
		final MappedHit plus=NcrnaMitoPilot.mapHit(9,13,'+',10,4),minus=NcrnaMitoPilot.mapHit(13,9,'-',10,4),dup=NcrnaMitoPilot.mapHit(12,14,'+',10,4),ordinary=NcrnaMitoPilot.mapHit(2,7,'+',10,0);
		require(plus.from==9 && plus.to==3 && plus.strand=='+' && plus.crossesOrigin && !plus.appendedDuplicate,"Origin crossing plus hit preserves five physical bases");
		require(minus.from==3 && minus.to==9 && minus.strand=='-' && minus.crossesOrigin && !minus.appendedDuplicate,"Origin crossing minus hit keeps transcript-oriented endpoints");
		require(dup.from==2 && dup.to==4 && dup.appendedDuplicate && !dup.crossesOrigin,"A wholly appended hit is a duplicate representation");
		require(ordinary.from==2 && ordinary.to==7 && !ordinary.crossesOrigin && !ordinary.appendedDuplicate,"Linear ordinary coordinates remain unchanged");
		require(NcrnaMitoPilot.mapHit(3,3,'-',10,4).strand=='-',"A single-base reverse hit cannot infer its strand from equal endpoints");
		rejectGeometry(1,12,10,4);rejectGeometry(0,4,10,4);rejectGeometry(11,15,10,4);
		final String header="accession\ttaxid\torganelle_class\torganism\tstatus\tbases\ttopology\trrna_features\textracted\tunsupported\tgbff\tsequence_source\tsequence_path\tcontig_location\n";
		final String lh="locus_id\taccession\ttaxid\torganelle_class\trrna_index\traw_location\tproduct\tstatus\tpartial\tpseudo\tpseudogene\tstrand\tmin_coordinate\tmax_coordinate\tsegments_in_transcript_order\tcrosses_origin\textracted_length\tgbff\n";
		final ByteBuilder records=new ByteBuilder(header),reverse=new ByteBuilder(header),loci=new ByteBuilder(lh),genomes=new ByteBuilder(),candidates=new ByteBuilder();
		final String[] lines=new String[24];
		for(int i=0;i<24;i++){
			final String id="G"+i+".1",org=i==23?"mitochondrion:kinetoplast":"mitochondrion";final int n=i<12?1:0,length=i==10?31:20+i;
			final ByteBuilder row=new ByteBuilder();row.append(id).tab().append(100+i).tab().append(org).append("\tfixture\t").append(n>0?"CANDIDATES":"NO_RRNA")
				.tab().append(length).append("\tcircular\t").append(n).tab().append(n).append("\t0\tfixture.gbff\t").append(i==5?"FNA_CONTIG":"GBFF_ORIGIN")
				.append("\tfixture.fa\tNA\n");lines[i]=row.toString();records.append(lines[i]);
			genomes.append('>').append(id).append(" tid_").append(100+i).nl();for(int j=0;j<length;j++){genomes.append("ACGT".charAt(j&3));}genomes.nl();
			if(n>0){
				final int len=i==1?9:8;
				loci.append(id).append("_rrna_1\t").append(id).tab().append(100+i).tab().append(org).append("\t1\t1..").append(len).append("\t16S ribosomal RNA\tEXTRACTED\t")
					.append(i==7).tab().append(i==11).append("\tNA\t+\t1\t").append(len).append("\t1..").append(len).append(":+\tfalse\t").append(len).append("\tfixture.gbff\n");
				candidates.append('>').append(id).append("_rrna_1 tid_").append(100+i).nl().append(i==1?"ACGTACGTA":"ACGTACGT").nl();
			}
		}
		for(int i=23;i>=0;i--){reverse.append(lines[i]);}
		write(root+"/records.tsv",records.toString());write(root+"/reverse.tsv",reverse.toString());write(root+"/loci.tsv",loci.toString());
		write(root+"/genomes.fa",genomes.toString());write(root+"/candidates.fa",candidates.toString());
		run(root,"normal",root+"/records.tsv",root+"/candidates.fa",24,10000);
		run(root,"reverse",root+"/reverse.tsv",root+"/candidates.fa",24,10000);
		for(String suffix:new String[]{"selection.tsv","records.tsv","loci.tsv.gz","strata.tsv","summary.tsv","population.tsv","genomes.fa.gz","search.fa.gz","candidates.fa.gz","canary.search.fa.gz","canary.candidates.fa.gz"}){
			require(read(root+"/normal."+suffix).equals(read(root+"/reverse."+suffix)),"Source record order must not change the deterministic pilot: "+suffix);
		}
		final String selection=read(root+"/normal.selection.tsv"),strata=read(root+"/normal.strata.tsv"),summary=read(root+"/normal.summary.tsv");
		require(lines(selection)==17 && summary.contains("\n16\t8\t8\t"),"Pilot must retain eight annotated and eight unknown genomes");
		for(String id:new String[]{"G23.1","G7.1","G11.1","G5.1"}){require(selection.contains("\n"+id+"\t"),"Mandatory rare/partial/pseudo/CONTIG member missing: "+id);}
		require(selection.contains("G23.1\t123\tmitochondrion:kinetoplast\t"),"The exact RefSeq organelle qualifier must select the rare stratum and remain unmodified");
		for(String s:new String[]{"kinetoplast","partial","pseudo","CONTIG_fallback"}){require(strata.contains(s+"\t1\t1\tPRESENT\n"),"Mandatory stratum accounting must be exact: "+s);}
		require(headers(read(root+"/normal.genomes.fa.gz"))==16 && headers(read(root+"/normal.search.fa.gz"))==16 && headers(read(root+"/normal.candidates.fa.gz"))==8,"Genome and locus FASTA counts must match selected ledgers");
		require(read(root+"/normal.loci.tsv.gz").contains("\ttrue\tfalse\t") && read(root+"/normal.loci.tsv.gz").contains("\tfalse\ttrue\t"),"Partial/pseudo labels remain in the copied locus ledger");
		reject(root,"budget",root+"/records.tsv",root+"/candidates.fa",24,1,"budget");
		reject(root,"wrong_count",root+"/records.tsv",root+"/candidates.fa",25,10000,"record count");
		write(root+"/duplicate.tsv",records.toString()+lines[0]);reject(root,"duplicate",root+"/duplicate.tsv",root+"/candidates.fa",25,10000,"Duplicate census");
		write(root+"/missing.fa",candidates.toString().replace(">G0.1_rrna_1 tid_100\nACGTACGT\n",""));
		reject(root,"missing",root+"/records.tsv",root+"/missing.fa",24,10000,"Missing extracted source locus");
		write(root+"/short.fa",candidates.toString().replace(">G0.1_rrna_1 tid_100\nACGTACGT\n",">G0.1_rrna_1 tid_100\nACGT\n"));
		reject(root,"short",root+"/records.tsv",root+"/short.fa",24,10000,"Candidate FASTA membership/length");
		write(root+"/lost_kinetoplast.tsv",records.toString().replace("mitochondrion:kinetoplast","mitochondrion"));
		reject(root,"lost_kinetoplast",root+"/lost_kinetoplast.tsv",root+"/candidates.fa",24,10000,"Kinetoplast source count differs");
		run(root,"longtail",root+"/records.tsv",root+"/candidates.fa",24,10000,"selection=longtail");
		run(root,"longtail_reverse",root+"/reverse.tsv",root+"/candidates.fa",24,10000,"selection=longtail");
		final String longest=read(root+"/longtail.selection.tsv");
		require(lines(longest)==7 && longest.equals(read(root+"/longtail_reverse.selection.tsv")),"Longtail must be deterministic, bounded and independent of source ordering");
		for(String id:new String[]{"G10.1","G11.1","G23.1","G7.1","G5.1","G1.1"}){require(longest.contains("\n"+id+"\t"),"Longest group/stratum/candidate or accession tie-break missing: "+id);}
		require(read(root+"/longtail.summary.tsv").contains("\n6\t5\t1\t178\t202\t5\t41\t6\t202\t5\t41\t5\n"),"Longtail complete genome/candidate accounting has a hand-derived oracle");
		final String population=read(root+"/longtail.population.tsv");
		require(population.contains("\ngenomic\t24\t757\t853\t0.001706\t43\n") && population.contains("\ncandidate\t12\t97\t97\t0.000194\t9\n"),"Search denominators must include the full source census and both strands, not selected targets");
		require(population.equals(read(root+"/normal.population.tsv")),"Selection policy cannot change population denominators");
		boolean rejected=false;try{run(root,"longtail_budget",root+"/records.tsv",root+"/candidates.fa",24,177,"selection=longtail");}catch(IllegalArgumentException e){require(e.getMessage().contains("budget"),"Unexpected longtail budget rejection: "+e);rejected=true;}
		require(rejected && !new File(root+"/longtail_budget.summary.tsv").exists(),"Longtail cannot silently replace an oversized required longest member");
		System.out.println("MITO_PILOT_TEST_PASS selected=16 annotated=8 no_rrna=8 mandatory=4 permutation_exact=true circular_both_strands=true expectedkinetoplast=1 geometry_negatives=3 input_negatives=6 longtail=6 longest_candidate=true longest_tie=true population_both_strands=true longtail_budget_negative=true");
	}
	private static void run(String root,String name,String records,String candidates,int expected,int budget,String... extra){
		final java.util.ArrayList<String> args=new java.util.ArrayList<String>(java.util.Arrays.asList("records="+records,"loci="+root+"/loci.tsv","genomes="+root+"/genomes.fa","candidates="+candidates,
			"outprefix="+root+"/"+name,"expectedrecords="+expected,"expectedloci=12","expectedkinetoplast=1","pergroup=8","maxbases="+budget,"window=5","t=1","zl=9"));
		for(String arg:extra){args.add(arg);}NcrnaMitoPilot.main(args.toArray(new String[0]));
	}
	private static void reject(String root,String name,String records,String candidates,int expected,int budget,String message){
		boolean failed=false;try{run(root,name,records,candidates,expected,budget);}catch(IllegalArgumentException e){require(e.getMessage().contains(message),"Unexpected failure: "+e);failed=true;}
		require(failed && !new File(root+"/"+name+".summary.tsv").exists(),"Invalid pilot cannot publish a success summary: "+name);
	}
	private static void rejectGeometry(int from,int to,int len,int pad){boolean failed=false;try{NcrnaMitoPilot.mapHit(from,to,'+',len,pad);}catch(IllegalArgumentException e){failed=true;}require(failed,"Impossible geometry must fail");}
	private static int lines(String s){int n=0;for(int i=0;i<s.length();i++){if(s.charAt(i)=='\n'){n++;}}return n;}
	private static int headers(String s){int n=0;for(int i=0;i<s.length();i++){if(s.charAt(i)=='>' && (i==0 || s.charAt(i-1)=='\n')){n++;}}return n;}
	private static String read(String path){final ByteFile f=ByteFile.makeByteFile(path,true);final ByteBuilder b=new ByteBuilder();try{for(byte[] line;(line=f.nextLine())!=null;){b.append(line).nl();}}finally{require(!f.close(),"Fixture read failed");}return new String(b.array,0,b.length(),StandardCharsets.UTF_8);}
	private static void write(String path,String value){final ByteStreamWriter w=new ByteStreamWriter(FileFormat.testOutput(path,FileFormat.TEXT,null,true,false,false,false));w.start();w.print(value.getBytes(StandardCharsets.UTF_8));require(!w.poisonAndWait(),"Fixture write failed");}
	private static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
}
