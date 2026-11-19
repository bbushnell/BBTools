package test.benchmark;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

import benchmark.CallGenesLengthScorecard;
import benchmark.CallGenesPartialScorecard;
import benchmark.CallGenesBenchmarkSupport;
import benchmark.CallGenesStopAudit;
import benchmark.CallGenesStopAudit.ReferenceSet;
import gff.GffLine;
import shared.Shared;

/**
 * Literal saved-GFF oracles for length bands and partial-call normalization.
 * Uses synthetic references, calls and a frozen-locus table to check last-wins
 * binding, physical contig edges, duplicate occurrences and conserved counters.
 * Each exclusion policy is graded with and without the optional band sidecars.
 * Run in a dedicated Dori fixture process: thread and CDS-only settings are global.
 */
public final class CallGenesPartialBandScorecardTest{

	/**
	 * Create and grade the synthetic files beneath one new directory.
	 * @param args Exactly one path whose directory must not already exist.
	 * @throws Exception If fixture setup, grading or output reads fail.
	 */
	public static void main(String[] args)throws Exception{
		check(args.length==1,"One fresh fixture directory required");
		final Path root=Paths.get(args[0]);
		Files.createDirectory(root);
		Shared.setThreads(1);
		CallGenesBenchmarkSupport.cdsOnly();
		// Keep the literal reference/call rows together so the edge cases are visible.
		final String refs=g("c",400,549,0,0,".")+g("c",700,999,0,0,".")+g("c",1300,1602,0,0,".")+
			g("c",1800,1949,0,0,".")+g("c",2000,2299,0,0,".")+g("c",2150,2299,0,0,"partial=01")+
			g("s",1,150,0,2,"partial=10;start_range=.,1")+g("t",153,302,1,2,"partial=01;end_range=302,.")+
			g("e",1,300,0,0,".")+g("e",1,301,0,0,"partial=01;end_range=301,.")+
			g("p",1,301,0,0,"partial=01;end_range=301,.")+g("m",1,301,1,0,"partial=10;start_range=.,1")+
			g("c",3400,3459,0,0,".")+g("c",3500,4102,0,0,".")+g("c",4200,4652,0,0,".");
		final String calls=g("c",247,549,0,0,".")+g("c",247,549,0,0,".")+g("c",850,999,0,0,".")+
			g("c",1453,1602,0,0,".")+g("c",2000,2299,0,0,".")+g("s",3,150,0,0,".")+
			g("t",153,300,1,0,".")+g("e",1,300,0,0,".")+g("p",1,300,0,0,".")+g("m",2,301,1,0,".")+
			g("c",2600,2700,0,0,"partial=10")+g("c",2800,2952,0,0,".")+
			g("c",3000,3150,0,0,".")+g("c",3200,3352,0,1,".");
		write(root.resolve("references.gff"),refs);
		write(root.resolve("calls.gff"),calls);
		final StringBuilder fasta=new StringBuilder();
		for(String contig:new String[]{"c","s","t","e","p","m"}){
			final int length=contig.equals("c")?5000:contig.equals("s")||contig.equals("t")?302:301;
			fasta.append('>').append(contig).append('\n');
			for(int i=0;i<length;i++){fasta.append('A');}
			fasta.append('\n');
		}
		write(root.resolve("reference.fa"),fasta.toString());
		write(root.resolve("grades.tsv"),"genome\tarm\t"+String.join("\t",CallGenesBenchmarkSupport.GRADE_NAMES)+
			"\ng0001\tbonus400\t1\t8\t13\t7\t1\t8\t13\t6\n");
		write(root.resolve("expected.tsv"),CallGenesPartialScorecard.EXPECTED_HEADER+"\nfixture\t1\t14\t8\t7\t6\t13\t26\n");
		write(root.resolve("inputs.tsv"),CallGenesPartialScorecard.INPUT_HEADER+"\nfixture\tg0001\t"+
			root.resolve("reference.fa")+"\t"+root.resolve("references.gff")+"\t"+root.resolve("calls.gff")+"\t"+root.resolve("grades.tsv")+"\tbonus400\n");
		final ReferenceSet winners=CallGenesStopAudit.loadReference(root.resolve("references.gff").toString());
		check(winners.uniqueStops()==14 && winners.duplicateKeys==1,"Literal last winner hides one full reference");
		final Path loci=root.resolve("loci.tsv");
		writeLoci(loci,winners);
		final CallGenesLengthScorecard.ReferenceBands binding=CallGenesLengthScorecard.readSavedBands(loci,Collections.singleton("g0001")).bind("g0001",winners);
		check(binding.annotationBand(winners.lookup("c",2299))==6,"Original partial winner retains frozen unknown band");
		boolean foreign=false;
		try{
			binding.annotationBand(line(g("c",2000,2299,0,0,".")));
		}catch(IllegalArgumentException expected){
			foreign=true;
		}
		check(foreign,"Accessor cannot resurrect a shadowed reference at the same stop");
		final List<GffLine> query=GffLine.loadGffFile(root.resolve("calls.gff").toString(),null,false);
		final String rawBands=binding.grade(query).rows("fixture","g0001").toString();
		for(String exclusion:new String[]{"any","threeprime"}){
			final Path plain=root.resolve(exclusion+"_plain"),banded=root.resolve(exclusion+"_banded");
			run(root,plain,exclusion,null);
			run(root,banded,exclusion,loci);
			for(String file:new String[]{"scorecard.tsv","by_genome.tsv","partial_counts.tsv","dropped_calls.tsv",
				"normalized_dropped_calls.tsv","start_only_partial_starts.tsv","PASS"}){
				check(Arrays.equals(Files.readAllBytes(plain.resolve(file)),Files.readAllBytes(banded.resolve(file))),"Existing output unchanged: "+exclusion+"/"+file);
			}
			check(!Files.exists(plain.resolve("band_scorecard.tsv")),"Band sidecars require explicit saved loci");
			final List<String> rows=Files.readAllLines(banded.resolve("band_scorecard.tsv"),StandardCharsets.US_ASCII);
			check(rows.size()==19 && rows.get(0).split("\t").length==37,"Two scopes, seven bins plus ALL and40_99, additive band field");
			final HashMap<String,Integer> h=header(rows.get(0));
			final String retained=exclusion.equals("any")?"exclude_partial":"terminal_stop_eligible_by_annotation";
			assertGrade(rows,h,"all_genes","ALL",new long[]{1,8,13,7,1,8,13,6},10,5,4,18);
			assertGrade(rows,h,"all_genes","40_99",new long[]{1,6,5,1,1,6,8,3},6,1,1,7);
			if(exclusion.equals("any")){
				assertGrade(rows,h,retained,"ALL",new long[]{1,5,7,4,1,5,10,6},5,4,4,15);
				assertGrade(rows,h,retained,"40_99",new long[]{1,4,3,1,1,4,6,3},4,1,1,7);
			}else{
				assertGrade(rows,h,retained,"ALL",new long[]{1,7,9,4,1,7,12,6},7,4,4,15);
				assertGrade(rows,h,retained,"40_99",new long[]{1,6,5,1,1,6,8,3},6,1,1,7);
			}
			final String[] all=row(rows,"all_genes","ALL");
			check(all[h.get("normalized_correct")].equals("10"),"Both normalized stop-edge calls become matching unknown-band references");
			final String[] shortRow=row(rows,"all_genes","40_99");
			check(shortRow[h.get("duplicate_matching_occurrences")].equals("1"),"Duplicate query occurrences do not inflate unique-key recall");
			for(String raw:rawBands.trim().split("\n")){
				final String[] r=raw.split("\t"),actual=row(rows,"all_genes",r[2]);
				check(actual[h.get("reference_keys")].equals(r[3]),"Frozen annotation denominator matches accepted raw grader: "+r[2]);
				for(int i=0;i<8;i++){check(actual[h.get(CallGenesBenchmarkSupport.GRADE_NAMES[i])].equals(r[i+11]),"Raw band parity: "+r[2]+"/"+i);}
			}
			conservation(rows,h,"all_genes");
			conservation(rows,h,retained);
			final List<String> perGenome=Files.readAllLines(banded.resolve("band_by_genome.tsv"),StandardCharsets.US_ASCII);
			for(int i=1;i<rows.size();i++){
				check(rows.get(i).replace("fixture\tALL\t","fixture\tg0001\t").equals(perGenome.get(i)),"Single-genome sidecar aggregation");
			}
		}
		final Path wrong=root.resolve("wrong_loci.tsv");
		write(wrong,new String(Files.readAllBytes(loci),StandardCharsets.US_ASCII).replace("g0001","g0002"));
		boolean rejected=false;
		try{
			run(root,root.resolve("wrong_roster"),"any",wrong);
		}catch(IllegalArgumentException expected){
			rejected=true;
		}
		check(rejected && !Files.exists(root.resolve("wrong_roster/PASS")),"A different saved cohort fails before publication");
		System.out.println("PARTIAL_BAND_SCORECARD_TEST_PASS checks="+checks+" caller_runs=0 classifier_runs=0 hashes=0");
	}

	/** Literal eight-field exact oracle plus independent normalized summary counts. */
	private static void assertGrade(List<String> rows,HashMap<String,Integer> h,String scope,String band,long[] raw,
			long correct,long missed,long extras,long startErrors){
		final String[] row=row(rows,scope,band);
		for(int i=0;i<8;i++){
			check(Long.parseLong(row[h.get(CallGenesBenchmarkSupport.GRADE_NAMES[i])])==raw[i],"Literal raw "+scope+"/"+band+"/"+i);
		}
		check(Long.parseLong(row[h.get("normalized_correct")])==correct,"Normalized TP "+scope+"/"+band);
		check(Long.parseLong(row[h.get("normalized_missed")])==missed,"Normalized FN "+scope+"/"+band);
		check(Long.parseLong(row[h.get("normalized_extras")])==extras,"Normalized extras "+scope+"/"+band);
		check(Long.parseLong(row[h.get("normalized_start_errors")])==startErrors,"Normalized start errors "+scope+"/"+band);
	}

	/** Sum extensive fields across bands into ALL and 40_99; skip ratio columns. */
	private static void conservation(List<String> rows,HashMap<String,Integer> h,String scope){
		final String[] all=row(rows,scope,"ALL"),shortRow=row(rows,scope,"40_99");
		for(int col=4;col<all.length;col++){
			if(col==h.get("recall") || col==h.get("precision") || col==h.get("normalized_recall") || col==h.get("normalized_precision")){
				continue;
			}
			long sum=0,shortSum=0;
			for(int b=0;b<BANDS.length;b++){
				final long n=Long.parseLong(row(rows,scope,BANDS[b])[col]);
				sum+=n;
				if(b==1 || b==2){shortSum+=n;}
			}
			check(sum==Long.parseLong(all[col]) && shortSum==Long.parseLong(shortRow[col]),"Band conservation "+scope+" column"+col);
		}
	}

	/** Invoke the saved-GFF grader with an optional frozen-locus binding. */
	private static void run(Path root,Path out,String exclusion,Path loci)throws Exception{
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList("in="+root.resolve("inputs.tsv"),"expected="+root.resolve("expected.tsv"),"out="+out,"exclude="+exclusion));
		if(loci!=null){args.add("loci="+loci);}
		CallGenesPartialScorecard.main(args.toArray(new String[args.size()]));
	}

	/** Write the saved-locus schema and literal annotation bands for winning stops. */
	private static void writeLoci(Path path,ReferenceSet refs)throws Exception{
		final String[] fields=new String[42];
		final String[] names={"genome","contig","true_stop_1based","strand","reference_start_1based","reference_left_1based","reference_right_1based","phase","annotation_band"};
		for(int i=0;i<42;i++){
			fields[i]=i<names.length?names[i]:"unused"+i;
		}
		final StringBuilder out=new StringBuilder(String.join("\t",fields)).append('\n');
		for(GffLine ref:refs.records()){
			Arrays.fill(fields,"0");
			fields[0]="g0001";
			fields[1]=ref.seqid;
			fields[2]=""+ref.trueStop();
			fields[3]=""+ref.strand;
			fields[4]=""+ref.trueStart();
			fields[5]=""+ref.start;
			fields[6]=""+ref.stop;
			fields[7]=""+ref.phase;
			String band="unknown";
			if(ref.seqid.equals("c")){
				final int stop=ref.trueStop();
				band=stop==549||stop==1949?"40_64":stop==999?"65_99":stop==1602?"100_149":stop==3459?"lt40":stop==4102?"ge200":stop==4652?"150_199":"unknown";
			}else if(ref.seqid.equals("s")){
				band="40_64";
			}else if(ref.seqid.equals("t") || ref.seqid.equals("e") && ref.trueStop()==300){
				band="65_99";
			}
			fields[8]=band;
			out.append(String.join("\t",fields)).append('\n');
		}
		write(path,out.toString());
	}

	/** Index the emitted column names for named-field assertions. */
	private static HashMap<String,Integer> header(String text){
		final HashMap<String,Integer> h=new HashMap<String,Integer>();
		final String[] fields=text.split("\t");
		for(int i=0;i<fields.length;i++){h.put(fields[i],i);}
		return h;
	}

	/** Find one scope/band row, failing if the grader omitted it. */
	private static String[] row(List<String> rows,String scope,String band){
		for(int i=1;i<rows.size();i++){
			final String[] r=rows.get(i).split("\t");
			if(r[2].equals(scope) && r[3].equals(band)){return r;}
		}
		throw new AssertionError("Missing "+scope+"/"+band);
	}

	/** Encode one synthetic CDS as a GFF line with literal coordinates. */
	private static String g(String c,int left,int right,int strand,int phase,String attrs){
		return c+"\tfixture\tCDS\t"+left+"\t"+right+"\t.\t"+(strand==0?"+":"-")+"\t"+phase+"\t"+attrs+"\n";
	}

	/** Parse one synthetic ASCII GFF record. */
	private static GffLine line(String text){
		return new GffLine(text.trim().getBytes(StandardCharsets.US_ASCII));
	}

	/** Write the fixture's ASCII-only table and sequence text. */
	private static void write(Path path,String text)throws Exception{
		Files.write(path,text.getBytes(StandardCharsets.US_ASCII));
	}

	/** Count explicit assertions even when JVM assertions are disabled. */
	private static void check(boolean pass,String message){
		checks++;
		if(!pass){throw new AssertionError(message);}
	}

	private static final String[] BANDS={"lt40","40_64","65_99","100_149","150_199","ge200","unknown"};
	private static int checks;
}
