package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.Locale;

import fileIO.FileFormat;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaJointPolicy;
import prok.NcrnaSweepFrontier;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Finite saved-metadata fixtures with literal expected frontiers and input
 * coverage. Missing FASTA/GFF paths demonstrate that payloads are never opened.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaSweepFrontierTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaSweepFrontierTest.class, false);
		try{
			final Parser parser=new Parser();
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown test option: "+arg);
			}
			check(parser.out1!=null, "Supply out=fresh_fixture_directory");
			final NcrnaSweepFrontierTest test=new NcrnaSweepFrontierTest(new File(parser.out1).getAbsoluteFile());
			test.run();
			pp.outstream.println("NCRNA_SWEEP_FRONTIER_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private NcrnaSweepFrontierTest(final File root_){
		root=root_;
		check(root.mkdir(), "Frontier fixtures require a fresh directory");
	}

	private void run(){
		write("taxonomy.tsv", "input_id\tfasta\tdomain\tphylum\ttaxonomy_status\n"
			+taxon("1", "bacteria", "A")+taxon("2", "bacteria", "B")
			+taxon("3", "bacteria", "C")+taxon("4", "bacteria", "C")
			+taxon("5", "bacteria", "D")+taxon("6", "bacteria", "D")
			+taxon("7", "UNKNOWN", "UNKNOWN")+taxon("8", "bacteria", "UNCLASSIFIED"));
		search("search1", "first_runtime");
		search("search2", "second_runtime");
		round("search1", 0, false);
		round("search1", 1, true);
		round("search2", 0, false);
		check(new File(path("search1/round2")).mkdir(), "Save an incomplete round to test explicit gap reporting");
		write("searches.tsv", "group\tfamily\tsearch\troot\nA\tTPP\tfirst\t"+path("search1")
			+"\nB\tTPP\tsecond\t"+path("search2")+"\n");
		report("report");
		final String all=read("report.all.tsv"), front=read("report.frontier.tsv"), strata=read("report.strata.tsv");
		check(rows(all)==15, "Three saved rounds of five policies must retain all fifteen measured rows");
		check(front.contains("A\tTPP\tfirst\t0\tbroad\tfalse\t3\t2\t0\t1\t3\t0\t1\t80\t800"),
			"The breadth-only point must remain available even when its total-TP/work pair is dominated");
		check(row(front, "A\tTPP\tfirst\t0\tbroad\t").contains("\tfalse\ttrue\tUNIQUE_MEASUREMENT\t"),
			"Two TP at eight hundred bases is outside the TP frontier but preserves more shipped inputs");
		check(row(front, "A\tTPP\tfirst\t0\trecall\t").contains("\ttrue\ttrue\tCONFLICTING_TUPLE\t"),
			"Conflicting same-tuple outcomes must be marked without rerunning or merging them");
		check(row(all, "A\tTPP\tfirst\t1\tbroad\t").contains("\tCONSISTENT_DUPLICATE\tfirst/round0/broad\t"),
			"Only agreeing complete per-input vectors can deduplicate repeated tuples");
		check(front.contains("A\tTPP\tfirst\t0\tshipped\ttrue"), "A dominated shipped baseline still belongs in the report");
		check(row(front, "A\tTPP\tfirst\t0\tcheap\t").contains("\t99\t"),
			"Apparent FP must remain visible rather than silently filtering a frontier point");
		check(strata.contains("A\tTPP\tfirst\t0\trecall\tCORPUS_LOWER_QUARTILE\t1\t0\n")
			&& strata.contains("A\tTPP\tfirst\t0\trecall\tUNKNOWN\t0\t1\n")
			&& strata.contains("A\tTPP\tfirst\t0\trecall\tUNCLASSIFIED\t1\t0\n"),
			"First gains and last losses must use per-input TP with separate rare, unknown and unclassified strata");
		check(strata.contains("A\tTPP\tfirst\t0\twaste\tCORPUS_LOWER_QUARTILE\t1\t0\n"),
			"Dominance on total axes must not hide a unique policy's rare-lineage coverage evidence");
		check(read("report.taxonomy.tsv").contains("bacteria\tB\t1\tCORPUS_LOWER_QUARTILE")
			&& read("report.taxonomy.tsv").contains("bacteria\tC\t2\tOTHER_CLASSIFIED"),
			"Nearest-rank quartile includes all tied lineage labels and counts unselected corpus inputs");
		check(read("report.gaps.tsv").contains("A\tfirst\t2\tINCOMPLETE_ROUND_NOT_REPLAYED"),
			"Incomplete saved rounds are explicit and never executed");
		check(all.contains("\tUNMEASURED\n") && all.contains("-20.000000"),
			"Java computes work percentages while whole-caller time stays unmeasured");
		write("searches.tsv", "group\tfamily\tsearch\troot\nA\tTPP\tfirst\t"+path("search1")
			+"\nA\tTPP\tsecond\t"+path("search2")+"\n");
		rejected(() -> report("incompatible"), "Different runtime declarations cannot share a comparable frontier");
		write("searches.tsv", "group\tfamily\tsearch\troot\nB\tTPP\tsecond\t"+path("search2")+"\n");
		for(final String suffix:new String[]{"inputs", "summary"}){
			final String name="search2/round0/sweep."+suffix+".tsv";
			write(name, read(name).replace("aligned_model_bases", "unrecorded_model_work"));
		}
		report("unmeasured");
		check(rows(read("unmeasured.all.tsv"))==5 && rows(read("unmeasured.frontier.tsv"))==1
			&& read("unmeasured.all.tsv").contains("\tNA\tNA\tUNIQUE_MEASUREMENT\t"),
			"Unmeasured model work stays NA and cannot establish a work frontier; shipped remains visible");
		write("search2/round0/sweep.inputs.tsv", read("search2/round0/sweep.inputs.tsv")
			.replace("unrecorded_model_work", "aligned_model_bases").replace("\t100\t1000\t100\n", "\t100\tNA\t100\n"));
		write("search2/round0/sweep.summary.tsv", read("search2/round0/sweep.summary.tsv")
			.replace("unrecorded_model_work", "aligned_model_bases").replace("\t100\t1000\t300\n", "\t100\tNA\t300\n"));
		report("literal_na");
		check(row(read("literal_na.all.tsv"), "B\tTPP\tsecond\t0\tshipped\t").contains("\t100\tNA\t300\t"),
			"A present literal NA model-base counter remains unavailable rather than zero or a proxy estimate");
		write("searches.tsv", "group\tfamily\tsearch\troot\nA\tTPP\tfirst\t"+path("search1")+"\n");
		write("search1/round0/sweep.inputs.tsv", read("search1/round0/sweep.inputs.tsv")+"unexpected\trow\n");
		rejected(() -> report("extra"), "An extra per-input row must not disappear from reported totals");
	}

	private String taxon(final String id, final String domain, final String phylum){
		return id+"\t"+path("absent_"+id+".fa")+"\t"+domain+"\t"+phylum+"\tEXACT_SOURCE_JOIN\n";
	}

	private void search(final String name, final String identity){
		check(new File(path(name)).mkdir(), "Each declared search needs one fixture directory");
		write(name+"/identity.tsv", "runtime\t"+identity+"\n");
		write(name+"/caller.config", "in="+path("absent_1.fa")+"\nin="+path("absent_7.fa")+"\nin="+path("absent_8.fa")+"\n");
		write(name+"/familymap.tsv", "literal family projection\n");
		write(name+"/contract.txt", "schema=2\npatience=latched_eligible_raw_identity\nbackend=native\ntarget=TPP\ngrading=target\n"
			+"cachewindows=128\ncachelength=4096\nreusecandidates=true\noverlap=0.5\nsimd=true\nlauncher=/absent/sweep.sh\n");
		final StringBuilder refs=new StringBuilder("input_id\tfasta\tgff\n");
		// Caller order differs; the off-panel reference must not acquire rows or open payloads.
		for(final String id:new String[]{"8", "2", "7", "1"}){
			refs.append(id).append('\t').append(path("absent_"+id+".fa")).append('\t').append(path("absent_"+id+".gff")).append('\n');
		}
		write(name+"/reference.tsv", refs.toString());
	}

	/** Independent literal per-input patterns and work totals; no SUT-derived oracle. */
	private void round(final String search, final int round, final boolean conflict){
		final String directory=search+"/round"+round;
		check(new File(path(directory)).mkdir(), "Round fixtures are fresh");
		write(directory+"/complete.txt", "native\tEVALUATION_RETURNED_AND_SUMMARY_PARSED\n");
		final String[] names={"shipped", "cheap", "recall", "broad", "waste"};
		final int[][] matches={{0, 1, 0}, {0, 0, 0}, {1, 0, 1}, {1, 1, 0}, {1, 0, 1}};
		if(conflict){matches[2]=new int[]{1, 0, 0};}
		final int[] work={100, 30, 60, 80, 90};
		final StringBuilder policies=new StringBuilder(NcrnaJointPolicy.HEADER).append('\n');
		final StringBuilder inputs=new StringBuilder("input_id\tpolicy\t"+METRICS+"\n");
		final StringBuilder summary=new StringBuilder("policy\t"+METRICS+"\n");
		for(int p=0; p<names.length; p++){
			policies.append(names[p]).append("\ttpp\tadaptive\t").append(p+1).append("\t1\t1\t0\t0\t")
				.append(p==0 ? -1 : 1).append('\n');
			final int tp=p==0 || p==2 && conflict ? 1 : p==1 ? 0 : 2;
			final int lost=p==0 || p==3 ? 0 : 1;
			final int gained=tp-1+lost;
			summary.append(names[p]).append('\t').append(tp).append('\t').append(p==1 ? 99 : 0)
				.append('\t').append(3-tp).append("\t3\t").append(lost).append('\t').append(gained)
				.append('\t').append(work[p]).append('\t').append(work[p]*10).append("\t300\n");
		}
		final String[] ids={"1", "7", "8"};
		for(int i=0; i<ids.length; i++){
			for(int p=0; p<names.length; p++){
				final int tp=matches[p][i], baseline=matches[0][i];
				final int cost=i==0 ? work[p] : 0;
				inputs.append(ids[i]).append('\t').append(names[p]).append('\t').append(tp)
					.append('\t').append(p==1 && i==0 ? 99 : 0).append('\t').append(1-tp).append("\t1\t")
					.append(baseline>tp ? 1 : 0).append('\t').append(tp>baseline ? 1 : 0)
					.append('\t').append(cost).append('\t').append(cost*10).append("\t100\n");
			}
		}
		write(directory+"/policies.tsv", policies.toString());
		write(directory+"/sweep.inputs.tsv", inputs.toString());
		write(directory+"/sweep.summary.tsv", summary.toString());
	}

	private void report(final String out){
		NcrnaSweepFrontier.main(new String[]{"in="+path("searches.tsv"), "taxonomy="+path("taxonomy.tsv"), "out="+path(out)});
	}
	private String path(final String name){return new File(root, name).getPath();}
	private void write(final String name, final String text){
		try(Output out=new Output(FileFormat.testOutput(path(name), FileFormat.TXT, null, true, true, false, false))){
			out.writer.print(new ByteBuilder(text.getBytes(StandardCharsets.UTF_8)));
		}
	}
	private String read(final String name){
		final ByteBuilder text=new ByteBuilder();
		try(Input in=new Input(path(name))){
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){text.append(line).nl();}
		}
		return new String(text.toBytes(), StandardCharsets.UTF_8);
	}
	private static int rows(final String text){
		int rows=-1;
		int start=0;
		for(int i=0; i<text.length(); i++){
			if(text.charAt(i)=='\n'){
				if(i>start && text.charAt(start)!='#'){rows++;}
				start=i+1;
			}
		}
		return rows;
	}
	private static String row(final String text, final String prefix){
		final int start=text.indexOf(prefix);
		check(start>=0, "Missing expected report row: "+prefix);
		return text.substring(start, text.indexOf('\n', start));
	}
	private static void rejected(final Runnable action, final String why){
		try{action.run();}
		catch(IllegalArgumentException expected){return;}
		throw new AssertionError(why);
	}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
	private final File root;
	private static final String METRICS="TP\tFP\tFN\tscored_references\tlost_baseline_TP\tgained_TP\tkernel_logical_alignments\taligned_model_bases\tCDS";
}
