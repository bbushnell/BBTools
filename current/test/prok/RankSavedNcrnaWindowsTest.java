package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.Locale;

import cardinality.DynamicDemiLog;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaDdlSavedData;
import prok.NcrnaDdlTrainingOverlap;
import prok.NcrnaDdlLocusOverlap;
import prok.NcrnaMemberLengthBounds;
import prok.NcrnaSavedRankSummary;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.RankSavedNcrnaWindows;
import shared.Shared;
import structures.ByteBuilder;

/** Literal two-model saved caches, with independent ranks and attribution.
 * No production ranking code constructs the expected control tuples.
 * @author Brian, Raiden
 */
public final class RankSavedNcrnaWindowsTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, RankSavedNcrnaWindowsTest.class, false);
		final int previous=DynamicDemiLog.exponentBits();
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown saved-cache fixture parameter: "+arg);
			}
			check(parser.out1!=null && parser.in1==null && parser.in2==null && parser.out2==null
				&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"Use out=fresh_fixture_directory and normal runtime options");
			boolean assertions=false;
			assert(assertions=true);
			check(assertions, "Literal saved-cache fixtures require assertions");
			final File root=new File(parser.out1).getAbsoluteFile();
			check(root.mkdir(), "Saved-cache fixtures must not overwrite existing evidence");
			DynamicDemiLog.setExponent(4);
			run(root.getPath());
			pp.outstream.println("NCRNA_SAVED_WINDOW_GRID_TEST_PASS");
		}finally{
			DynamicDemiLog.setExponent(previous);
			Shared.closeStream(pp.outstream);
		}
	}

	private static void run(final String root){
		final Fixture good=new Fixture(root, "good");
		final RankSavedNcrnaWindows.Config c=good.config("one", 1);
		RankSavedNcrnaWindows.run(c);
		RankSavedNcrnaWindows.run(good.config("two", 2));
		for(final String suffix:new String[]{".methods.tsv", ".summary.tsv", ".loci.tsv.gz", ".orientation.tsv", ".last_losses.tsv"}){
			equal(good.dir+"/one"+suffix, good.dir+"/two"+suffix);
		}
		try(Input in=new Input(good.dir+"/one.summary.tsv")){
			check(text(in.reader.nextLine()).startsWith("seen_stratum\tmethod\t"), "The summary must retain its explicit seen-status column");
			for(int stratum=0; stratum<3; stratum++){
				for(int method=0; method<=96; method++){
					final String expected=(stratum==0 ? "NOT_SEEN" : stratum==1 ? "SEEN" : "UNKNOWN")+'\t'+method
						+(stratum==2 ? "\t0\t0\tNA\t0\t0\t0\t-1\t-1\t-1" : "\t1\t1\t1\t1\t1\t1\t1\t1\t1");
					check(text(in.reader.nextLine()).equals(expected+"\tSAVED_LOCUS_MINIMUM_NOT_CALLER_ALIGNMENTS"),
						"Literal A-member loci rank first; ambiguity/gaps/zero-native do not enter the common denominator");
				}
			}
			check(in.reader.nextLine()==null, "No pooled or duplicate rows may extend the three complete strata");
		}
		try(Input in=new Input(good.dir+"/one.last_losses.tsv")){
			check(in.reader.nextLine()!=null && in.reader.nextLine()==null, "Both common inputs retain their only saved reachable locus");
		}
		loci(good.dir+"/one.loci.tsv.gz");
		try(Input in=new Input(good.dir+"/two.work.tsv")){
			boolean found=false;
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				if(text(row).equals("workers_with_scored_windows\t2")){found=true;}
			}
			check(found, "Both workers must score an initial window of the shared-parent fixture, not merely start");
		}
		final String[] faults={"inventory", "slot", "missing_sequence", "duplicate_sequence", "event", "link",
			"unique", "missing_method", "duplicate_method", "tuple", "strand"};
		final String[] reasons={"counts must match", "input IDs", "lose its sequence", "one exact-length", "event metadata",
			"rank/tie/model parity", "cached unique link", "all original method rows", "only once", "rank/tie/model parity", "strand"};
		check(faults.length==reasons.length, "Every adversarial cache has one intended rejection");
		for(int i=0; i<faults.length; i++){
			final Fixture bad=new Fixture(root, faults[i]);
			try{
				if(faults[i].equals("tuple") || faults[i].equals("link")){RankSavedNcrnaWindows.run(bad.config("rejected", 2));}
				else{NcrnaDdlSavedData.load(bad.dir+"/panels.tsv", bad.dir+"/selection.tsv", 2, 3);}
			}catch(IllegalArgumentException | IllegalStateException failure){
				Throwable cause=failure;
				while(cause.getCause()!=null){cause=cause.getCause();}
				check(cause.getMessage().contains(reasons[i]), "Wrong rejection for "+faults[i]+": "+cause);
				check(!new File(bad.dir+"/rejected.summary.tsv").exists(), "An invalid cache cannot publish a success-looking summary");
				continue;
			}
			throw new AssertionError("Invalid saved cache was accepted: "+faults[i]);
		}
	}

	private static void loci(final String path){
		int rows=0;
		final parse.LineParser1 p=new parse.LineParser1('\t');
		try(Input in=new Input(path)){
			check(in.reader.nextLine()!=null, "Per-locus evidence requires a header");
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				check(p.terms()==9, "Each locus row names its original source disposition");
				final int method=p.parseInt(5);
				final int rank=p.parseInt(6);
				final boolean beta=p.termEquals("beta", 0);
				check(beta || Arrays.equals(p.parseByteArray(0), PANEL_BYTES), "Panel names must round-trip their exact UTF-8 bytes");
				final boolean p0=p.termEquals("p0", 3), zero=p.termEquals("p4", 3);
				check(rank==(p0 ? 1 : zero && method>0 ? 2 : -1),
					"Two events must merge by minimum; ambiguous/gap/no-window remain unranked and native-zero gets a sketch rank");
				check(!beta || p0, "Panel-local parent IDs and event IDs cannot leak into another panel");
				rows++;
			}
		}
		check(rows==6*97, "All six targets and every method survive reporting, including unrankable targets");
	}

	/** Shared tiny cache producer for the saved-grid and hybrid integration tests. */
	static final class Fixture{
		Fixture(final String root, final String fault_){
			this(root, fault_, null);
		}
		Fixture(final String root, final String fault_, final byte[] padding_){
			fault=fault_;
			padding=padding_;
			dir=root+'/'+fault;
			check(new File(dir).mkdir(), "Each literal saved-cache case needs a fresh directory");
			members();
			write("panels.tsv", NcrnaSavedRankSummary.PANEL_HEADER+'\n'
				+PANEL+"\t"+dir+"/alpha\t"+dir+"/alpha.labels.tsv\t.\t.\t2\t2\t5\t1\t1\t1\n"
				+"beta\t"+dir+"/beta\t"+dir+"/beta.labels.tsv\t.\t.\t1\t1\t1\t1\t0\t0\n");
			write("selection.tsv", INVENTORY+'\n'
				+PANEL+"\t0\tinput_a\t/fixture/a.fa\t10\t5\t1\t1\t1\t"+(fault.equals("inventory") ? "0\t2" : "1\t1")+'\n'
				+PANEL+"\t1\tinput_empty\t/fixture/empty.fa\t11\t0\t0\t0\t0\t0\t0\n"
				+"beta\t0\tinput_b\t/fixture/b.fa\t12\t1\t1\t0\t0\t0\t0\n");
			panel(false);
			panel(true);
		}
		RankSavedNcrnaWindows.Config config(final String output, final int threads){
			final RankSavedNcrnaWindows.Config c=new RankSavedNcrnaWindows.Config();
			c.panels=dir+"/panels.tsv";
			c.selection=dir+"/selection.tsv";
			c.output=dir+'/'+output;
			c.members=dir+"/members.fa";
			c.ledger=dir+"/ledger.tsv";
			c.assignments=dir+"/assignments.tsv";
			c.census=dir+"/models.tsv";
			c.summary=dir+"/members.summary.tsv";
			c.references=dir+"/ref.fa";
			c.models=2;
			c.groups=3;
			c.threads=threads;
			return c;
		}
		private void members(){
			final ByteBuilder ledger=new ByteBuilder("locus_id\tsource_family\tmember_family\tbuild_group\trouting_rule\trouting_bits\tlength\tinput_header\n");
			final ByteBuilder assignments=new ByteBuilder("#fullRecordId\tsourceTid\tbuilderGroup\tclusterKey\tmodelLabel\tstatus\n");
			final ByteBuilder members=new ByteBuilder(), references=new ByteBuilder();
			final ByteBuilder census=new ByteBuilder(NcrnaMemberLengthBounds.HEADER+'\n');
			for(int model=0; model<2; model++){
				final String name="cobalamin_RF00174_fixture"+model;
				final String raw="RF00174_"+model+" source=tid_1_species status=ANNOTATION_FULL_COPY";
				final String header=raw+" member_family=RF00174 build_group=RF00174";
				ledger.append("RF00174_").append(model).append("\tRF00174\tRF00174\tRF00174\tACCEPTED_RF00174_PRODUCT\tNA\t20\t").append(raw).nl();
				assignments.append(header).append("\t1\tunknown\tRF00174_fixture\t").append(name).append("\tkept\n");
				members.append('>').append(header).nl().append(model==0 ? A : G).nl();
				references.append('>').append(name).append(" n=1\n").append(model==0 ? A : G).nl();
				census.append(name).append("\tRF00174\t20\t1\t20\t20\t20\t20\t20\t1\tPURE\t1\t1\t0\t0\t0\t0\t0\t0\n");
			}
			write("ledger.tsv", ledger.toString());
			write("assignments.tsv", assignments.toString());
			write("members.fa", members.toString());
			write("ref.fa", references.toString());
			write("models.tsv", census.toString());
			write("members.summary.tsv", "emitted_family\tinput_members\toriginal_kept\toriginal_small\toriginal_orphan\tadded_memberships\tconstruction_memberships\tmodels\tpure\tmixed\tempty\tpurity_cut\tgroup_min\tgroup_p10\tgroup_median\tgroup_p90\tgroup_max\n"
				+"RF00174\t2\t2\t0\t0\t0\t2\t2\t2\t0\t0\t0.8\t20\t20\t20\t20\t20\n");
		}
		private void panel(final boolean beta){
			final String prefix=beta ? "beta" : "alpha", id=beta ? "input_b" : "input_a", fasta=beta ? "/fixture/b.fa" : "/fixture/a.fa";
			write(prefix+".inputs.tsv", NcrnaDdlTrainingOverlap.SAMPLE_HEADER+'\n'
				+"0\t"+id+'\t'+fasta+"\t1\t1000\t"+(beta ? "1\t0\t1" : "10\t0\t5")+"\t0\t0\t0\t0\t0\n"
				+(beta ? "" : "1\tinput_empty\t/fixture/empty.fa\t1\t1000\t0\t0\t0\t0\t0\t0\t0\t0\n"));
			final ByteBuilder parents=new ByteBuilder(NcrnaSavedRankSummary.PARENT_HEADER+'\n');
			final ByteBuilder labels=new ByteBuilder(NcrnaSavedRankSummary.STRATA_HEADER+'\n');
			final ByteBuilder loci=new ByteBuilder(NcrnaSavedRankSummary.LOCUS_HEADER+'\n');
			for(int parent=0; parent<(beta ? 1 : 5); parent++){
				final int start=parent==0 ? 1 : (parent-1)*100+101;
				parents.append("0\t").append(id).tab().append(fasta).append("\tp").append(parent).append("\tcontig\t")
					.append(start).tab().append(start+19).append("\t+\tfalse\ttrue\ttrue\tfalse\tRF00174\t.\n");
				labels.append("0\t").append(id).append("\tp").append(parent).tab()
					.append(beta && padding==null || parent==1 ? "SEEN" : parent==2 ? "UNKNOWN" : "NOT_SEEN").nl();
				final int count=beta ? 1 : parent==0 ? 7 : parent==2 ? 0 : 1;
				final boolean unique=parent==0 || parent==4;
				for(int method=0; method<148; method++){
					final int rank=parent==0 ? 1 : parent==4 && method>0 ? 2 : -1;
					final String reason=parent==2 ? "NO_CAPTURED_WINDOW" : parent==3 ? "NO_COMPATIBLE_MODEL"
						: parent==1 ? "AMBIGUOUS_SPANS_ONLY" : rank<0 ? "ZERO_NATIVE_INDEX_SUPPORT" : "RANKED";
					loci.append("0\t").append(id).append("\tp").append(parent).tab().append(method).tab().append(count)
						.tab().append(parent==3 ? 0 : count).tab().append(parent==1 ? 1 : 0).tab().append(unique ? count : 0)
						.tab().append(rank).tab().append(reason).append("\tfalse\tfalse\n");
				}
			}
			// Another scored subtype keeps the target's original two-parent ambiguity.
			if(!beta){parents.append("0\tinput_a\t/fixture/a.fa\tother_subtype\tcontig\t101\t120\t+\tfalse\tfalse\ttrue\tfalse\tRF01482\t.\n");}
			write(prefix+".parents.tsv.gz", parents.toString());
			write(prefix+".labels.tsv", labels.toString());
			write(prefix+".loci.tsv.gz", loci.toString());
			final ByteBuilder windows=new ByteBuilder(NcrnaDdlLocusOverlap.WINDOW_HEADER+'\n');
			final ByteBuilder ranks=new ByteBuilder(RANKS+'\n'), links=new ByteBuilder(LINKS+'\n'), sequences=new ByteBuilder();
			final ByteBuilder oracle=new ByteBuilder(ORACLE+'\n');
			// Eleven windows exceed the worker's eight-window claim block, exercising
			// cross-worker minima for the same parent and an independent second panel.
			for(int event=0; event<(beta ? 1 : 10); event++){
				final int length=20+(padding==null ? 0 : padding.length);
				final int parent=beta || event<2 || event>=5 ? 0 : event==2 ? 1 : event;
				final int start=parent==0 ? 0 : (parent-1)*100+100;
				final boolean unique=parent==0 || parent==4, supported=parent!=3;
				final boolean a=beta || event==1 || event==3;
				final int ordinal=!beta && event==0 && fault.equals("event") ? 1 : 0;
				final String name="0_"+event;
				final String lead="0\t"+id+'\t'+event+"\t0\tcontig\t";
				for(int model=0; model<2; model++){
					if(beta && model==1 && fault.equals("native_missing")){continue;}
					final int count=!beta && event==0 && model==0 && fault.equals("native_count") ? 65536
						: !beta && event==0 && model==0 && fault.equals("native_tuple") ? 1 : (a==(model==0) ? 12 : 0);
					final boolean accepted=model==0 && supported;
					final String modelName=!beta && event==0 && model==0 && fault.equals("native_name")
						? "cobalamin_RF00174_fixture0 n=2" : "cobalamin_RF00174_fixture"+model+" n=1";
					final String row=lead+model+'\t'+modelName+'\t'+count+'\t'+(count>0 ? 1 : 0)+'\t'+accepted+'\t'
						+(accepted ? start+1 : -1)+'\t'+(accepted ? start+20 : -1)+"\t0\t"+(accepted ? "1" : ".")+'\t'
						+(accepted ? unique ? 1 : 2 : 0)+'\n';
					oracle.append(row);
					if(!beta && event==0 && model==0 && fault.equals("native_duplicate")){oracle.append(row);}
				}
				windows.append(lead).append("contig\t").append(!beta && event==0 && fault.equals("strand") ? 1 : 0)
					.append("\t0\t").append(start).tab().append(start+length-1).tab().append(start).tab().append(start+length-1)
					.tab().append(1000+(padding==null ? 0 : padding.length)).append("\tFIXTURE\t").append(length).nl();
				if(beta || event!=0 || !fault.equals("missing_sequence")){
					sequences.append('>').append(name).nl().append(a || padding!=null && parent==0 ? A : G);
					if(padding!=null){sequences.append(padding);}
					sequences.nl();
				}
				if(!beta && event==0 && fault.equals("duplicate_sequence")){sequences.append('>').append(name).nl().append(G).nl();}
				for(int method=0; method<148; method++){
					if(!beta && event==0 && method==147 && fault.equals("missing_method")){continue;}
					final int rank=!unique || method==0 && !a ? -1 : a ? 1 : 2;
					final String row="0\t"+(!beta && event==0 && fault.equals("slot") ? "wrong_id" : id)+'\t'+event+'\t'+ordinal
						+"\tcontig\tp"+parent+'\t'+method+'\t'+length+"\t20\t"+(supported ? 1 : 0)+'\t'+(unique ? 1 : 0)+'\t'+rank+'\t'
						+(!beta && event==0 && method==52 && fault.equals("tuple") ? 1 : rank)+'\t'+rank+'\t'+(rank<0 ? -1 : 0)+'\t'
						+(!supported ? "NO_COMPATIBLE_MODEL" : !unique ? "AMBIGUOUS_SPANS_ONLY" : rank<0 ? "ZERO_NATIVE_INDEX_SUPPORT" : "RANKED")+'\n';
					ranks.append(row);
					if(!beta && event==0 && method==147 && fault.equals("duplicate_method")){ranks.append(row);}
				}
				if(supported){
					links.append(lead).append(!beta && event==0 && fault.equals("link") ? 1 : 0).append("\tp").append(parent)
						.append("\ttrue\t").append(unique ? 1 : 2).tab().append(!beta && event==2 && fault.equals("unique") || unique).nl();
					if(!unique){links.append(lead).append("0\tother_subtype\tfalse\t2\tfalse\n");}
				}
			}
			write(prefix+".windows.tsv.gz", windows.toString());
			write(prefix+".ranks.tsv.gz", ranks.toString());
			write(prefix+".links.tsv.gz", links.toString());
			write(prefix+".windows.fa.gz", sequences.toString());
			write(prefix+".oracle.tsv.gz", oracle.toString());
		}
		private void write(final String path, final String value){
			try(Output out=new Output(dir+'/'+path)){out.writer.print(value.getBytes(StandardCharsets.UTF_8));}
		}
		final String dir, fault;
		final byte[] padding;
	}

	private static void equal(final String a, final String b){
		try(Input left=new Input(a); Input right=new Input(b)){
			for(byte[] row=left.reader.nextLine(); row!=null; row=left.reader.nextLine()){
				check(Arrays.equals(row, right.reader.nextLine()), "Worker count cannot change saved rank evidence: "+a);
			}
			check(right.reader.nextLine()==null, "Worker count cannot add evidence rows");
		}
	}
	private static String text(final byte[] row){
		check(row!=null, "Literal fixture output cannot end before the expected row");
		return new String(row, StandardCharsets.UTF_8);
	}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
	private static final String A="AAAAAAAAAAAAAAAAAAAA", G="GGGGGGGGGGGGGGGGGGGG";
	private static final String PANEL="α-panel";
	private static final byte[] PANEL_BYTES=PANEL.getBytes(StandardCharsets.UTF_8);
	private static final String INVENTORY="panel\tinput_slot\tinput_id\tfasta\tresolved_assembly_group\ttarget_loci\tcommon_loci\tnative_zero_support\tno_compatible_model\tno_window\tambiguous_only";
	private static final String RANKS="input_slot\tinput_id\tevent\tcontig_ordinal\tseqid\tparent\tmethod\twindow_length\tlocus_length\tcompatible_models\tunique_models\trank\ttie_first\ttie_last\tmodel_index\treason";
	private static final String LINKS="input_slot\tinput_id\tevent\tcontig_ordinal\tseqid\tmodel_index\tparent\ttarget\tspan_parent_links\tuniquely_attributed";
	private static final String ORACLE="input_slot\tinput_id\tevent\tcontig_ordinal\tseqid\tmodel_index\tmodel\tnative_shared_keys\tnative_shortlist_rank\tfully_accepted\taccepted_start1\taccepted_stop1\tpartial\tscore\tscored_parent_links";
}
