package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.Locale;
import java.util.function.Consumer;

import cardinality.DynamicDemiLog;
import idaligner.Tracer;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaDdlCohorts;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.RankDualNcrnaWindows;
import shared.Shared;
import structures.ByteBuilder;

/** Frozen bins/roles, hierarchical statistics and complete saved-only reporting.
 * Source genome paths are nonexistent; the cache is a literal synthetic witness.
 * @author Brian, Raiden
 */
public final class RankDualNcrnaWindowsTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, RankDualNcrnaWindowsTest.class, false);
		final int previous=DynamicDemiLog.exponentBits();
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown saved dual fixture parameter: "+arg);
			}
			boolean assertions=false;
			assert(assertions=true);
			check(assertions && parser.out1!=null && parser.out2==null && parser.in1==null && parser.in2==null
				&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"Use -ea and out=fresh_fixture_directory without sampling or external inputs");
			final File root=new File(parser.out1).getAbsoluteFile();
			check(root.mkdir(), "Saved dual fixtures cannot overwrite existing evidence");
			DynamicDemiLog.setExponent(4);
			NcrnaDdlMomentsTest.main(new String[]{"run=t", "t=1"});
			literal();
			integration(root.getPath());
			NcrnaDdlLearnedWeightFixture.run(root.getPath());
			pp.outstream.println("NCRNA_DUAL_SAVED_RANKS_TEST_PASS");
		}finally{
			DynamicDemiLog.setExponent(previous);
			Shared.closeStream(pp.outstream);
		}
	}

	private static void literal(){
		check(RankDualNcrnaWindows.invalidBase("ACGTNacgtn".getBytes(StandardCharsets.US_ASCII))==-1
			&& RankDualNcrnaWindows.invalidBase("ACGR".getBytes(StandardCharsets.US_ASCII))==3,
			"Unsupported IUPAC symbols cannot silently enter the shared ACGTN word encoder");
		final float[] q=NcrnaDdlCohorts.quantiles(new float[]{1,.5f,.25f,.25f});
		check(Arrays.equals(q, new float[]{.25f,.25f,.5f}) && NcrnaDdlCohorts.bin(.25f, q)==0,
			"Nearest-rank quartiles include every tied identity at the first-quartile boundary");
		check(Float.isNaN(NcrnaDdlCohorts.quantiles(new float[0])[0]) && NcrnaDdlCohorts.bin(Float.NaN, q)==4,
			"Missing identity or an empty cutoff population cannot become the low-identity class");
		final int[] groups={2,4,6,8,10};
		final int[] sides=NcrnaDdlCohorts.split(groups);
		int a=0;
		for(final int side:sides){
			check(side==0 || side==1, "Each original group belongs to exactly one half");
			if(side==0){a++;}
		}
		check(a==2 && Arrays.equals(groups, new int[]{2,4,6,8,10}), "The seeded split preserves the group roster and assigns floor(n/2) groups to A");
		check(NcrnaDdlCohorts.role(1, 1, 0)==NcrnaDdlCohorts.EXCLUDED_SEEN
			&& NcrnaDdlCohorts.role(1, -1, 0)==NcrnaDdlCohorts.AUXILIARY_SEEN
			&& NcrnaDdlCohorts.role(0, 1, 0)==NcrnaDdlCohorts.JUDGE
			&& NcrnaDdlCohorts.role(0, 1, 1)==NcrnaDdlCohorts.TRAIN,
			"Judging groups never fit through a SEEN parent; direction swapping changes their role without splitting a group");
		final double[] event=new double[8];
		final int signal=RankDualNcrnaWindows.normalizedMoments(new char[]{0,4,2,6}, new char[]{2,6,4,8}, 0, 4, new long[]{3}, event);
		check(signal==2 && Arrays.equals(event, new double[]{2,8,4,20,4,20,6,40}),
			"Every model has equal weight inside its own event/class; second moments retain model variation");
		check(RankDualNcrnaWindows.best(new char[]{0,1}, new char[]{1,0}, 0, 2, 1+1e-12, new long[]{2})==1
			&& RankDualNcrnaWindows.best(new char[]{0,1}, new char[]{1,0}, 0, 2, 1, new long[]{2})==2,
			"Empirical weights keep double precision; exact ties retain original model order");
		check(RankDualNcrnaWindows.empiricalWeight(7,1,2,5,1,1)==.75
			&& RankDualNcrnaWindows.empiricalWeight(1,2,1,3,1,1)==0
			&& Double.isNaN(RankDualNcrnaWindows.empiricalWeight(7,1,0,5,1,1))
			&& Double.isNaN(RankDualNcrnaWindows.empiricalWeight(7,1,2,1,2,1)),
			"Weight is the ratio of standardized separations; missing spread/nonpositive long separation is unavailable, not clipped");
	}

	private static void integration(final String root){
		final RankSavedNcrnaWindowsTest.Fixture f=new RankSavedNcrnaWindowsTest.Fixture(root, "dual");
		coherent(f);
		baseline(f, "baseline.tsv.gz", false);
		RankDualNcrnaWindows.run(config(f, "one", 1));
		RankDualNcrnaWindows.run(config(f, "two", 2));
		for(final String suffix:SUFFIXES){equal(f.dir+"/one"+suffix, f.dir+"/two"+suffix);}
		DynamicDemiLog.setExponent(6);
		final boolean global=Tracer.GLOBAL;
		Tracer.GLOBAL=true;
		try{
			RankDualNcrnaWindows.main(new String[]{"panels="+f.dir+"/panels.tsv", "selection="+f.dir+"/selection.tsv",
				"baseline="+f.dir+"/baseline.tsv.gz", "identities="+f.dir+"/one", "members="+f.dir+"/members.fa", "ledger="+f.dir+"/ledger.tsv",
				"assignments="+f.dir+"/assignments.tsv", "census="+f.dir+"/models.tsv", "summary="+f.dir+"/members.summary.tsv",
				"ref="+f.dir+"/ref.fa", "out="+f.dir+"/cli", "expectedmodels=2", "expectedgroups=3", "t=2"});
		}finally{Tracer.GLOBAL=global;}
		check(DynamicDemiLog.exponentBits()==6, "The actual CLI must set its native exponent and restore its embedding caller's previous setting");
		DynamicDemiLog.setExponent(4);
		for(final String suffix:SUFFIXES){equal(f.dir+"/one"+suffix, f.dir+"/cli"+suffix);}
		RankCompactNcrnaWindowsTest.run(f);
		savedFault(f, "parent", ".identity.tsv", rows->{change(rows, 2, "wrong_input");}, "parent binding");
		savedFault(f, "model", ".identity.tsv", rows->{change(rows, 8, "wrong_model");}, "model/event/gene witnesses");
		savedFault(f, "event", ".identity.tsv", rows->{change(rows, 9, "999999");}, "model/event/gene witnesses");
		savedFault(f, "value", ".identity.tsv", rows->{change(rows, 13, "0.123456789");}, "decimal disagrees");
		savedFault(f, "bin", ".identity.tsv", rows->{change(rows, 15, "Q4");}, "bins/group sides differ");
		savedFault(f, "status", ".identity.tsv", rows->{change(rows, 17, "NOT_NATIVE_RANKABLE");}, "status disagrees");
		savedFault(f, "short", ".identity.tsv", rows->{rows.remove(rows.size()-1);}, "omit an original parent");
		savedFault(f, "extra", ".identity.tsv", rows->{rows.add(rows.get(1));}, "extra or duplicate original parent");
		savedFault(f, "group", ".groups.tsv", rows->{
			final String[] fields=rows.get(1).split("\t", -1);
			change(rows, 1, fields[1].equals("A") ? "B" : "A");
		}, "roster/side/seed differs");
		savedFault(f, "cut", ".identity_bins.tsv", rows->{
			final String[] fields=rows.get(1).split("\t", -1);
			change(rows, 16, Integer.toString(Integer.parseInt(fields[16])+1));
		}, "cutoff population/availability/float bits differ");
		final LineParser1 p=new LineParser1('\t');
		int methods=0, inputRows=0, fitRows=0;
		try(Input in=new Input(f.dir+"/one.methods.tsv")){
			check(in.reader.nextLine()!=null, "Methods require their complete geometry/weight schema");
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){methods++;}
		}
		check(methods==67, "Keep three references, eight competing methods and eight times seven partition methods");
		try(Input in=new Input(f.dir+"/one.inputs.tsv")){
			check(in.reader.nextLine()!=null, "The original selected-input inventory requires a header");
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){inputRows++;}
		}
		check(inputRows==3, "The selected zero-target input cannot disappear from reporting");
		try(Input in=new Input(f.dir+"/one.statistics_loci.tsv.gz")){
			check(in.reader.nextLine()!=null, "Every partition retains event/model counts and fit eligibility per parent");
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				if(!p.termEquals("p0", 4)){continue;}
				final boolean alpha=Arrays.equals(p.parseByteArray(1), ALPHA);
				check(p.parseInt(7)==(alpha ? 7 : 1) && p.termEquals(alpha ? "false" : "true", 12) && p.termEquals("true", 13),
					"Fit A excludes the judging NOT_SEEN group; fit B includes it; outside-roster SEEN auxiliary remains eligible");
				fitRows++;
			}
		}
		check(fitRows==16, "Both common fixture parents must retain all eight partition fit records");
		try(Input in=new Input(f.dir+"/one.loci.tsv.gz")){
			check(in.reader.nextLine()!=null, "Every method/parent retains original identity and disposition");
			int rows=0;
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){rows++;}
			check(rows==67*6, "Every target, including ambiguous/missing/native-zero parents, survives every method");
		}
		baseline(f, "bad.tsv.gz", true);
		final RankDualNcrnaWindows.Config bad=config(f, "rejected", 1);
		bad.baseline=f.dir+"/bad.tsv.gz";
		try{RankDualNcrnaWindows.run(bad);}
		catch(IllegalArgumentException invalid){
			check(invalid.getMessage().contains("Native minima") && !new File(f.dir+"/rejected.identity.tsv").exists(),
				"An exchanged native baseline must fail before identity measurement or rank publication");
			return;
		}
		throw new AssertionError("A changed original native minimum was accepted");
	}

	/** Small literal fixture tables only; the product reader never uses split or
	 * retains field Strings. Each fault must fail at its intended binding check. */
	private static void savedFault(final RankSavedNcrnaWindowsTest.Fixture f, final String tag,
			final String target, final Consumer<ArrayList<String>> mutation, final String reason){
		final String prefix=f.dir+"/bad_saved_"+tag;
		for(final String suffix:new String[]{".identity.tsv", ".identity_bins.tsv", ".groups.tsv"}){
			final ArrayList<String> rows=new ArrayList<String>();
			try(Input in=new Input(f.dir+"/one"+suffix)){
				for(byte[] row; (row=in.reader.nextLine())!=null;){rows.add(new String(row, StandardCharsets.UTF_8));}
			}
			if(suffix.equals(target)){mutation.accept(rows);}
			try(Output out=new Output(prefix+suffix)){
				for(final String row:rows){out.writer.println(row.getBytes(StandardCharsets.UTF_8));}
			}
		}
		final RankDualNcrnaWindows.Config c=config(f, "reject_saved_"+tag, 1);
		c.identities=prefix;
		try{RankDualNcrnaWindows.run(c);}
		catch(IllegalArgumentException expected){
			check(expected.getMessage().contains(reason) && !new File(c.output+".identity.tsv").exists(),
				"Malformed frozen evidence must fail before any new report,at its intended check: "+tag+" actual="+expected.getMessage());
			return;
		}
		throw new AssertionError("Malformed frozen identity evidence was accepted: "+tag);
	}

	private static void change(final ArrayList<String> rows, final int field, final String value){
		check(rows.size()>1, "Each mutation needs the original first data row");
		final String[] fields=rows.get(1).split("\t", -1);
		fields[field]=value;
		rows.set(1, String.join("\t", fields));
	}

	/** Preserve the old toy ranking cache before making its gene witnesses agree. */
	static void coherent(final RankSavedNcrnaWindowsTest.Fixture f){
		final String path=f.dir+"/alpha.windows.fa.gz", original=f.dir+"/alpha.original.windows.fa.gz";
		check(new File(path).renameTo(new File(original)), "Archive the toy FASTA before changing its synthetic gene witnesses");
		try(Input in=new Input(original); Output out=new Output(path)){
			boolean p0=false;
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				if(row.length>0 && row[0]=='>'){
					final int event=Parse.parseInt(row, 3, row.length);
					p0=event<2 || event>=5;
				}else if(p0){Arrays.fill(row, (byte)'A');}
				out.writer.println(row);
			}
		}
	}

	private static void baseline(final RankSavedNcrnaWindowsTest.Fixture f, final String file, final boolean bad){
		baseline(f, file, bad, false);
	}

	static void baseline(final RankSavedNcrnaWindowsTest.Fixture f, final String file, final boolean bad, final boolean learned){
		try(Output out=new Output(f.dir+'/'+file)){
			out.writer.println("panel\tinput_slot\tinput_id\tparent\tseen_stratum\tmethod\told_rank\tquery_fixed_rank\tall_fixed_rank\tsource_disposition");
			final ByteBuilder b=new ByteBuilder();
			for(int panel=0; panel<2; panel++){
				for(int parent=0; parent<(panel==0 ? 5 : 1); parent++){
					final String state=panel==1 || parent==1 ? "SEEN" : parent==2 ? "UNKNOWN" : "NOT_SEEN";
					final String reason=parent==0 || learned && parent==1 ? "RANKED" : parent==1 ? "AMBIGUOUS_SPANS_ONLY" : parent==2 ? "NO_CAPTURED_WINDOW"
						: parent==3 ? "NO_COMPATIBLE_MODEL" : "ZERO_NATIVE_INDEX_SUPPORT";
					for(int method=0; method<19; method++){
						final int rank=bad && panel==0 && parent==0 && method==0 ? 2 : parent==0 || learned && parent==1 ? 1 : parent==4 && method>0 ? 2 : -1;
						b.clear().append(panel==0 ? ALPHA : BETA)
							.append("\t0\t").append(panel==0 ? "input_a" : "input_b").append("\tp").append(parent).tab().append(state).tab()
							.append(method<3 ? REFERENCES[method] : "canonical_control_"+CONTROLS[method-3])
							.tab().append(rank).tab().append(rank).tab().append(rank).tab().append(reason).nl();
						out.writer.print(b);
					}
				}
			}
		}
	}

	static RankDualNcrnaWindows.Config config(final RankSavedNcrnaWindowsTest.Fixture f, final String output, final int threads){
		final RankDualNcrnaWindows.Config c=new RankDualNcrnaWindows.Config();
		c.panels=f.dir+"/panels.tsv";
		c.selection=f.dir+"/selection.tsv";
		c.baseline=f.dir+"/baseline.tsv.gz";
		c.members=f.dir+"/members.fa";
		c.ledger=f.dir+"/ledger.tsv";
		c.assignments=f.dir+"/assignments.tsv";
		c.census=f.dir+"/models.tsv";
		c.summary=f.dir+"/members.summary.tsv";
		c.references=f.dir+"/ref.fa";
		c.output=f.dir+'/'+output;
		c.models=2;
		c.groups=3;
		c.threads=threads;
		return c;
	}

	private static void equal(final String a, final String b){
		try(Input left=new Input(a); Input right=new Input(b)){
			for(byte[] row=left.reader.nextLine(); row!=null; row=left.reader.nextLine()){
				check(Arrays.equals(row, right.reader.nextLine()), "Worker count changed a deterministic saved report: "+a);
			}
			check(right.reader.nextLine()==null, "An extra worker cannot append duplicate saved report rows");
		}
	}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
	private static final int[] CONTROLS={52,54,55,57,58,60,61,63,73,75,76,78,79,81,82,84};
	private static final String[] REFERENCES={"native", "masked10_raw", "hybrid"};
	private static final byte[] ALPHA="α-panel".getBytes(StandardCharsets.UTF_8), BETA="beta".getBytes(StandardCharsets.UTF_8);
	private static final String[] SUFFIXES={".identity.tsv", ".identity_bins.tsv", ".groups.tsv", ".models.tsv", ".methods.tsv", ".loci.tsv.gz", ".summary.tsv", ".last_losses.tsv", ".halves.tsv.gz", ".moments.tsv", ".statistics_loci.tsv.gz", ".weights.tsv", ".inputs.tsv", ".work.tsv"};
}
