package test.prok;

import java.io.File;
import java.util.Arrays;
import java.util.Locale;

import fileIO.ByteStreamWriter;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaFeatureTrialReport;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Hand-derived paired seed deltas and refusal of incomplete/duplicated reports.
 * @author Brian Bushnell, Ganyu
 */
public final class NcrnaFeatureTrialReportTest{
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaFeatureTrialReportTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String a=pair[0].toLowerCase(Locale.ROOT);
				final String b=pair.length>1 ? pair[1] : null;
				check(parser.parse(arg, a, b), "Unknown feature-report fixture argument");
			}
			check(parser.out1!=null && parser.in1==null && parser.in2==null && parser.out2==null
				&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"Supply out=fresh_fixture_directory and native runtime options");
			final File root=new File(parser.out1);
			check(root.mkdir(), "Fixture directory must be fresh");
			run(root);
			pp.outstream.println("NCRNA_FEATURE_TRIAL_REPORT_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}
	private static void run(final File root){
		final ByteBuilder manifest=new ByteBuilder("variant\tseed\tend\tsummary\n");
		final ByteBuilder conditional=new ByteBuilder("variant\tseed\tend\tsummary\n");
		String lengthFirst=null;
		for(final String variant:new String[]{"baseline", "length", "gc", "identity"}){
			for(int seed=58; seed<=60; seed++){
				for(final String end:new String[]{"5prime", "3prime"}){
					if(variant.equals("identity") && (seed!=58 || !end.equals("5prime"))){continue;}
					final String path=path(root, variant+'_'+seed+'_'+end+".tsv");
					final int exact=variant.equals("length") ? 61-seed : variant.equals("gc") ? 3 : 2;
					final double rms=variant.equals("length") ? 1.8+(seed-58)*.1 : variant.equals("gc") ? 1.7 : 2;
					String saved=summary(end, exact, rms);
					if(variant.equals("identity")){
						saved=summary(end, 1, 3).replace("\t4\t66.66666667\t", "\t2\t33.33333333\t")
							.replace("\t2\t1\t3\t0.66666667", "\t0\t0\t6\tNA");
					}
					write(path, saved);
					write(path+".conditional.tsv", saved.replace("better_fraction\n", "better_fraction_if_changed\n").replace("\tNA\n", "\t\n"));
					manifest.append(variant).tab().append(seed).tab().append(end).tab().append(path).nl();
					conditional.append(variant).tab().append(seed).tab().append(end).tab().append(path).append(".conditional.tsv\n");
					if(variant.equals("length") && seed==58 && end.equals("5prime")){lengthFirst=path;}
				}
			}
		}
		final String input=path(root, "manifest.tsv");
		write(input, manifest.toString());
		new NcrnaFeatureTrialReport(new String[]{"in="+input, "out="+path(root, "good"), "expected=6,6"}).process();
		write(path(root, "conditional.tsv"), conditional.toString());
		new NcrnaFeatureTrialReport(new String[]{"in="+path(root, "conditional.tsv"), "out="+path(root, "conditional"), "expected=6,6"}).process();
		for(final String suffix:new String[]{".all.tsv", ".paired.tsv", ".spread.tsv"}){
			try(Input legacy=new Input(path(root, "good"+suffix)); Input current=new Input(path(root, "conditional"+suffix))){
				for(byte[] row; (row=legacy.reader.nextLine())!=null;){
					check(Arrays.equals(row, current.reader.nextLine()), "Old/new fraction schemas must produce identical normalized reports");
					final LineParser1 p=new LineParser1('\t').set(row);
					for(int i=0; i<p.terms(); i++){
						check(!p.termEquals("NA", i) && !p.termEquals("NaN", i) && !p.termEquals("Infinity", i),
							"No-change fractions must be omitted, and measured results must remain finite");
					}
				}
				check(current.reader.nextLine()==null, "Normalization cannot add scientific rows");
			}
		}
		int count=0;
		try(Input file=new Input(path(root, "good.spread.tsv"))){
			file.reader.nextLine();
			final LineParser1 p=new LineParser1('\t');
			for(byte[] line=file.reader.nextLine(); line!=null; line=file.reader.nextLine()){
				p.set(line);
				check(p.terms()==12, "Seed spread has the fixed descriptive schema");
				if(p.termEquals("length", 0)){
					check(p.parseInt(3)==3 && close(p.parseDouble(5), 0) && close(p.parseDouble(6), -100/6.0)
						&& close(p.parseDouble(7), 100/6.0) && close(p.parseDouble(8), -.1)
						&& close(p.parseDouble(9), -.2) && close(p.parseDouble(10), 0),
						"Length exact deltas +1,0,-1 cancel; RMS deltas -.2,-.1,0 average -.1");
				}else if(p.termEquals("gc", 0)){
					check(p.parseInt(3)==3 && close(p.parseDouble(5), 100/6.0) && close(p.parseDouble(8), -.3),
						"Three paired GC seeds each gain one exact trial and reduce RMS .3");
				}else{
					check(p.termEquals("identity", 0) && p.termEquals("5prime", 1) && p.parseInt(3)==1
						&& close(p.parseDouble(5), -100/6.0) && close(p.parseDouble(8), 1), "The unchanged identity canary stays one seed with raw errors, one fewer exact trial and one RMS base above baseline");
				}
				count++;
			}
		}
		check(count==10, "Two policies for each length/GC end and the sole identity start");
		try(Input file=new Input(path(root, "good.paired.tsv"))){
			file.reader.nextLine();
			final LineParser1 p=new LineParser1('\t').set(file.reader.nextLine());
			check(p.termEquals("length", 0) && p.parseInt(1)==58 && p.parseInt(5)==3 && p.parseInt(6)==2
				&& p.parseInt(7)==1 && close(p.parseDouble(9), 100/6.0) && close(p.parseDouble(16), -.2),
				"Per-seed deltas use the paired baseline, preserve counts and state sign in units");
			int omitted=0;
			for(byte[] row; (row=file.reader.nextLine())!=null;){
				p.set(row);
				if(p.termEquals("identity", 0)){
					check(p.parseLong(17)==0 && p.parseLong(19)==0 && p.parseLong(21)==6 && p.length(23)==0
						&& close(p.parseDouble(24), 2.0/3), "Only the no-change member fraction is absent; the measured baseline fraction remains");
					omitted++;
				}
			}
			check(omitted==2, "Both one-seed identity policies exercise conditional fraction omission");
		}
		final String all=manifest.toString();
		write(input, all.substring(0, all.lastIndexOf("identity\t")));
		fails(input, path(root, "missing"));
		write(input, all+all.substring(all.lastIndexOf("identity\t")));
		fails(input, path(root, "duplicate"));
		write(input, all);
		write(lengthFirst, summary("5prime", 3, 1.8).replace("do_nothing\t5prime\t6\t1\t16.66666667", "do_nothing\t5prime\t6\t2\t33.33333333"));
		fails(input, path(root, "changed_baseline"));
		write(lengthFirst, summary("5prime", 3, 1.8).replace("0.66666667", "0.5"));
		fails(input, path(root, "bad_fraction"));
		check(!new File(path(root, "bad_fraction.all.tsv")).exists(), "Malformed inputs fail before any report publication");
		write(lengthFirst, summary("5prime", 3, 1.8).replace("\t0\t0\t6\tNA", "\t0\t0\t6\t0"));
		fails(input, path(root, "invented_fraction"));
		for(final String value:new String[]{"NaN", "Infinity"}){
			write(lengthFirst, summary("5prime", 3, 1.8).replace("\t1.80000000\t", "\t"+value+"\t"));
			fails(input, path(root, "bad_rms_"+value));
		}
	}
	private static String summary(final String end, final int exact, final double rms){
		final ByteBuilder b=new ByteBuilder("policy\tend\tusable_jitter_trials\texact\texact_pct\twithin3\twithin3_pct\trms_after_bp\tbetter\tworse\ttied\tbetter_fraction\n");
		b.append("do_nothing\t").append(end).append("\t6\t1\t16.66666667\t2\t33.33333333\t3\t0\t0\t6\tNA\n");
		for(final String policy:new String[]{"margin0", "margin0.1"}){
			b.append(policy).tab().append(end).append("\t6\t").append(exact).tab().append(100*exact/6.0, 8)
				.append("\t4\t66.66666667\t").append(rms, 8).append("\t2\t1\t3\t0.66666667\n");
		}
		return b.toString();
	}
	private static void fails(final String input, final String output){
		try{new NcrnaFeatureTrialReport(new String[]{"in="+input, "out="+output, "expected=6,6"}).process();}
		catch(IllegalArgumentException expected){return;}
		throw new AssertionError("Incomplete or contradictory saved summaries must be rejected");
	}
	private static boolean close(final double a, final double b){return Math.abs(a-b)<1e-7;}
	private static String path(final File root, final String name){return new File(root, name).getPath();}
	private static void write(final String path, final String text){
		try(Output output=new Output(new ByteStreamWriter(path, true, false, false))){output.writer.print(text);}
	}
	private static void check(final boolean valid, final String message){if(!valid){throw new AssertionError(message);}}
}
