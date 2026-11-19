package test.prok;

import java.io.File;
import java.util.Locale;

import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaModelEndpointOffsets;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;

/** Literal median/sign/support/geometry/full-error controls, without fitting.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaModelEndpointOffsetsTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaModelEndpointOffsetsTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String a=pair[0].toLowerCase(Locale.ROOT);
				final String b=pair.length>1 ? pair[1] : null;
				check(parser.parse(arg, a, b), "Unknown model-offset fixture argument");
			}
			check(parser.out1!=null && parser.in1==null && parser.in2==null && parser.out2==null
				&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"Supply only out=fresh_directory and native runtime options");
			final File root=new File(parser.out1).getAbsoluteFile();
			check(root.mkdir(), "Finite fixture directory must be fresh");
			run(root);
			compare(root);
			pp.outstream.println("NCRNA_MODEL_ENDPOINT_OFFSETS_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private static void run(final File root){
		final String refA=write(root, "a.fa", ">A n=5\nAAAA\n");
		final String refB=write(root, "b.fa", ">B n=2\nCCCC\n");
		final String refC=write(root, "c.fa", ">C n=0\nGGGG\n");
		final String bindings=write(root, "sources.tsv", "input_ordinal0\tassembly_id\tinput_id\tfasta\tcanonical_new_gff\n"
			+"0\ttrain_assembly\ttrain\t/train.fa\t/train.gff\n"
			+"1\touter_assembly\touter_alias\t/other_alias.fa\t/outer.gff\n");
		final String audit=write(root, "outer.tsv", "# HELD OUT (panel genomes excluded from tables and training): exclusion bindings; not a fitted result\n"
			+"panel\trow1\tfasta\tassembly\tstatus\n"
			+"/new1000\t1\t/outer.fa\touter_assembly\tEXCLUDED_ASSEMBLY\n"
			+"/old300\t1\t/second.fa\tsecond_outer_assembly\tEXCLUDED_ASSEMBLY\n");
		final String trainText=HEADER
			+"train\tt0\tmain\tA n=5\tt0\t102\t202\t100\t199\t+\t2\t3\n"
			+"train\tt1\tmain\tA n=5\tt1\t102\t202\t100\t199\t+\t2\t3\n"
			+"train\tt2\tmain\tA n=5\tt2\t102\t202\t100\t199\t+\t2\t3\n"
			+"train\tt3\tmain\tA n=5\tt3\t104\t204\t100\t199\t+\t4\t5\n"
			+"train\tt4\tmain\tA n=5\tt4\t108\t208\t100\t199\t+\t8\t9\n"
			+"train\tt5\tsecond\tB n=2\tt5\t304\t402\t300\t399\t-\t-3\t-4\n"
			+"train\tt6\tsecond\tB n=2\tt6\t302\t400\t300\t399\t-\t-1\t-2\n"
			+"outer_alias\touter\tmain\tA n=5\touter\t1100\t1299\t1000\t1199\t+\t100\t100\n";
		final String training=write(root, "train.tsv", trainText);
		final String panelText=HEADER
			+"test\te0\tmain\tA n=5\te0\t102\t202\t100\t199\t+\t2\t3\n"
			+"test\te1\tsecond\tB n=2\te1\t304\t402\t300\t399\t-\t-3\t-4\n"
			+"test\te2\tmain\tA n=5\te2\t700\t799\t700\t799\t+\t0\t0\n"
			+"test\te3\tthird\tC n=0\te3\t500\t599\t500\t599\t+\t0\t0\n"
			+"test\te4\tmain\tA n=5\te4\t112\t212\t100\t199\t+\t12\t13\n"
			+"test\te5\tmain\tA n=5\te5\t1\t100\t1\t100\t+\t0\t0\n";
		final String panel=write(root, "panel.tsv", panelText);
		final String[] args={"in="+training, "refs="+refA+','+refB+','+refC,
			"members=main,second,third", "family=clan", "evaluation="+panel+','+panel,
			"panels=new1000,old300", "bindings="+bindings, "exclude="+audit,
			"expectedmodels=3", "expected=8,7,6,6", "thresholds=1,3,10", "out="+path(root, "good")};
		new NcrnaModelEndpointOffsets(args).process();
		final String prefix=path(root, "good");
		LineParser1 p=find(prefix+".models.tsv", "main", "A", "5prime");
		check(p.terms()==10 && p.parseInt(3)==5 && p.parseInt(4)==2 && p.parseInt(5)==2
			&& p.parseInt(6)==4 && p.parseInt(7)==2 && p.parseInt(8)==4,
			"Literal errors 2,2,2,4,8 give median2, IQR2 and four within3; excluded alias never teaches");
		p=find(prefix+".models.tsv", "second", "B", "5prime");
		check(p.parseInt(3)==2 && p.parseInt(4)==-3 && p.parseInt(7)==2,
			"Even support uses the stated lower nearest-rank median, preserving negative signs");
		p=find(prefix+".models.tsv", "third", "C", "3prime");
		check(p.parseInt(3)==0 && p.termEquals("NA", 4) && p.termEquals("NA", 9),
			"A zero-support catalog model remains explicit and has no learned median");
		p=find(prefix+".support.tsv", "clan", null, null);
		check(p.terms()==12 && p.parseInt(1)==3 && p.parseInt(2)==1 && p.parseInt(3)==1
			&& p.parseInt(4)==1 && p.parseInt(9)==0 && p.parseInt(10)==2 && p.parseInt(11)==5,
			"The complete catalog support distribution includes zero, two and five matches");
		p=findResult(prefix+".single.tsv", "clan", 1, "5prime");
		check(p.terms()==28 && p.parseLong(4)==6 && p.parseLong(5)==5 && p.parseLong(6)==4
			&& p.parseLong(7)==1 && p.parseLong(8)==3 && p.parseLong(10)==4
			&& p.parseLong(16)==5 && p.parseLong(18)==6
			&& p.parseLong(23)==3 && p.parseLong(24)==1 && p.parseLong(25)==2
			&& Math.abs(p.parseDouble(20)-Math.sqrt(157/6.0))<1e-7
			&& Math.abs(p.parseDouble(21)-Math.sqrt(104/6.0))<1e-7 && p.parseDouble(22)<0,
			"Correct signs fix plus/minus calls, preserve real tail distances, harm a correct call, and hold invalid coordinates");
		p=findResult(prefix+".both.tsv", "clan", 1, "both");
		check(p.parseLong(8)==3 && p.parseLong(10)==4 && p.parseLong(16)==5 && p.parseLong(18)==6
			&& Math.abs(p.parseDouble(20)-Math.sqrt(351/12.0))<1e-7
			&& Math.abs(p.parseDouble(21)-Math.sqrt(213/12.0))<1e-7,
			"Joint tolerances require both ends of one locus; pooled RMS counts two ends per locus");
		p=findResult(prefix+".single.tsv", "second", 3, "5prime");
		check(p.parseLong(4)==1 && p.parseLong(5)==0 && p.parseLong(6)==0 && p.parseLong(25)==1,
			"The two-match model is unsupported at minimum n3 and remains unshifted");
		p=findResult(prefix+".single.tsv", "clan", 10, "3prime");
		check(p.parseLong(5)==0 && p.parseLong(6)==0 && p.parseLong(25)==6 && p.termEquals("NA", 26),
			"Unsupported models hold, with unavailable improvement fraction rather than fake zero success");
		final String[] breakdown=java.util.Arrays.copyOf(args, args.length+2);
		breakdown[10]="thresholds=3";
		breakdown[11]="out="+path(root, "breakdown");
		breakdown[12]="breakdown=main";
		breakdown[13]="breakdowncounts=4,4";
		new NcrnaModelEndpointOffsets(breakdown).process();
		p=findBreakdown(path(root, "breakdown.breakdown.tsv"), "3_5", "both");
		check(p.termEquals("SUPPORTED", 0) && p.termEquals("1_2", 2) && p.parseLong(7)==3
			&& p.parseLong(11)==1 && p.parseLong(13)==1 && p.parseLong(26)==2 && p.parseLong(27)==1 && p.parseLong(28)==0,
			"Paired bins use max absolute median3 and max IQR2: three valid shifts conserve one exact, two better and one worse");
		p=findBreakdown(path(root, "breakdown.breakdown.tsv"), "0", "both");
		check(p.parseLong(7)==1 && p.parseLong(10)==1 && p.parseLong(11)==1 && p.parseLong(13)==1,
			"Invalid joint intervals apply zero shift and remain separate from the nonzero median bin");
		check(!new File(path(root, "breakdown.models.tsv")).exists(), "Breakdown mode never republishes original result tables");
		final String[] centered=java.util.Arrays.copyOf(args, args.length+1);
		centered[10]="thresholds=3";
		centered[11]="out="+path(root, "centered");
		centered[12]="centeredhistograms=t";
		new NcrnaModelEndpointOffsets(centered).process();
		p=find(path(root, "centered.offsets.tsv"), "main", "A", null);
		check(p.parseInt(2)==5 && p.parseInt(3)==3 && p.parseInt(4)==2 && p.parseInt(5)==3,
			"The compact resource retains ordered model identity, exact support and the applied signed medians");
		p=find(path(root, "centered.offsets.tsv"), "second", "B", null);
		check(p.parseInt(2)==2 && p.parseInt(4)==0 && p.parseInt(5)==0, "Unsupported models keep explicit zero shifts in the resource");
		p=find(path(root, "centered.training.tsv"), "main", null, null);
		check(p.parseLong(1)==5 && p.parseLong(2)==5 && p.parseLong(3)==0,
			"The outer alias is excluded again before residual distributions are counted");
		try(Input input=new Input(path(root, "centered.main.histogram.tsv"))){
			input.reader.nextLine();
			final String[] expected={"main\t5prime\tsigned\t0\t3", "main\t5prime\tsigned\t2\t1", "main\t5prime\tsigned\t6\t1",
				"main\t3prime\tsigned\t0\t3", "main\t3prime\tsigned\t2\t1", "main\t3prime\tsigned\t6\t1"};
			for(final String line:expected){
				check(java.util.Arrays.equals(input.reader.nextLine(), line.getBytes(java.nio.charset.StandardCharsets.US_ASCII)),
					"Residual histogram must preserve literal signed errors after subtracting the training median");
			}
			check(input.reader.nextLine()==null, "Centered histogram has no invented tail mass");
		}
		minimumLengths(root, centered);
		refit(root, training, bindings, refA+','+refB+','+refC, audit, panel);
		write(panel, panelText+panelText.substring(panelText.indexOf('\n')+1));
		args[11]="out="+path(root, "duplicate");
		fails(()->new NcrnaModelEndpointOffsets(args).process());
		write(panel, panelText.replace("\t-3\t-4\n", "\t3\t-4\n"));
		args[11]="out="+path(root, "bad_sign");
		fails(()->new NcrnaModelEndpointOffsets(args).process());
		write(panel, panelText.replace("A n=5", "UNKNOWN"));
		args[11]="out="+path(root, "unknown_model");
		fails(()->new NcrnaModelEndpointOffsets(args).process());
		write(panel, panelText);
		write(training, trainText.substring(0, trainText.lastIndexOf("outer_alias")));
		args[11]="out="+path(root, "truncated");
		fails(()->new NcrnaModelEndpointOffsets(args).process());
		check(!new File(path(root, "truncated.models.tsv")).exists(), "Accepted input-count failure precedes all output publication");
	}

	/** A former outer assembly must teach the product refit, never silently remain excluded. */
	private static void refit(final File root, final String training, final String bindings, final String refs,
			final String audit, final String panel){
		final String prefix=path(root, "refit");
		final String[] args={"refit=t", "centeredhistograms=t", "in="+training, "refs="+refs,
			"members=main,second,third", "family=clan", "bindings="+bindings, "expectedmodels=3",
			"expected=8,8", "thresholds=6", "out="+prefix};
		new NcrnaModelEndpointOffsets(args).process();
		LineParser1 p=find(prefix+".offsets.tsv", "main", "A", null);
		check(p.parseInt(2)==6 && p.parseInt(3)==6 && p.parseInt(4)==2 && p.parseInt(5)==3,
			"Including the former outer assembly raises main-model support from five to six and enables its learned offsets");
		p=find(prefix+".models.tsv", "main", "A", "5prime");
		check(p.parseInt(3)==6 && p.parseInt(6)==8, "The all-data model distribution includes the former outer error100");
		p=find(prefix+".main.histogram.tsv", "main", "5prime", "signed");
		check(p.terms()==5, "All-data residuals retain the existing signed histogram schema");
		try(Input input=new Input(prefix+".population.tsv")){
			input.reader.nextLine();
			p=new LineParser1('\t').set(input.reader.nextLine());
			check(p.parseLong(0)==8 && p.parseLong(1)==8 && p.parseLong(2)==0
				&& p.termEquals("ALL_DATA_REFIT_NO_HELDOUT_NO_EVALUATION", 5), "Product refit has no held-out population or evaluation claim");
		}
		check(!new File(prefix+".single.tsv").exists() && !new File(prefix+".both.tsv").exists(),
			"All-data refit never manufactures held-out replay results");
		final String[] excluded=java.util.Arrays.copyOf(args, args.length+1);
		excluded[10]="out="+path(root, "refit_excluded");
		excluded[11]="exclude="+audit;
		fails(()->new NcrnaModelEndpointOffsets(excluded).process());
		excluded[10]="out="+path(root, "refit_evaluation");
		excluded[11]="evaluation="+panel;
		fails(()->new NcrnaModelEndpointOffsets(excluded).process());
		excluded[10]="out="+path(root, "refit_lost");
		excluded[8]="expected=8,7";
		excluded[11]="refit=t";
		fails(()->new NcrnaModelEndpointOffsets(excluded).process());
	}

	/** The same medians may yield a forbidden short span; preserve both original errors. */
	private static void minimumLengths(final File root, final String[] legacy){
		final String[] args=java.util.Arrays.copyOf(legacy, legacy.length+1);
		args[10]="thresholds=1";
		args[11]="out="+path(root, "minimum");
		args[13]="minlengths=101,101,1";
		new NcrnaModelEndpointOffsets(args).process();
		LineParser1 p=find(path(root, "minimum.training.tsv"), "main", null, null);
		check(p.terms()==9 && p.parseLong(1)==5 && p.parseLong(2)==5 && p.parseLong(3)==5
			&& p.parseLong(5)==101 && p.parseLong(6)==0 && p.parseLong(7)==5 && p.parseLong(8)==0,
			"Five valid 101-base raw spans center to 100 bases: every pair holds, without losing support or rows");
		p=find(path(root, "minimum.training.tsv"), "second", null, null);
		check(p.parseLong(1)==2 && p.parseLong(3)==2 && p.parseLong(7)==2 && p.parseLong(8)==2,
			"Both reverse-strand raw spans are already short and remain explicit after supported centering");
		p=find(path(root, "minimum.main.histogram.tsv"), "main", "5prime", "signed");
		check(p.parseInt(3)==2 && p.parseLong(4)==3, "Minimum-length holds retain original start errors rather than centered zeros");
		p=find(path(root, "minimum.second.histogram.tsv"), "second", "3prime", "signed");
		check(p.parseInt(3)==-4 && p.parseLong(4)==1, "Reverse-strand holds retain original stop errors with their sign");
		args[11]="out="+path(root, "minimum_equal");
		args[13]="minlengths=100,100,1";
		new NcrnaModelEndpointOffsets(args).process();
		p=find(path(root, "minimum_equal.training.tsv"), "main", null, null);
		check(p.parseLong(3)==0 && p.parseLong(7)==0, "A centered span exactly equal to the family minimum is valid");
		p=find(path(root, "minimum_equal.training.tsv"), "second", null, null);
		check(p.parseLong(3)==0 && p.parseLong(8)==2, "Reverse-strand centering can repair an originally short span without deleting its audit count");
		args[11]="out="+path(root, "minimum_bad_count");
		args[13]="minlengths=100,100";
		fails(()->new NcrnaModelEndpointOffsets(args));
		args[13]="minlengths=100,0,1";
		fails(()->new NcrnaModelEndpointOffsets(args));
		args[13]="minlengths=100,100,1";
		args[12]="centeredhistograms=f";
		fails(()->new NcrnaModelEndpointOffsets(args));
		check(!new File(path(root, "minimum_bad_count.training.tsv")).exists(),
			"Invalid minimum-length bindings must fail before publishing an export");
	}

	/** Literal old/pooled axes distinguish stable support from newly learned medians. */
	private static void compare(final File root){
		final String header="member\tmodel\tend\tn\tmedian_error_bp\tq25_bp\tq75_bp\tiqr_bp\twithin3_of_median\twithin3_of_median_pct\n";
		final String oldText=header
			+"main\tA\t5prime\t5\t2\t2\t2\t0\t5\t100\n"
			+"main\tA\t3prime\t5\t3\t3\t3\t0\t5\t100\n"
			+"second\tB\t5prime\t2\t-3\t-3\t-3\t0\t2\t100\n"
			+"second\tB\t3prime\t2\t-4\t-4\t-4\t0\t2\t100\n"
			+"third\tC\t5prime\t0\tNA\tNA\tNA\tNA\t0\tNA\n"
			+"third\tC\t3prime\t0\tNA\tNA\tNA\tNA\t0\tNA\n";
		final String newText=header
			+"main\tA\t5prime\t10\t4\t4\t4\t0\t10\t100\n"
			+"main\tA\t3prime\t10\t0\t0\t0\t0\t10\t100\n"
			+"second\tB\t5prime\t5\t-2\t-2\t-2\t0\t5\t100\n"
			+"second\tB\t3prime\t5\t-4\t-4\t-4\t0\t5\t100\n"
			+"third\tC\t5prime\t4\t7\t7\t7\t0\t4\t100\n"
			+"third\tC\t3prime\t4\t-7\t-7\t-7\t0\t4\t100\n";
		final String oldPath=write(root, "stability_before.tsv", oldText);
		final String newPath=write(root, "stability_after.tsv", newText);
		final String prefix=path(root, "stability");
		final String[] args={"compare="+oldPath+','+newPath, "members=main,second,third",
			"family=clan", "expectedmodels=3", "thresholds=3", "out="+prefix};
		new NcrnaModelEndpointOffsets(args).process();
		LineParser1 p=find(prefix+".stability.tsv", "main", "A", "3prime");
		check(p.parseLong(3)==5 && p.parseLong(4)==10 && p.parseLong(7)==-3,
			"Saved-table comparison preserves signed new-minus-old median changes");
		p=find(prefix+".stability.tsv", "third", "C", "5prime");
		check(p.termEquals("NA", 5) && p.termEquals("NA", 7) && p.parseInt(10)==0 && p.parseInt(11)==1,
			"A previously absent median is NA, not a fabricated zero or a stability observation");
		p=find(prefix+".stability_summary.tsv", "clan", "5prime", null);
		check(p.terms()==12 && p.parseLong(3)==3 && p.parseLong(4)==1 && p.parseLong(5)==3
			&& p.parseLong(6)==1 && p.parseLong(7)==2 && p.parseLong(8)==0 && p.parseLong(9)==1,
			"One old-supported model changes by two; two newly supported models never enter stable-support bins");
		p=find(prefix+".stability_summary.tsv", "clan", "3prime", null);
		check(p.parseLong(9)==0 && p.parseLong(10)==1 && p.parseLong(11)==0,
			"A signed minus-three median change belongs in the three-to-five absolute-change bin");
		args[5]="out="+path(root, "stability_bad_axis");
		write(newPath, newText.replace("third\tC", "third\tD"));
		fails(()->new NcrnaModelEndpointOffsets(args).process());
		check(!new File(path(root, "stability_bad_axis.stability.tsv")).exists(), "Axis mismatch must precede result publication");
		write(newPath, newText.replace("prime\t10\t", "prime\t4\t"));
		args[5]="out="+path(root, "stability_lost_support");
		fails(()->new NcrnaModelEndpointOffsets(args).process());
		write(newPath, newText);
		write(oldPath, oldText.replace("prime\t0\tNA", "prime\t0\t0"));
		args[5]="out="+path(root, "stability_fake_zero");
		fails(()->new NcrnaModelEndpointOffsets(args).process());
	}

	private static LineParser1 findResult(final String file, final String cohort, final int threshold, final String end){
		try(Input input=new Input(file)){
			input.reader.nextLine();
			final LineParser1 p=new LineParser1('\t');
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				p.set(line);
				if(p.termEquals("new1000", 0) && p.termEquals(cohort, 1) && p.parseInt(2)==threshold && p.termEquals(end, 3)){
					return new LineParser1('\t').set(line.clone());
				}
			}
		}
		throw new AssertionError("Missing finite replay result");
	}
	private static LineParser1 findBreakdown(final String file, final String median, final String end){
		try(Input input=new Input(file)){
			input.reader.nextLine();
			final LineParser1 p=new LineParser1('\t');
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				p.set(line);
				if(p.termEquals("new1000", 3) && p.termEquals(median, 1) && p.termEquals(end, 6)){
					return new LineParser1('\t').set(line.clone());
				}
			}
		}
		throw new AssertionError("Missing literal median/spread breakdown bin");
	}
	private static LineParser1 find(final String file, final String first, final String second, final String third){
		try(Input input=new Input(file)){
			input.reader.nextLine();
			final LineParser1 p=new LineParser1('\t');
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				p.set(line);
				if(p.termEquals(first, 0) && (second==null || p.termEquals(second, 1)) && (third==null || p.termEquals(third, 2))){
					return new LineParser1('\t').set(line.clone());
				}
			}
		}
		throw new AssertionError("Missing finite model/support result");
	}
	private static String path(final File root, final String name){return new File(root, name).getPath();}
	private static String write(final File root, final String name, final String text){
		final String file=path(root, name);
		write(file, text);
		return file;
	}
	private static void write(final String file, final String text){
		// Negative cases deliberately replace this fixture's own input records.
		try(Output output=new Output(new fileIO.ByteStreamWriter(file, true, false, false))){output.writer.print(text);}
	}
	private static void fails(final Runnable action){
		try{action.run();}catch(IllegalArgumentException expected){return;}
		throw new AssertionError("Malformed witness input must fail loudly");
	}
	private static void check(final boolean valid, final String message){if(!valid){throw new AssertionError(message);}}
	private static final String HEADER="input_id\treference_id\temitted_member\tmodel\tseqid\tcall_start\tcall_stop\tref_start\tref_stop\tstrand\td5\td3\n";
}
