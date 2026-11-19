package test.prok;

import java.io.File;
import java.util.Locale;

import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.EvaluateNcrnaCenteredEndpoints;
import prok.NcrnaFixedEndpointTrials;
import prok.NcrnaLearningPopulation;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Literal original-frame error arithmetic, holds and invalid-proposal controls.
 * No model, training data or caller is needed.
 * @author Brian, Ganyu
 */
public final class EvaluateNcrnaCenteredEndpointsTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, EvaluateNcrnaCenteredEndpointsTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				check(parser.parse(arg, pair[0].toLowerCase(Locale.ROOT), pair.length>1 ? pair[1] : null), "Unknown finite fixture option");
			}
			check(parser.out1!=null && parser.in1==null && parser.in2==null && parser.out2==null
				&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1, "Use out=fresh_directory and native runtime flags");
			final File root=new File(parser.out1);
			check(root.mkdir(), "Finite fixture directory must be fresh");
			run(root, false, false);
			run(root, true, false);
			run(root, false, true);
			run(root, true, true);
			pp.outstream.println("NCRNA_CENTERED_ENDPOINT_COMPARISON_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Every expected metric below is derived from literal distances, not the tested scorer. */
	private static void run(final File root, final boolean holdsOnly, final boolean stop){
		final String prefix=new File(root, (holdsOnly ? "holds" : "mixed")+(stop ? "3" : "5")).getPath();
		final ByteBuilder nativeRows=new ByteBuilder(NcrnaLearningPopulation.NATIVE_HEADER).nl();
		final ByteBuilder ledger=new ByteBuilder(NcrnaFixedEndpointTrials.AUDIT_HEADER).nl();
		final ByteBuilder baseline=scoreHeader(stop), centered=scoreHeader(stop);
		final int[] rawStart={25,25,25,25,25,1,25}, rawStop={44,44,44,44,34,20,44};
		final int[] m5={5,3,5,-50,0,0,30}, m3={5,3,5,0,0,0,0};
		final int[] oldMove={-3,-2,-5,-5,-10,-10,0}, newMove={0,-2,0,0,0,10,0};
		final String[] status={"VECTORIZED", "VECTORIZED", "CENTER_MISSING_FEATURE_HOLD", "CENTER_CROSSED_HOLD",
			"CENTER_BELOW_MINIMUM_HOLD", "VECTORIZED", "CENTER_OUTSIDE_CONTEXT_HOLD"};
		if(stop){
			// Reflect around the middle of the 60-base context: all error magnitudes
			// remain literal, but the scorer must use stop coordinates and opposite signs.
			for(int i=0; i<7; i++){
				final int start=rawStart[i], median=m5[i];
				rawStart[i]=59-rawStop[i];
				rawStop[i]=59-start;
				m5[i]=-m3[i];
				m3[i]=-median;
				oldMove[i]=-oldMove[i];
				newMove[i]=-newMove[i];
			}
		}
		int n=0, usable=0;
		for(int i=0; i<7; i++){
			if(holdsOnly && i!=3 && i!=4 && i!=6){continue;}
			n++;
			nativeRows.append('g').append(i).tab().append('g').append(i)
				.append("\tsource\tassembly\t1\tBacteria\tphylum\tcontig\t+\t100\t119\t20\t20\t20\n");
			final ByteBuilder identity=identity(i, rawStart[i], rawStop[i]);
			ledger.append(identity).tab().append(n).tab().append(rawStart[i]-m5[i]).tab().append(rawStop[i]-m3[i])
				.tab().append(m5[i]).tab().append(m3[i]).tab().append(status[i]).tab().append(i==4 ? 0 : 10)
				.append("\t10\t").append(i==4 ? 0 : 1).append("\t16\n");
			score(baseline, identity, stop ? 39-rawStop[i] : 20-rawStart[i], oldMove[i]);
			if(status[i].equals("VECTORIZED")){
				usable++;
				score(centered, identity(i, rawStart[i]-m5[i], rawStop[i]-m3[i]),
					stop ? 39-rawStop[i]+m3[i] : 20-rawStart[i]+m5[i], newMove[i]);
			}
		}
		baseline.append("#COMPLETE\t").append(n).nl();
		centered.append("#COMPLETE\t").append(usable).nl();
		write(prefix+".loci.tsv", nativeRows);
		write(prefix+".ledger.tsv", ledger);
		write(prefix+".old.tsv", baseline);
		if(!holdsOnly){write(prefix+".new.tsv", centered);}
		final String[] args=holdsOnly ? new String[7] : new String[8];
		args[0]="in="+prefix+".old.tsv";
		args[1]="ledger="+prefix+".ledger.tsv";
		args[2]="loci="+prefix+".loci.tsv";
		args[3]=stop ? "end=3prime" : "end=5prime";
		args[4]="expected="+n;
		args[5]="usable="+usable;
		args[6]="out="+prefix+".result";
		if(!holdsOnly){args[7]="centered="+prefix+".new.tsv";}
		new EvaluateNcrnaCenteredEndpoints(args).process();
		LineParser1 p=find(prefix+".result.summary.tsv", "centered_margin0");
		check(p.terms()==33 && p.parseLong(2)==n && p.parseLong(24)==1 && p.parseLong(28)==1,
			"Every original trial survives, with already-short raw spans and unsupported models explicit");
		geometry(prefix, ledger, n, stop, holdsOnly);
		if(holdsOnly){
			check(p.parseLong(3)==0 && p.parseDouble(8)==5 && p.parseLong(25)==0 && p.parseLong(27)==3
				&& p.parseLong(29)==0 && p.parseLong(30)==0 && p.parseLong(31)==3,
				"All-held assays require no centered score stream and preserve all raw endpoints");
			return;
		}
		check(p.parseLong(3)==3 && p.parseLong(4)==3 && p.parseLong(5)==7 && close(p.parseDouble(8), Math.sqrt(156/7.0))
			&& close(p.parseDouble(9), Math.sqrt(511/7.0)) && close(p.parseDouble(12), Math.sqrt(904/7.0))
			&& p.parseDouble(13)<0 && p.parseLong(20)==2 && p.parseLong(23)==2
			&& p.parseLong(25)==3 && p.parseLong(26)==1 && p.parseLong(27)==3
			&& p.parseLong(29)==2 && p.parseLong(30)==2 && p.parseLong(31)==3,
			"Centered scoring uses original-frame errors, keeps a final too-short network proposal, and partitions center/network/raw paths");
		p=find(prefix+".result.summary.tsv", "center_only");
		check(close(p.parseDouble(8), Math.sqrt(440/7.0)) && p.parseLong(3)==2 && p.parseLong(29)==4,
			"Valid missing-feature trials retain the center, while invalid centers return to the original pair");
		p=find(prefix+".result.summary.tsv", "uncentered_margin0");
		check(p.parseLong(20)==1 && p.parseLong(22)==1 && p.parseLong(23)==0,
			"Historical uncentered proposals are counted outside context without adding a new hold");
		final String complete=centered.toString();
		write(prefix+".new.tsv", new ByteBuilder(complete.substring(0, complete.lastIndexOf("#COMPLETE"))));
		args[6]="out="+prefix+".truncated";
		fails(()->new EvaluateNcrnaCenteredEndpoints(args).process());
		check(!new File(prefix+".truncated.summary.tsv").exists(), "Incomplete score streams never publish a summary");
		write(prefix+".new.tsv", new ByteBuilder(complete.replace("g1\tsource", "WRONG\tsource")));
		args[6]="out="+prefix+".wrong_identity";
		fails(()->new EvaluateNcrnaCenteredEndpoints(args).process());
		write(prefix+".new.tsv", new ByteBuilder(complete));
		write(prefix+".ledger.tsv", new ByteBuilder(ledger.toString().replace("CENTER_MISSING_FEATURE_HOLD", "CENTER_CROSSED_HOLD")));
		args[6]="out="+prefix+".wrong_hold";
		fails(()->new EvaluateNcrnaCenteredEndpoints(args).process());
	}

	/** Equal raw-invalid and center-invalid marginals can hide a repaired raw
	 * span and a different newly short center; verify their actual intersection.
	 */
	private static void geometry(final String prefix, final ByteBuilder ledger, final int n, final boolean stop, final boolean holdsOnly){
		final String[] args={"geometryonly=t", "ledger="+prefix+".ledger.tsv", "loci="+prefix+".loci.tsv",
			stop ? "end=3prime" : "end=5prime", "expected="+n, "out="+prefix+".audit"};
		new EvaluateNcrnaCenteredEndpoints(args).process();
		String path=prefix+".audit.geometry.tsv";
		check(geometryCount(path, "BELOW_MINIMUM", "BELOW_MINIMUM", n)==1
			&& geometryCount(path, "VALID", "CROSSED", n)==1 && geometryCount(path, "VALID", "OUTSIDE_CONTEXT", n)==1
			&& geometryCount(path, "BELOW_MINIMUM", "VALID", n)==0 && geometryCount(path, "VALID", "BELOW_MINIMUM", n)==0,
			"The original geometry table has one inherited invalid span, two new invalid centers and no repair");
		final String original=ledger.toString();
		String changed=original.replace("\t25\t34\t0\t0\tCENTER_BELOW_MINIMUM_HOLD\t0\t10\t0\t16",
			"\t15\t34\t10\t0\tVECTORIZED\t10\t10\t1\t16");
		if(!holdsOnly){
			changed=changed.replace(stop ? "\t20\t39\t-5\t-5\tCENTER_MISSING_FEATURE_HOLD" : "\t20\t39\t5\t5\tCENTER_MISSING_FEATURE_HOLD",
				stop ? "\t20\t34\t-5\t0\tCENTER_BELOW_MINIMUM_HOLD" : "\t20\t34\t5\t10\tCENTER_BELOW_MINIMUM_HOLD");
		}
		check(!changed.equals(original), "Literal control must actually repair the initially short span");
		write(prefix+".ledger.tsv", new ByteBuilder(changed));
		args[5]="out="+prefix+".repaired_audit";
		new EvaluateNcrnaCenteredEndpoints(args).process();
		path=prefix+".repaired_audit.geometry.tsv";
		check(geometryCount(path, "BELOW_MINIMUM", "BELOW_MINIMUM", n)==0
			&& geometryCount(path, "BELOW_MINIMUM", "VALID", n)==1
			&& geometryCount(path, "VALID", "BELOW_MINIMUM", n)==(holdsOnly ? 0 : 1),
			"Joint geometry counts identify a repaired raw span and a different newly short center, without network scores");
		write(prefix+".ledger.tsv", new ByteBuilder(original));
		check(!new File(prefix+".audit.summary.tsv").exists() && !new File(prefix+".audit.trials.tsv").exists(),
			"Geometry-only mode must never regenerate endpoint metrics or predictions");
	}
	private static long geometryCount(final String path, final String raw, final String center, final int expected){
		long found=-1, total=0;
		int rows=0;
		try(Input input=new Input(path)){
			input.reader.nextLine();
			final LineParser1 p=new LineParser1('\t');
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				p.set(line);
				check(p.terms()==5 && p.parseLong(3)>=0, "Geometry cross-tab uses five fields and nonnegative counts");
				rows++;
				total+=p.parseLong(3);
				if(p.termEquals(raw, 1) && p.termEquals(center, 2)){
					check(found<0, "Each geometry pair occurs exactly once");
					found=p.parseLong(3);
				}
			}
		}
		check(rows==16 && total==expected && found>=0, "The full geometry cross-tab preserves every original trial including zero cells");
		return found;
	}

	private static ByteBuilder identity(final int i, final int start, final int stop){
		return new ByteBuilder().append('g').append(i).append("\tsource\tassembly\tevaluation\t").append(i)
			.append("\tmodel").append(i).append("\t0\t").append(start).tab().append(stop).append("\t20\t39");
	}
	private static ByteBuilder scoreHeader(final boolean stop){
		final ByteBuilder b=new ByteBuilder("#NCRNA_FULL_ERROR_SCORES_V1\n#coding\tFP32_RAW_SIGNED_INT\n#matrix\tfixture\n#rows\tfixture\n#network\tfixture\n#end\t")
			.append(stop ? "3prime" : "5prime").nl()
			.append(NcrnaFixedEndpointTrials.ROW_HEADER).append("\ttrue_class\ttrue_minus_raw_bp\tmargin0_class\tmargin01_class\targmax22_class");
		for(int i=0; i<22; i++){b.append("\tscore_bits_").append(i);}
		return b.nl();
	}
	private static void score(final ByteBuilder b, final ByteBuilder identity, final int correction, final int shift){
		final int target=correction<-10 || correction>10 ? 21 : correction+10;
		final int selected=shift+10;
		b.append(identity).tab().append(target).tab().append(correction).tab().append(selected).tab().append(selected).tab().append(selected);
		for(int i=0; i<22; i++){b.tab().append(Float.floatToRawIntBits(i==selected ? 1f : 0f));}
		b.nl();
	}
	private static LineParser1 find(final String path, final String policy){
		try(Input input=new Input(path)){
			input.reader.nextLine();
			final LineParser1 p=new LineParser1('\t');
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				p.set(line);
				if(p.termEquals(policy, 0)){return new LineParser1('\t').set(line.clone());}
			}
		}
		throw new AssertionError("Missing finite policy row: "+policy);
	}
	private static void write(final String path, final ByteBuilder text){
		try(Output output=new Output(new fileIO.ByteStreamWriter(path, true, false, false))){output.writer.print(text);}
	}
	private static boolean close(final double a, final double b){return Math.abs(a-b)<1e-7;}
	private static void fails(final Runnable action){
		try{action.run();}catch(IllegalArgumentException expected){return;}
		throw new AssertionError("Malformed saved input must fail loudly");
	}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
}
