package prok;

import java.io.File;
import java.util.Locale;

import fileIO.ReadWrite;
import map.ObjectIntMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/** Combines already-qualified real held-call gate rows without learning a
 * transform. Folds0..3 teach; fold4 monitors. RegressionTrainer fits scaling
 * from teaching rows only. Outer rows can never enter this input route.
 * All physical-call audit rows, including holds, are checked before acceptance.
 * @author Brian, Ganyu
 */
public final class NcrnaCenteredGateInputs{
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaCenteredGateInputs.class, false);
		try{
			new NcrnaCenteredGateInputs(pp.args).process();
			pp.outstream.println("NCRNA_CENTERED_GATE_INPUTS_COMPLETE_REAL_FOLDS_NO_FIT");
		}finally{Shared.closeStream(pp.outstream);}
	}

	public NcrnaCenteredGateInputs(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			final String a=pair[0].toLowerCase(Locale.ROOT);
			final String b=pair.length>1 ? pair[1] : null;
			if(a.equals("member")){member=b;}
			else if(a.equals("end")){end=b;}
			else{require(parser.parse(arg, a, b), "Unknown centered-gate input argument: "+arg);}
		}
		manifest=parser.in1;
		output=parser.out1;
		require(manifest!=null && output!=null && ("RF00174".equals(member) || "RF01482".equals(member))
			&& ("5prime".equals(end) || "3prime".equals(end)), "Bind a real-fold manifest, fresh output, cobalamin member and explicit end");
		require(!parser.overwrite && !parser.append && parser.in2==null && parser.out2==null
			&& parser.maxReads<0 && parser.samplerate==1 && !new File(output).exists()
			&& Tools.testInputFiles(false, true, manifest), "Complete real gate populations cannot be sampled or overwritten");
	}

	/*--------------------------------------------------------------*/
	/*----------------         Fold Assembly        ----------------*/
	/*--------------------------------------------------------------*/

	public void process(){
		require(new File(output).mkdir(), "Cannot create fresh centered-gate input directory");
		final boolean oldBgzf=ReadWrite.USE_BGZF;
		final boolean oldForce=ReadWrite.FORCE_BGZIP;
		ReadWrite.USE_BGZF=true;
		ReadWrite.FORCE_BGZIP=true;
		try(Output training=new Output(path("train.tsv.gz")); Output monitor=new Output(path("monitor.tsv.gz"));
				Output accounting=new Output(path("datasets.tsv")); Input input=new Input(manifest)){
			training.writer.println(MATRIX_HEADER);
			monitor.writer.println(MATRIX_HEADER);
			accounting.writer.println("member\tend\tfold\tdataset\tmatched\tproposals\tpositives\tdeterministic_holds\trole\tqualification");
			RrnaResourceIO.header(input.reader, MANIFEST_HEADER, manifest);
			final LineParser1 p=new LineParser1('\t');
			final ObjectSet<String> datasets=new ObjectSet<String>(String.class), prefixes=new ObjectSet<String>(String.class);
			for(byte[] row; (row=input.reader.nextLine())!=null;){
				p.set(row);
				require(p.terms()==5 && p.length(1)>0 && p.termStartsWith("/", 2), "Every dataset needs explicit fold, name, absolute prefix and accepted counts");
				final int fold=p.parseInt(0);
				final String dataset=p.parseString(1);
				final String prefix=p.parseString(2);
				final long calls=p.parseLong(3), proposals=p.parseLong(4);
				require(fold>=0 && fold<5 && calls>=0 && proposals>=0 && proposals<=calls,
					"Only internal held folds0..4 enter gate fitting; outer data are forbidden even when empty");
				require(datasets.add(fold+"\t"+dataset) && prefixes.add(prefix), "A physical dataset/fold and its matrix prefix may enter once only");
				readDataset(prefix, dataset, fold, calls, proposals, fold==4 ? monitor : training, accounting);
				files[fold]++;
			}
			for(int fold=0; fold<5; fold++){
				require(files[fold]>0 && matched[fold]>0 && proposed[fold]>0,
					"Every internal fold must contribute real matched calls and actual gate rows: fold="+fold);
			}
		}finally{
			ReadWrite.USE_BGZF=oldBgzf;
			ReadWrite.FORCE_BGZIP=oldForce;
		}
		try(Output summary=new Output(path("summary.tsv"))){
			summary.writer.println("member\tend\tfold\tdatasets\tmatched\tproposals\tpositives\tholds\tinputs\tschema\trole");
			final ByteBuilder row=new ByteBuilder();
			for(int fold=0; fold<5; fold++){
				summary.writer.print(row.clear().append(member).tab().append(end).tab().append(fold).tab().append(files[fold])
					.tab().append(matched[fold]).tab().append(proposed[fold]).tab().append(positives[fold])
					.tab().append(matched[fold]-proposed[fold]).tab().append(35).tab().append(NcrnaCenteredGateFeatures.SCHEMA)
					.tab().append(fold==4 ? "MONITOR_NOT_NESTED_INDEPENDENT" : "TRAINING_ONLY_SCALING_INPUT").nl());
			}
		}
	}

	/** Streaming merge preserves original matrix bytes; no rescale, redraw,
	 * filtering or clipping. Audit and matrix closure, labels, displacements and
	 * assembly-qualified uniqueness are independent checks on the accepted producer. */
	private void readDataset(final String prefix, final String dataset, final int fold,
			final long expected, final long expectedProposals, final Output output, final Output accounting){
		assert(fold>=0 && fold<5) : "Outer assemblies must not enter gate training or monitoring";
		final LineParser1 p=new LineParser1('\t'), v=new LineParser1('\t');
		final ByteBuilder text=new ByteBuilder();
		long n=0, moves=0, better=0;
		int lastRow=0;
		boolean complete=false;
		try(Input audit=new Input(prefix+".rows.tsv.gz"); Input matrix=new Input(prefix+".vectors.tsv.gz")){
			readHeader(audit, dataset, member, fold, end);
			final boolean legacyMaximum=EvaluateNcrnaRealEndpoints.readGeometryHeader(audit.reader, true, prefix);
			RrnaResourceIO.header(matrix.reader, MATRIX_HEADER, prefix);
			for(byte[] row; (row=audit.reader.nextLine())!=null;){
				p.set(row);
				if(p.termEquals("#COMPLETE", 0)){
					require(p.terms()==3 && p.parseLong(1)==expected && p.parseLong(2)==expectedProposals
						&& n==expected && moves==expectedProposals && audit.reader.nextLine()==null && matrix.reader.nextLine()==null,
						"Gate audit/matrix must close at their accepted full-call and proposal counts: "+prefix);
					complete=true;
					break;
				}
				require(p.terms()==28 && p.parseInt(0)>lastRow && p.termEquals(member, 4) && p.termEquals(end, 20),
					"Real gate audit requires increasing original cache rows, exact member and endpoint");
				require(NcrnaTextIO.maximumLength(p, 15, legacyMaximum)>=p.parseInt(14),
					"Gate audit must preserve a positive minimum and its optional upper bound");
				lastRow=p.parseInt(0);
				final String assembly=p.parseString(2);
				if(assemblyFolds.contains(assembly)){
					require(assemblyFolds.get(assembly)==fold, "One assembly cannot enter multiple held folds, even with distinct reference IDs");
				}else{assemblyFolds.put(assembly, fold);}
				require(seen.add(assembly+'\t'+p.parseString(3)),
					"An assembly-qualified real reference may enter only one held fold; duplicate calls cannot inflate gate training");
				final int e=end.equals("5prime") ? 0 : 1;
				final int raw=p.parseInt(10+e), truth=p.parseInt(12+e), center=p.parseInt(16+e);
				final long proposal=p.parseLong(23);
				final int neural=p.parseInt(24);
				final long displacement=p.parseLong(25);
				final long gateRow=p.parseLong(22);
				require(raw>=0 && truth>=0 && neural>=-10 && neural<=10 && displacement==(long)center-raw,
					"Raw/truth coordinates and signed applied median must retain their real-frame meaning");
				if(gateRow>0){
					require(p.termEquals("PROPOSAL", 21) && proposal!=raw && proposal==(long)center+neural && gateRow==++moves,
						"Gate matrices contain exactly the dense list of distinct whole-move proposals");
					final boolean improves=NcrnaCenteredGateFeatures.improves(raw, proposal, truth);
					require(p.parseInt(26)==(improves ? 1 : 0), "Audit target must mean strictly closer to RAW, including equal-distance negatives");
					final byte[] vector=matrix.reader.nextLine();
					require(vector!=null, "Gate matrix ended before the complete proposal audit");
					v.set(vector);
					require(v.terms()==36 && v.termEquals(improves ? "1" : "0", 35), "Gate matrix shape/target must match the real-call audit");
					for(int i=0; i<35; i++){require(Float.isFinite(v.parseFloat(i)), "Gate training cannot silently replace nonfinite features");}
					require(v.parseFloat(32)==Math.abs(displacement) && v.parseFloat(33)==displacement
						&& v.parseFloat(34)==(neural==0 ? 0 : 1), "Centered gate appendages must match actual median displacement and neural movement");
					output.writer.print(text.clear().append(vector).nl());
					if(improves){better++;}
				}else{
					require(gateRow==0 && proposal==raw && p.parseInt(26)==0
						&& (p.termEquals("FIXED", 21) || p.termEquals("INVALID_CENTER_RAW_HOLD", 21)
						|| p.termEquals("MISSING_FEATURE_HOLD", 21) || p.termEquals("SAME_AS_RAW_HOLD", 21)),
						"Deterministic holds stay in accounting and never become fitted examples");
				}
				n++;
			}
		}
		require(complete, "Truncated real gate audit has no accepted completion: "+prefix);
		matched[fold]+=n;
		proposed[fold]+=moves;
		positives[fold]+=better;
		accounting.writer.print(new ByteBuilder().append(member).tab().append(end).tab().append(fold).tab().append(dataset)
			.tab().append(n).tab().append(moves).tab().append(better).tab().append(n-moves)
			.tab().append(fold==4 ? "MONITOR" : "TRAIN").append("\tREAL_HELD_CALL_ROWS_NO_OUTER_NO_TRANSFORM\n"));
	}

	/** Shared framing for the real-call adapter's audit and score streams. Returns
	 * the exact first-stage network path for consumers that need score provenance. */
	public static String readHeader(final Input input, final String dataset, final String member, final int fold, final String end){
		RrnaResourceIO.header(input.reader, "#NCRNA_CENTERED_REAL_GATE_V1", dataset);
		RrnaResourceIO.header(input.reader, "#schema\t"+NcrnaCenteredGateFeatures.SCHEMA, dataset);
		RrnaResourceIO.header(input.reader, "#panel\t"+dataset, dataset);
		RrnaResourceIO.header(input.reader, "#member\t"+member, dataset);
		RrnaResourceIO.header(input.reader, "#fold\t"+fold, dataset);
		final LineParser1 p=new LineParser1('\t');
		String network=null;
		for(final String key:new String[]{"#geometry", "#folds", "#offsets", "#network"}){
			final byte[] row=input.reader.nextLine();
			require(row!=null, "Truncated real-gate provenance header");
			p.set(row);
			require(p.terms()==2 && p.termEquals(key, 0) && p.termStartsWith("/", 1), "Real gate provenance names each absolute upstream resource");
			if(key.equals("#network")){network=p.parseString(1);}
		}
		final byte[] seed=input.reader.nextLine();
		require(seed!=null, "Real gate provenance requires a seed label");
		p.set(seed);
		require(p.terms()==2 && p.termEquals("#seed", 0) && p.termEquals("58", 1),
			"Only the declared auxiliary seed58 supplies out-of-fold gate teaching rows");
		RrnaResourceIO.header(input.reader, "#end\t"+end, dataset);
		return network;
	}

	private String path(final String file){return new File(output, file).getPath();}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------       Fields / Constants     ----------------*/
	/*--------------------------------------------------------------*/

	private String manifest, output, member, end;
	private final long[] matched=new long[5], proposed=new long[5], positives=new long[5];
	private final int[] files=new int[5];
	private final ObjectSet<String> seen=new ObjectSet<String>(String.class);
	private final ObjectIntMap<String> assemblyFolds=new ObjectIntMap<String>(String.class);
	public static final String MANIFEST_HEADER="fold\tdataset\tprefix\texpected_calls\texpected_proposals";
	public static final String MATRIX_HEADER="#dims\t35\t1\t0";
}
