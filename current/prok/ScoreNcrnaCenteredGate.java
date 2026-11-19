package prok;

import java.io.File;
import java.util.Locale;

import fileIO.ReadWrite;
import map.ObjectSet;
import map.ObjectIntMap;
import ml.CellNet;
import ml.CellNetParser;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/** Scores each distinct whole-move proposal once and preserves every real call
 * for later paired coordinate replay. No first-stage inference, fitting, cutoff
 * selection or filtering is performed. Scores are stored as raw FP32 bits.
 * @author Brian, Ganyu
 */
public final class ScoreNcrnaCenteredGate{
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, ScoreNcrnaCenteredGate.class, false);
		try{
			new ScoreNcrnaCenteredGate(pp.args).process(null);
			pp.outstream.println("NCRNA_CENTERED_GATE_SCORES_COMPLETE_NO_CUTOFF_NO_FIRST_STAGE_INFERENCE");
		}finally{Shared.closeStream(pp.outstream);}
	}

	public ScoreNcrnaCenteredGate(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			final String a=pair[0].toLowerCase(Locale.ROOT);
			final String b=pair.length>1 ? pair[1] : null;
			if(a.equals("member")){member=b;}
			else if(a.equals("end")){end=b;}
			else if(a.equals("net")){network=b;}
			else if(a.equals("seed")){seed=Integer.parseInt(b);}
			else if(a.equals("partition")){
				require("internal".equals(b) || "outer".equals(b), "partition=internal or outer is required");
				partition=b;
			}else{require(parser.parse(arg, a, b), "Unknown centered gate-score option: "+arg);}
		}
		manifest=parser.in1;
		output=parser.out1;
		require(manifest!=null && output!=null && network!=null && partition!=null && seed>=58 && seed<=60
			&& ("RF00174".equals(member) || "RF01482".equals(member)) && ("5prime".equals(end) || "3prime".equals(end)),
			"Bind real-call manifest, gate network, member, end, seed58/59/60, explicit partition and fresh output directory");
		require(!parser.overwrite && !parser.append && parser.in2==null && parser.out2==null
			&& parser.maxReads<0 && parser.samplerate==1 && !new File(output).exists()
			&& Tools.testInputFiles(false, true, manifest, network), "Score complete populations once into fresh outputs");
	}

	/*--------------------------------------------------------------*/
	/*----------------        Cached Scoring        ----------------*/
	/*--------------------------------------------------------------*/

	/** An injected finite-test network cannot alter input validation or counts. */
	public void process(final CellNet fixture){
		final CellNet net=fixture==null ? CellNetParser.load(network, false) : fixture;
		require(net.numInputs()==35 && net.numOutputs()==1, "The centered whole-move gate requires35inputs and one score");
		require(new File(output).mkdir(), "Cannot create fresh gate-score output directory");
		final boolean bgzf=ReadWrite.USE_BGZF, force=ReadWrite.FORCE_BGZIP;
		ReadWrite.USE_BGZF=true;
		ReadWrite.FORCE_BGZIP=true;
		try(Output cache=new Output(path("scores.tsv.gz")); Output accounting=new Output(path("datasets.tsv")); Input in=new Input(manifest)){
			cache.writer.print(new ByteBuilder(MAGIC).append("\n#schema\t").append(NcrnaCenteredGateFeatures.SCHEMA)
				.append("\n#member\t").append(member).append("\n#end\t").append(end).append("\n#seed\t").append(seed)
				.append("\n#partition\t").append(partition).append("\n#network\t").append(network).nl());
			cache.writer.println(OPTIONAL_HEADER);
			cache.writer.print(cacheHeader().nl());
			accounting.writer.println("dataset\tfold\tmatched\tproposals\tfirst_stage_scored\tgate_seed\tfirst_stage_seed");
			RrnaResourceIO.header(in.reader, NcrnaCenteredGateInputs.MANIFEST_HEADER, manifest);
			final LineParser1 p=new LineParser1('\t');
			final ObjectSet<String> used=new ObjectSet<String>(String.class);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				require(p.terms()==5 && p.length(1)>0 && p.termStartsWith("/", 2), "Gate scoring needs exact dataset/prefix/count bindings");
				final int fold=p.parseInt(0);
				final String dataset=p.parseString(1), prefix=p.parseString(2);
				final long expected=p.parseLong(3), proposals=p.parseLong(4);
				require(partition.equals("outer") ? fold==-1 : fold>=0 && fold<5, "Scoring cannot mix internal and outer populations");
				require(expected>=0 && proposals>=0 && proposals<=expected && used.add(fold+"\t"+dataset), "Dataset counts/identities must be valid and unique");
				scoreDataset(net, prefix, dataset, fold, expected, proposals, cache, accounting);
			}
			require(total>0, "A complete gate-score result needs a nonempty real-call population");
			cache.writer.print(new ByteBuilder("#COMPLETE\t").append(total).tab().append(scored).nl());
		}finally{
			ReadWrite.USE_BGZF=bgzf;
			ReadWrite.FORCE_BGZIP=force;
		}
	}

	private void scoreDataset(final CellNet net, final String prefix, final String dataset, final int fold,
			final long expected, final long expectedProposals, final Output cache, final Output accounting){
		final int firstSeed=fold==-1 ? seed : 58;
		final LineParser1 p=new LineParser1('\t'), v=new LineParser1('\t'), s=new LineParser1('\t');
		final ByteBuilder row=new ByteBuilder();
		final float[] features=new float[35], heads=new float[22];
		long n=0, moves=0, inferred=0;
		int last=0;
		boolean complete=false;
		try(Input audit=new Input(prefix+".rows.tsv.gz"); Input matrix=new Input(prefix+".vectors.tsv.gz"); Input scores=new Input(prefix+".scores.tsv.gz")){
			final String first=readHeader(audit, dataset, fold, firstSeed);
			require(first.equals(readHeader(scores, dataset, fold, firstSeed)), "Audit and saved first-stage outputs must bind the identical network");
			final boolean legacyMaximum=EvaluateNcrnaRealEndpoints.readGeometryHeader(audit.reader, true, prefix);
			row.clear().append("row1\tstatus");
			for(int i=0; i<22; i++){row.append("\tscore_bits_").append(i);}
			RrnaResourceIO.header(scores.reader, row.toString(), prefix);
			RrnaResourceIO.header(matrix.reader, NcrnaCenteredGateInputs.MATRIX_HEADER, prefix);
			for(byte[] line; (line=audit.reader.nextLine())!=null;){
				p.set(line);
				if(p.termEquals("#COMPLETE", 0)){
					require(p.terms()==3 && p.parseLong(1)==n && n==expected && p.parseLong(2)==moves && moves==expectedProposals
						&& audit.reader.nextLine()==null && matrix.reader.nextLine()==null, "Gate-score inputs must close at the accepted call and proposal counts");
					final byte[] tail=scores.reader.nextLine();
					require(tail!=null, "Saved first-stage scores need a complete terminal marker");
					s.set(tail);
					require(s.terms()==3 && s.termEquals("#COMPLETE", 0) && s.parseLong(1)==n && s.parseLong(2)==inferred
						&& scores.reader.nextLine()==null, "Saved first-stage inference accounting must match actual available ends");
					complete=true;
					break;
				}
				require(p.terms()==28 && p.parseInt(0)>last && p.termEquals(member, 4) && p.termEquals(end, 20), "Real gate rows retain ordered original identity, member and endpoint");
				final int maximum=NcrnaTextIO.maximumLength(p, 15, legacyMaximum);
				require(p.parseInt(14)>0 && maximum>=p.parseInt(14), "Gate rows retain the actual family length bounds");
				last=p.parseInt(0);
				final String assembly=p.parseString(2);
				if(assemblyFolds.contains(assembly)){
					require(assemblyFolds.get(assembly)==fold, "One assembly cannot enter different scored folds");
				}else{assemblyFolds.put(assembly, fold);}
				require(seen.add(assembly+'\t'+p.parseString(3)), "Duplicate assembly-qualified references cannot inflate a scored cohort");
				final byte[] scoreLine=scores.reader.nextLine();
				require(scoreLine!=null, "Missing first-stage row");
				s.set(scoreLine);
				final String status=p.parseString(21);
				require((s.terms()==2 || s.terms()==24) && s.parseInt(0)==last && s.termEquals(status, 1), "First-stage row identity/status must match the complete call audit");
				final boolean available=s.terms()==24;
				if(available){
					for(int i=0; i<22; i++){
						heads[i]=Float.intBitsToFloat(s.parseInt(i+2));
						require(Float.isFinite(heads[i]), "Saved first-stage evidence must be finite");
					}
					inferred++;
				}
				final int e=end.equals("5prime") ? 0 : 1;
				final int raw=p.parseInt(10+e), truth=p.parseInt(12+e), center=p.parseInt(16+e);
				final long proposal=p.parseLong(23), gateRow=p.parseLong(22);
				final int neural=available ? NcrnaCenteredGateFeatures.adjustment(heads) : 0;
				require(p.parseInt(24)==neural && p.parseLong(25)==(long)center-raw,
					"Cached neural and median displacements must retain the qualified whole-move semantics");
				float score=0;
				if(gateRow>0){
					require(available && status.equals("PROPOSAL") && proposal==(long)center+neural && proposal!=raw && gateRow==++moves,
						"Only actual distinct-from-raw candidates may invoke the gate");
					final boolean better=NcrnaCenteredGateFeatures.improves(raw, proposal, truth);
					final byte[] vector=matrix.reader.nextLine();
					require(vector!=null, "Proposal matrix ended prematurely");
					v.set(vector);
					require(v.terms()==36 && v.termEquals(better ? "1" : "0", 35) && p.parseInt(26)==(better ? 1 : 0), "Matrix and audit targets must compare with RAW");
					for(int i=0; i<35; i++){
						features[i]=v.parseFloat(i);
						require(Float.isFinite(features[i]), "Gate input cannot contain a fabricated nonfinite value");
						if(i<22){require(Float.floatToRawIntBits(features[i])==Float.floatToRawIntBits(heads[i]), "Gate features must preserve the exact cached first-stage outputs");}
					}
					require(features[32]==Math.abs((long)center-raw) && features[33]==(long)center-raw && features[34]==(neural==0 ? 0 : 1), "Gate displacement inputs disagree with their original real call");
					net.applyInput(features);
					net.feedForward();
					score=net.getOutput(0);
					require(Float.isFinite(score), "A scored whole-move proposal requires finite gate evidence");
				}else{
					require(gateRow==0 && proposal==raw && p.parseInt(26)==0
						&& (available ? status.equals("SAME_AS_RAW_HOLD") && (long)center+neural==raw
						: status.equals("FIXED") || status.equals("INVALID_CENTER_RAW_HOLD") || status.equals("MISSING_FEATURE_HOLD")),
						"Deterministic holds stay raw and cannot silently acquire a gate inference");
				}
				// Normalize old numeric unbounded controls before copying their text into a new cache.
				p.setBounds(15);
				row.clear().append(dataset).tab().append(fold).tab().append(line, 0, p.a());
				NcrnaTextIO.appendMaximumLength(row, maximum);
				row.append(line, p.b(), line.length-p.b()).tab();
				if(gateRow>0){row.append(Float.floatToRawIntBits(score));}
				for(int i=0; i<22; i++){
					row.tab();
					if(available){row.append(Float.floatToRawIntBits(heads[i]));}
				}
				cache.writer.print(row.nl());
				n++;
			}
		}
		require(complete, "Truncated real-gate data cannot produce accepted cached scores");
		total+=n;
		scored+=moves;
		accounting.writer.print(row.clear().append(dataset).tab().append(fold).tab().append(n).tab().append(moves)
			.tab().append(inferred).tab().append(seed).tab().append(firstSeed).nl());
	}

	private String readHeader(final Input in, final String dataset, final int fold, final int firstSeed){
		RrnaResourceIO.header(in.reader, "#NCRNA_CENTERED_REAL_GATE_V1", dataset);
		RrnaResourceIO.header(in.reader, "#schema\t"+NcrnaCenteredGateFeatures.SCHEMA, dataset);
		RrnaResourceIO.header(in.reader, "#panel\t"+dataset, dataset);
		RrnaResourceIO.header(in.reader, "#member\t"+member, dataset);
		RrnaResourceIO.header(in.reader, "#fold\t"+fold, dataset);
		final LineParser1 p=new LineParser1('\t');
		String first=null;
		for(final String key:new String[]{"#geometry", "#folds", "#offsets", "#network"}){
			final byte[] row=in.reader.nextLine();
			require(row!=null, "Truncated real gate provenance");
			p.set(row);
			require(p.terms()==2 && p.termEquals(key, 0) && p.termStartsWith("/", 1), "Real gate provenance requires each absolute upstream resource");
			if(key.equals("#network")){first=p.parseString(1);}
		}
		RrnaResourceIO.header(in.reader, "#seed\t"+firstSeed, dataset);
		RrnaResourceIO.header(in.reader, "#end\t"+end, dataset);
		return first;
	}

	public static ByteBuilder cacheHeader(){
		return cacheHeader(false);
	}
	/** Old versions share the numeric-maximum header; only new caches permit its absence. */
	public static ByteBuilder cacheHeader(final boolean legacyMaximum){
		final ByteBuilder row=new ByteBuilder("dataset\tfold\t")
			.append(legacyMaximum ? EvaluateNcrnaRealEndpoints.LEGACY_GATE_HEADER : EvaluateNcrnaRealEndpoints.GATE_HEADER).append("\tgate_score_bits");
		for(int i=0; i<22; i++){row.append("\tfirst_score_bits_").append(i);}
		return row;
	}
	private String path(final String file){return new File(output, file).getPath();}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}
	/*--------------------------------------------------------------*/
	/*----------------       Fields / Constants     ----------------*/
	/*--------------------------------------------------------------*/
	private String manifest, output, member, end, network, partition;
	private int seed=-1;
	private long total, scored;
	private final ObjectSet<String> seen=new ObjectSet<String>(String.class);
	private final ObjectIntMap<String> assemblyFolds=new ObjectIntMap<String>(String.class);
	public static final String MAGIC="#NCRNA_CENTERED_GATE_SCORES_V3";
	public static final String NUMERIC_MAXIMUM_MAGIC="#NCRNA_CENTERED_GATE_SCORES_V2";
	public static final String LEGACY_MAGIC="#NCRNA_CENTERED_GATE_SCORES_V1";
	/** Scores are absent when no inference occurred; status and proposal row determine presence. */
	public static final String OPTIONAL_HEADER="#Optional scores: gate_score_bits is empty unless status=PROPOSAL; first_score_bits_* are all empty unless first-stage inference occurred.";
}
