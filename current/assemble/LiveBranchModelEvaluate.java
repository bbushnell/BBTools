package assemble;

import java.io.File;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.Locale;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import ml.CellNet;
import ml.CellNetParser;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;
import structures.LongList;

/**
 * Scores one live-branch corpus partition without resampling or mixing controls
 * into validation candidates. Emits all attainable {@code score>=cutoff} decisions.
 * Validates every sidecar row while consuming features only for the selected partition.
 * Validation permits grey candidates only; validation_controls permits controls only;
 * training may contain either role. Inference reuses one model and feature vector serially.
 * Retains packed scores by role and delegates in-place sorting/tied-score curves to
 * SameKModelEvaluate. Predictions preserve selected-row order; scores.tsv retains
 * the seven provenance columns and adds score_bits and score. Separate grey.tsv and
 * control.tsv curves share summary.tsv, including a row for an empty role.
 * Partial outputs remain after failure; complete.txt is written only after the other
 * writers finish successfully.
 * @author Fischl
 */
public final class LiveBranchModelEvaluate{

	/**
	 * Runs serial inference through standard BBTools parsing and closes any owned diagnostic stream.
	 * Standard output/error remain open; shared parser settings are not restored after the run.
	 * @param args Feature input, full sidecar, model, partition and new output directory options
	 */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, LiveBranchModelEvaluate.class, false);
		try{new LiveBranchModelEvaluate(pp.args, pp.outstream).run();}
		finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Parses options and requires readable feature, sidecar and model paths and a new output path.
	 * Partition names are case-sensitive: train, validation or validation_controls.
	 * Header/provenance validation and output-directory creation are deferred to run.
	 * @param args Expanded options, including standard parser settings
	 * @param stream Borrowed diagnostic stream; main owns its closing policy
	 */
	private LiveBranchModelEvaluate(final String[] args, final PrintStream stream){
		outstream=stream;
		final Parser parser=new Parser();
		for(String arg : args){
			final String[] parts=arg.split("=", 2);
			String name=parts[0].toLowerCase(Locale.ROOT);
			while(name.startsWith("-")){name=name.substring(1);}
			final String value=parts.length>1 ? parts[1] : null;
			if(name.equals("rows")){
				rows=value;
			}else if(name.equals("net")){
				network=value;
			}else if(name.equals("partition")){
				partition=Arrays.asList(PARTITIONS).indexOf(value);
			}else if(name.equals("simdff")){
				simdFF=Parse.parseBoolean(value);
			}else if(!parser.parse(arg, name, value)){throw new IllegalArgumentException("Unknown evaluation option: "+arg);}
		}
		input=parser.in1; output=parser.out1;
		if(input==null || rows==null || network==null || output==null || partition<0 || parser.in2!=null || parser.out2!=null ||
				!Tools.testInputFiles(false, true, input, rows, network) || new File(output).exists()){
			throw new IllegalArgumentException("Required in=partition.tsv rows=rows.tsv net=model.bbnet partition=validation out=NEW_directory.");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Checks exact headers, one-based partition-local ordinals, identities, roles and binary targets.
	 * Sidecar partitions may be interleaved; only selected rows consume feature rows, in order.
	 * Selected rows additionally require matching feature targets, finite nonnegative values
	 * and a K encoding within 1e-6 of the provenance K. Model scores must be finite in [0,1].
	 * Packed score bits normalize signed zero; printed score text retains the model's value.
	 * The selected partition must be nonempty. Emits separate grey/control curves with tied
	 * groups kept together and a reject-all row, even for a role with no selected rows.
	 * Retains one packed long per selected row until both role curves have been written.
	 * Failed prediction-output acquisition attempts both input closes and rethrows the
	 * original failure with cleanup failures suppressed. Partial output ownership is separate.
	 */
	private void run(){
		final CellNet net=CellNetParser.load(network, false);
		if(net==null || net.numInputs()!=12 || net.numOutputs()!=1){throw new IllegalArgumentException("Expected a12-input, one-output model.");}
		net.simdFF=simdFF;
		if(!new File(output).mkdir()){throw new IllegalArgumentException("Cannot create new evaluation directory.");}
		final LongList[] scores={new LongList(), new LongList()};
		final long[] ordinals=new long[3];
		final float[] vector=new float[12];
		final LineParser1 mp=new LineParser1('\t'), dp=new LineParser1('\t');
		final ByteBuilder buffer=new ByteBuilder(), word=new ByteBuilder(), genome=new ByteBuilder();
		final ByteFile meta=ByteFile.makeByteFile(rows, true);
		final ByteFile data;
		try{data=ByteFile.makeByteFile(input, true);}
		catch(RuntimeException|Error failure){
			try{
				if(meta.close()){failure.addSuppressed(new IllegalStateException("Evaluation metadata cleanup failed."));}
			}catch(RuntimeException|Error cleanup){
				if(cleanup!=failure){failure.addSuppressed(cleanup);}
			}
			throw failure;
		}
		final ByteStreamWriter predictions;
		try{predictions=open(output+"/scores.tsv");}
		catch(RuntimeException|Error failure){
			// ASM024: output acquisition can fail while both inputs are already owned.
			try{
				try{
					if(meta.close()){failure.addSuppressed(new IllegalStateException("Evaluation metadata cleanup failed."));}
				}catch(RuntimeException|Error cleanup){
					if(cleanup!=failure){failure.addSuppressed(cleanup);}
				}
			}finally{
				try{
					try{
						if(data.close()){failure.addSuppressed(new IllegalStateException("Evaluation feature input cleanup failed."));}
					}catch(RuntimeException|Error cleanup){
						if(cleanup!=failure){failure.addSuppressed(cleanup);}
					}
				}finally{throw failure;}// Preserve the primary even if recording a cleanup failure fails.
			}
		}
		try{
			require(meta.nextLine(), ROW_HEADER);
			require(data.nextLine(), "#dims\t12\t1");
			require(data.nextLine(), LiveBranchTrainingCorpus.SCHEMA);
			predictions.print(ROW_HEADER+"\tscore_bits\tscore\n");
			for(byte[] line=meta.nextLine(); line!=null; line=meta.nextLine()){
				mp.set(line);
				if(mp.terms()!=7){throw new IllegalArgumentException("Wrong live-branch provenance width.");}
				int rowPartition=-1;
				for(int i=0; i<3; i++){if(mp.termEquals(PARTITIONS[i], 0)){rowPartition=i;}}
				final int role=mp.termEquals("grey", 5) ? 0 : mp.termEquals("control", 5) ? 1 : -1;
				final int k=mp.parseInt(3), label=mp.parseInt(6);
				mp.appendTerm(genome.clear(), 2); mp.appendTerm(word.clear(), 4);
				if(rowPartition<0 || role<0 || mp.parseLong(1)!=++ordinals[rowPartition] ||
						(label!=0 && label!=1) || !validIdentity(genome, word, k) ||
						(rowPartition==1 && role!=0) || (rowPartition==2 && role!=1)){
					throw new IllegalArgumentException("Invalid partition/ordinal/genome/K/word/role/target in sidecar.");
				}
				if(rowPartition!=partition){continue;}
				final byte[] values=data.nextLine();
				if(values==null){throw new IllegalArgumentException("Missing paired feature row.");}
				dp.set(values);
				if(dp.terms()!=13 || dp.parseInt(12)!=label){throw new IllegalArgumentException("Feature width/target disagrees with provenance.");}
				for(int i=0; i<12; i++){
					vector[i]=dp.parseFloat(i);
					if(!Float.isFinite(vector[i]) || vector[i]<0){throw new IllegalArgumentException("Invalid normalized feature.");}
				}
				final float encodedK=(float)(0.125*Math.log(k)*(1/Math.log(2)));
				if(Math.abs(vector[9]-encodedK)>1e-6){throw new IllegalArgumentException("Feature K disagrees with provenance.");}
				final float score=net.applyInput(vector).feedForward();
				final long packed=SameKModelEvaluate.pack(score, label);
				scores[role].add(packed);
				predictions.print(buffer.clear().append(line).tab().append(packed>>>1).tab().appendSlow(score).nl());
			}
			if(data.nextLine()!=null || ordinals[partition]==0){throw new IllegalArgumentException("Unpaired extra rows or empty selected partition.");}
		}finally{
			//TODO: Probable bug ASM024 - unchecked close failures can skip later cleanup and mask the processing failure.
			if(meta.close() | data.close() | predictions.poisonAndWait()){throw new IllegalStateException("Incomplete input/prediction output.");}
		}
		final ByteStreamWriter summary=open(output+"/summary.tsv");
		try{
			summary.print("partition\trole\trows\tpositive\tnegative\tdistinct_scores\tsimd_ff\tsimd_fma\n");
			for(int role=0; role<2; role++){
				final ByteStreamWriter curve=open(output+"/"+ROLES[role]+".tsv");
				final long[] totals;
				try{totals=SameKModelEvaluate.writeCurve(scores[role], curve);}
				finally{if(curve.poisonAndWait()){throw new IllegalStateException("Curve write failed.");}}
				summary.print(buffer.clear().append(PARTITIONS[partition]).tab().append(ROLES[role]).tab()
						.append(scores[role].size).tab().append(totals[0]).tab().append(totals[1]).tab().append(totals[2])
						.tab().append(simdFF).tab().append(Shared.SIMD_FMA).nl());
			}
		}finally{if(summary.poisonAndWait()){throw new IllegalStateException("Summary write failed.");}}
		final String message="LIVE_BRANCH_EVALUATE_PASS partition="+PARTITIONS[partition]+
				" grey="+scores[0].size+" controls="+scores[1].size;
		final ByteStreamWriter complete=open(output+"/complete.txt");
		complete.print(message+"\n");
		if(complete.poisonAndWait()){throw new IllegalStateException("Completion marker write failed.");}
		outstream.println(message);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Requires a g-plus-four-digits genome ID and an uppercase ACGT word of length K,
	 * with K in [24,127]. The word must be lexicographically no greater than its reverse
	 * complement. Checks syntax/canonical orientation, not membership in a reference genome.
	 * @param genome Genome identifier scratch, not modified
	 * @param word DNA word scratch, not modified
	 * @param k Expected word length in bases
	 * @return Whether both identities satisfy the sidecar's syntax and orientation rules
	 */
	private static boolean validIdentity(final ByteBuilder genome, final ByteBuilder word, final int k){
		assert(genome!=null && word!=null) : "Sidecar identities are required to associate scores with upstream junctions.";
		if(genome.length()!=5 || genome.array[0]!='g' || k<24 || k>127 || word.length()!=k){return false;}
		for(int i=1; i<5; i++){if(genome.array[i]<'0' || genome.array[i]>'9'){return false;}}
		int order=0;
		for(int i=0; i<k; i++){
			final byte a=word.array[i], b=word.array[k-1-i];
			final int x=a<0 ? -1 : AminoAcid.baseToNumber[a], y=b<0 ? -1 : AminoAcid.baseToNumber[b];
			if(x<0 || y<0 || a!=AminoAcid.numberToBase[x]){return false;}
			if(order==0){order=Integer.compare(x, 3-y);}
		}
		return order<=0;
	}

	/**
	 * Starts an owned writer for a new file inside this run's output directory.
	 * @param path Output filename; append and overwrite are disabled
	 * @return Started writer whose caller must finish with poisonAndWait
	 */
	private static ByteStreamWriter open(final String path){
		final ByteStreamWriter writer=new ByteStreamWriter(path, false, false, false);
		writer.start();
		return writer;
	}

	/** Rejects a missing or byte-different UTF-8 header before any inference. */
	private static void require(final byte[] actual, final String expected){
		if(!Arrays.equals(actual, expected.getBytes(StandardCharsets.UTF_8))){throw new IllegalArgumentException("Expected header: "+expected);}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Borrowed diagnostic stream; main closes it only when it is not stdout/stderr. */
	private final PrintStream outstream;
	/** Selected features, full provenance, saved model and a new output directory. */
	private String input, rows, network, output;
	/** Required PARTITIONS index; -1 means no supported partition was selected. */
	private int partition=-1;
	/** Per-model feed-forward switch, separate from the global SIMD/FMA settings. */
	private boolean simdFF=false;

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	/** Accepted names and ordinal-counter order; the two validation partitions have distinct roles. */
	private static final String[] PARTITIONS={"train", "validation", "validation_controls"};
	/** Role indices used for retained score populations and output filenames. */
	private static final String[] ROLES={"grey", "control"};
	/** Exact seven-column sidecar header shared with the analytic fixture. */
	public static final String ROW_HEADER="partition\trow\tgenome\tk\tword\trole\ttarget";
}
