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
 * into validation candidates. Emits all attainable score>=cutoff decisions.
 * @author Fischl
 */
public final class LiveBranchModelEvaluate {

	/** Runs inference through standard BBTools parsing. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, LiveBranchModelEvaluate.class, false);
		try{new LiveBranchModelEvaluate(pp.args, pp.outstream).run();}
		finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Requires explicit partition identity, matching provenance, and new outputs. */
	private LiveBranchModelEvaluate(final String[] args, final PrintStream stream){
		outstream=stream;
		final Parser parser=new Parser();
		for(String arg : args){
			final String[] parts=arg.split("=", 2);
			String name=parts[0].toLowerCase(Locale.ROOT);
			while(name.startsWith("-")){name=name.substring(1);}
			final String value=parts.length>1 ? parts[1] : null;
			if(name.equals("rows")){rows=value;}
			else if(name.equals("net")){network=value;}
			else if(name.equals("partition")){partition=Arrays.asList(PARTITIONS).indexOf(value);}
			else if(name.equals("simdff")){simdFF=Parse.parseBoolean(value);}
			else if(!parser.parse(arg, name, value)){throw new IllegalArgumentException("Unknown evaluation option: "+arg);}
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

	/** Checks the full sidecar and selected feature stream before publishing completion. */
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
		final ByteFile meta=ByteFile.makeByteFile(rows, true), data=ByteFile.makeByteFile(input, true);
		final ByteStreamWriter predictions=open(output+"/scores.tsv");
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

	/** Validates canonical DNA and the frozen genome naming convention without per-row Strings. */
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

	/** Opens only new files inside this run's newly created output directory. */
	private static ByteStreamWriter open(final String path){
		final ByteStreamWriter writer=new ByteStreamWriter(path, false, false, false);
		writer.start();
		return writer;
	}

	/** Rejects incompatible schemas before any inference. */
	private static void require(final byte[] actual, final String expected){
		if(!Arrays.equals(actual, expected.getBytes(StandardCharsets.UTF_8))){throw new IllegalArgumentException("Expected header: "+expected);}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final PrintStream outstream;
	private String input, rows, network, output;
	private int partition=-1;
	private boolean simdFF=false;
	private static final String[] PARTITIONS={"train", "validation", "validation_controls"};
	private static final String[] ROLES={"grey", "control"};
	static final String ROW_HEADER="partition\trow\tgenome\tk\tword\trole\ttarget";
}
