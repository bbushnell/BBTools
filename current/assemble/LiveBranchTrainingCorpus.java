package assemble;

import java.io.File;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectIntMap;
import parse.LineParser1;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Joins frozen per-genome live-branch exports without splitting individual rows.
 * Training retains the upstream sampled controls; validation candidates and
 * controls have separate tables so abundant controls cannot dominate checkpoint
 * selection. The upstream frozen cohort defines folds; this tool does not prove
 * homology isolation or generate reference truth.
 * @author Fischl
 */
public final class LiveBranchTrainingCorpus {

	/** Standard config expansion and native output-stream lifecycle. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, LiveBranchTrainingCorpus.class, false);
		try{new LiveBranchTrainingCorpus(pp.args, pp.outstream).run();}
		finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Requires the selected manifest, expected genome count and a new output directory. */
	private LiveBranchTrainingCorpus(final String[] args, final PrintStream stream){
		outstream=stream;
		final Parser parser=new Parser();
		for(String arg : args){
			final int eq=arg.indexOf('=');
			String name=(eq<0 ? arg : arg.substring(0, eq)).toLowerCase(Locale.ROOT);
			while(name.startsWith("-")){name=name.substring(1);}
			final String value=eq<0 ? null : arg.substring(eq+1);
			if(name.equals("expected")){expected=Integer.parseInt(value);}
			else if(!parser.parse(arg, name, value)){throw new IllegalArgumentException("Unknown live training option: "+arg);}
		}
		input=parser.in1; output=parser.out1;
		if(input==null || output==null || expected<1 || parser.in2!=null || parser.out2!=null ||
				new File(output).exists() || !Tools.testInputFiles(false, true, input)){
			throw new IllegalArgumentException("Required in=selected.tsv expected=848 out=NEW_directory.");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Validates all small routing/terminal records, then streams each vector/provenance pair once. */
	private void run(){
		final ArrayList<Sample> samples=readManifest();
		for(Sample sample : samples){preflight(sample);}
		if(!new File(output).mkdir()){throw new IllegalArgumentException("Cannot create output directory: "+output);}
		final ByteStreamWriter[] writers=new ByteStreamWriter[5];
		final long[] emitted=new long[3];
		boolean error=false;
		try{
			for(int i=0; i<3; i++){
				writers[i]=open(output+"/"+PARTITIONS[i]+".tsv");
				writers[i].print("#dims\t12\t1\n"+SCHEMA+"\n");
			}
			writers[3]=open(output+"/rows.tsv");
			writers[3].print("partition\trow\tgenome\tk\tword\trole\ttarget\n");
			writers[4]=open(output+"/census.tsv");
			writers[4].print("genome\tfold\tk\tgrey_safe\tgrey_unsafe\tgrey_unknown\tcontrol_safe_available\tcontrol_unsafe_available\tcontrol_unknown\tgrey_emitted\tcontrols_emitted\n");
			for(Sample sample : samples){convert(sample, writers, emitted);}
		}finally{
			for(ByteStreamWriter writer : writers){if(writer!=null){error=writer.poisonAndWait() | error;}}
		}
		if(error){throw new IllegalStateException("Incomplete live training output.");}
		final String summary="LIVE_BRANCH_TRAINING_PASS genomes="+samples.size()+" train="+emitted[0]+
				" validation="+emitted[1]+" validation_controls="+emitted[2];
		final ByteStreamWriter complete=open(output+"/complete.txt");
		complete.print(summary+"\n");
		if(complete.poisonAndWait()){throw new IllegalStateException("Cannot save training completion.");}
		outstream.println(summary);
	}

	/** The hash-bound upstream manifest is authoritative; TEST is deliberately not accepted here. */
	private ArrayList<Sample> readManifest(){
		final ArrayList<Sample> samples=new ArrayList<Sample>();
		final ObjectIntMap<String> ids=new ObjectIntMap<String>(String.class), dirs=new ObjectIntMap<String>(String.class);
		final ByteFile file=ByteFile.makeByteFile(input, true);
		final LineParser1 fields=new LineParser1('\t');
		try{
			require(file.nextLine(), MANIFEST_HEADER);
			for(byte[] line=file.nextLine(); line!=null; line=file.nextLine()){
				fields.set(line);
				if(fields.terms()!=10){throw new IllegalArgumentException("Expected ten selected-manifest columns.");}
				final String id=fields.parseString(0), fold=fields.parseString(1), action=fields.parseString(5);
				final int k=fields.parseInt(2), depth=fields.parseInt(3);
				final long seed=fields.parseLong(4);
				final String run=fields.parseString(6), labels=fields.parseString(7);
				if(!id.matches("g[0-9]{4}") || (!fold.equals("train") && !fold.equals("validation")) ||
						k<24 || k>127 || depth<15 || depth>150 || seed<0 ||
						(!action.equals("new") && !action.equals("reuse")) ||
						!new File(run).isAbsolute() || !new File(labels).isAbsolute() ||
						ids.get(id)>=0 || dirs.get(labels)>=0){
					throw new IllegalArgumentException("Invalid/duplicate selected genome or unsupported fold: "+id);
				}
				ids.put(id, samples.size()); dirs.put(labels, samples.size());
				samples.add(new Sample(id, fold.equals("train"), k, seed, run, labels));
			}
		}finally{if(file.close()){throw new IllegalStateException("Cannot read selected manifest.");}}
		if(samples.size()!=expected){throw new IllegalArgumentException("Selected genome count "+samples.size()+" != "+expected);}
		return samples;
	}

	/** Binds full census to the exact TRAIN sampling policy or unsampled VALID export. */
	private static void preflight(final Sample sample){
		if(!Tools.testInputFiles(false, true, sample.run+"/exit_code.txt", sample.labels+"/vectors.tsv", sample.labels+"/rows.tsv")){
			throw new IllegalArgumentException("Missing completed source files for "+sample.id);
		}
		require(singleLine(sample.run+"/exit_code.txt"), "0");
		final ByteFile census=ByteFile.makeByteFile(sample.labels+"/census.tsv", true);
		final LineParser1 fields=new LineParser1('\t');
		try{
			require(census.nextLine(), "role\tsafe\tunsafe\tunknown");
			for(int role=0; role<ROLES.length; role++){
				final byte[] row=census.nextLine();
				if(row==null){throw new IllegalArgumentException("Truncated role census.");}
				fields.set(row);
				if(fields.terms()!=4 || !fields.termEquals(ROLES[role], 0)){throw new IllegalArgumentException("Wrong role census order/width.");}
				for(int col=0; col<3; col++){
					final long value=fields.parseLong(col+1);
					if(value<0){throw new IllegalArgumentException("Negative role count.");}
					sample.census[role][col]=value;
				}
			}
			if(census.nextLine()!=null){throw new IllegalArgumentException("Extra role census row.");}
		}finally{if(census.close()){throw new IllegalStateException("Cannot read role census.");}}
		final String[] terminal=new String(singleLine(sample.labels+"/complete.txt"), StandardCharsets.US_ASCII).split(" ");
		if(terminal.length!=(sample.train ? 11 : 9) || !terminal[0].equals("LIVE_BRANCH_CORPUS_PASS")){
			throw new IllegalArgumentException("Wrong completion schema or sampled validation export: "+sample.id);
		}
		final long encounters=number(terminal[1], "encounters="), nodes=number(terminal[2], "nodes=");
		final long[] totals=new long[3];
		for(long[] row : sample.census){for(int i=0; i<3; i++){totals[i]=Math.addExact(totals[i], row[i]);}}
		if(encounters<nodes || nodes!=Math.addExact(Math.addExact(totals[0], totals[1]), totals[2]) ||
				number(terminal[3], "safe=")!=totals[0] || number(terminal[4], "unsafe=")!=totals[1] ||
				number(terminal[5], "unknown=")!=totals[2]){throw new IllegalArgumentException("Completion and full census disagree.");}
		sample.grey=number(terminal[7], "grey="); sample.controls=number(terminal[8], "controls=");
		final long available=Math.addExact(sample.census[2][0], sample.census[2][1]);
		final long expectedControls=sample.train ? LiveBranchCorpus.controlLimit(Math.toIntExact(sample.grey), Math.toIntExact(available), 0.1) : available;
		if(sample.grey!=Math.addExact(sample.census[0][0], sample.census[0][1]) || sample.controls!=expectedControls ||
				number(terminal[6], "training=")!=Math.addExact(sample.grey, sample.controls)){
			throw new IllegalArgumentException("Completion disagrees with grey population or control policy: "+sample.id);
		}
		if(sample.train && (!terminal[9].equals("controlfraction=0.1") || number(terminal[10], "controlseed=")!=sample.seed)){
			throw new IllegalArgumentException("Training control seed/fraction disagrees with the frozen manifest.");
		}
	}

	/** Streams rows and values together; sorted canonical identities prevent duplicate node weighting. */
	private static void convert(final Sample sample, final ByteStreamWriter[] writers, final long[] emitted){
		final ByteFile data=ByteFile.makeByteFile(sample.labels+"/vectors.tsv", true);
		final ByteFile meta=ByteFile.makeByteFile(sample.labels+"/rows.tsv", true);
		final LineParser1 dp=new LineParser1('\t'), mp=new LineParser1('\t');
		final ByteBuilder word=new ByteBuilder(), previous=new ByteBuilder(), out=new ByteBuilder();
		final long[][] counts=new long[2][2];
		long ordinal=0;
		try{
			require(data.nextLine(), "#dims\t12\t1");
			require(data.nextLine(), "#schema=livebranch12-v1; one row per known grey-zone/control node; "+
					(sample.train ? "controlfraction=0.1; controlseed="+sample.seed : "no balancing"));
			require(meta.nextLine(), "row\tword\trole\ttarget");
			final float encodedK=(float)(0.125*Math.log(sample.k)*(1/Math.log(2)));
			for(byte[] line=meta.nextLine(); line!=null; line=meta.nextLine()){
				mp.set(line);
				if(mp.terms()!=4 || mp.parseLong(0)!=++ordinal){throw new IllegalArgumentException("Wrong provenance width/ordinal at "+sample.id);}
				final int role=mp.termEquals("grey", 2) ? 0 : mp.termEquals("control", 2) ? 1 : -1;
				final int target=mp.parseInt(3);
				mp.appendTerm(word.clear(), 1);
				if(role<0 || (target!=0 && target!=1) || !canonical(word, sample.k) || compare(previous, word)>=0){
					throw new IllegalArgumentException("Invalid/duplicate canonical node, role or target at "+sample.id);
				}
				previous.clear().append(word);
				final byte[] values=data.nextLine();
				if(values==null){throw new IllegalArgumentException("Missing paired vector at "+sample.id);}
				dp.set(values);
				if(dp.terms()!=13 || dp.parseInt(12)!=target){throw new IllegalArgumentException("Wrong vector width or target/provenance mismatch.");}
				for(int i=0; i<12; i++){
					final float value=dp.parseFloat(i);
					if(!Float.isFinite(value) || value<0){throw new IllegalArgumentException("Invalid live feature at "+sample.id);}
				}
				if(Math.abs(dp.parseFloat(9)-encodedK)>1e-6){throw new IllegalArgumentException("Feature K disagrees with the selected genome.");}
				counts[role][target]++;
				final int partition=sample.train ? 0 : role==0 ? 1 : 2;
				writers[partition].print(out.clear().append(values).nl());
				writers[3].print(out.clear().append(PARTITIONS[partition]).tab().append(++emitted[partition]).tab().append(sample.id)
						.tab().append(sample.k).tab().append(word).tab().append(ROLES[role==0 ? 0 : 2]).tab().append(target).nl());
			}
			if(data.nextLine()!=null || counts[0][0]!=sample.census[0][1] || counts[0][1]!=sample.census[0][0] ||
					counts[1][0]+counts[1][1]!=sample.controls || counts[1][0]>sample.census[2][1] || counts[1][1]>sample.census[2][0]){
				throw new IllegalArgumentException("Incomplete inputs or label counts disagree with the role census at "+sample.id);
			}
		}finally{if(data.close() | meta.close()){throw new IllegalStateException("Incomplete live corpus input.");}}
		out.clear().append(sample.id).tab().append(sample.train ? "train" : "validation").tab().append(sample.k);
		for(long count : sample.census[0]){out.tab().append(count);}
		for(long count : sample.census[2]){out.tab().append(count);}
		writers[4].print(out.tab().append(sample.grey).tab().append(sample.controls).nl());
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Exact canonical DNA identity, independent of the model's numeric features. */
	private static boolean canonical(final ByteBuilder word, final int k){
		if(word.length()!=k){return false;}
		int order=0;
		for(int i=0; i<k; i++){
			final byte a=word.array[i], b=word.array[k-1-i];
			final int x=a<0 ? -1 : AminoAcid.baseToNumber[a], y=b<0 ? -1 : AminoAcid.baseToNumber[b];
			if(x<0 || y<0 || a!=AminoAcid.numberToBase[x]){return false;}
			if(order==0){order=Integer.compare(x, 3-y);}
		}
		return order<=0;
	}

	/** Strict lexical order is the exporter's persistent unique-node order. */
	private static int compare(final ByteBuilder a, final ByteBuilder b){
		assert(a!=null && b!=null) : "Streaming duplicate detection requires both retained identities.";
		for(int i=0; i<Math.min(a.length(), b.length()); i++){
			if(a.array[i]!=b.array[i]){return Byte.compare(a.array[i], b.array[i]);}
		}
		return Integer.compare(a.length(), b.length());
	}

	/** A terminal record must be complete and singular, not an arbitrary log prefix. */
	private static byte[] singleLine(final String path){
		final ByteFile file=ByteFile.makeByteFile(path, true);
		try{
			final byte[] line=file.nextLine();
			if(line==null || file.nextLine()!=null){throw new IllegalArgumentException("Expected one terminal record: "+path);}
			return line;
		}finally{if(file.close()){throw new IllegalStateException("Cannot read terminal record: "+path);}}
	}

	/** Named metadata fields are cold; all per-vector parsing remains byte-oriented. */
	private static long number(final String field, final String prefix){
		if(!field.startsWith(prefix)){throw new IllegalArgumentException("Expected field "+prefix);}
		final long value=Long.parseLong(field.substring(prefix.length()));
		if(value<0){throw new IllegalArgumentException("Negative count/seed: "+field);}
		return value;
	}

	/** Rejects similarly shaped but incompatible upstream formats. */
	private static void require(final byte[] line, final String expected){
		if(!Arrays.equals(line, expected.getBytes(StandardCharsets.US_ASCII))){throw new IllegalArgumentException("Expected header/record: "+expected);}
	}

	/** All paths are beneath a newly created output directory. */
	private static ByteStreamWriter open(final String path){
		final ByteStreamWriter writer=new ByteStreamWriter(path, false, false, false);
		writer.start(); return writer;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Cold routing metadata and full role counts for one genome, never retained vector rows. */
	private static final class Sample {
		/** Binds the immutable selected-manifest route before reading any vector. */
		Sample(final String id_, final boolean train_, final int k_, final long seed_, final String run_, final String labels_){
			id=id_; run=run_; labels=labels_;
			train=train_; k=k_; seed=seed_;
		}
		final String id, run, labels;
		final boolean train;
		final int k;
		final long seed;
		final long[][] census=new long[4][3];
		long grey, controls;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final PrintStream outstream;
	private final String input, output;
	private int expected=848;
	static final String SCHEMA="#schema=livebranch12-v1; genome folds fixed; validation candidates and controls separated";
	static final String MANIFEST_HEADER="genome\tfold\tk\tdepth\tseed\taction\trun\tlabels\treads\treference";
	private static final String[] PARTITIONS={"train", "validation", "validation_controls"};
	private static final String[] ROLES={"grey", "hard", "control", "excluded"};
}
