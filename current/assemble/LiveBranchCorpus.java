package assemble;

import java.io.File;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Locale;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectIntMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Deduplicates an observation-only live-node trace and labels its preferred
 * route against every source-reference occurrence. Unknown truth is never a
 * negative label. Encounter multiplicity is provenance, not training weight.
 * One invocation labels one genome/library/K; genome-level folds belong to the
 * external frozen cohort, not to row-level random splitting in this tool.
 * @author Fischl
 */
public final class LiveBranchCorpus {

	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/

	/** Standard config expansion; reference truth never enters the assembler process. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, LiveBranchCorpus.class, false);
		try{new LiveBranchCorpus(pp.args, pp.outstream).run();}
		finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Requires existing trace/reference files and a fresh output directory. */
	private LiveBranchCorpus(final String[] args, final PrintStream stream){
		outstream=stream;
		final Parser parser=new Parser();
		for(String arg : args){
			final int eq=arg.indexOf('=');
			String name=(eq<0 ? arg : arg.substring(0, eq)).toLowerCase(Locale.ROOT);
			while(name.startsWith("-")){name=name.substring(1);}
			final String value=eq<0 ? null : arg.substring(eq+1);
			if(name.equals("ref")){reference=value;}
			else if(name.equals("circular")){circular=Parse.parseBoolean(value);}
			else if(name.equals("controlfraction")){controlFraction=Double.parseDouble(value);}
			else if(name.equals("controlseed")){controlSeed=Long.parseLong(value);}
			else if(!parser.parse(arg, name, value)){throw new IllegalArgumentException("Unknown live corpus option: "+arg);}
		}
		input=parser.in1; output=parser.out1;
		controlLimit(0, 0, controlFraction);// Validate even an empty corpus before creating output.
		if(input==null || reference==null || output==null || parser.in2!=null || parser.out2!=null ||
				new File(output).exists() || !Tools.testInputFiles(false, true, input, reference)){
			throw new IllegalArgumentException("Required in=trace.tsv ref=source.fa out=NEW_directory; circular=f unless explicitly known.");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Validates the complete trace before reading truth or creating any outputs. */
	private void run(){
		final ArrayList<Node> nodes=readTrace();
		Collections.sort(nodes, (a, b)->a.word.compareTo(b.word));
		final ArrayList<byte[]> queries=new ArrayList<byte[]>(nodes.size());
		for(Node node : nodes){queries.add(node.word.bases);}
		final ArrayList<byte[]> refs=new ArrayList<byte[]>();
		for(Contig contig : TadpoleGraph.readContigs(reference)){refs.add(contig.bases);}
		if(refs.isEmpty()){throw new IllegalArgumentException("Truth reference is empty.");}
		final boolean[] topology=new boolean[refs.size()];
		Arrays.fill(topology, circular);
		final KmerBranchOracle.Result[] truth=KmerBranchOracle.scan(queries, refs, topology, k);
		write(nodes, truth);
	}

	/**
	 * Retains one immutable key and count vector per node, with a reusable lookup
	 * key for duplicate encounters. Counts/flags/levels must agree on duplicates;
	 * arrival and entry provenance may differ across orientations and schedules.
	 */
	private ArrayList<Node> readTrace(){
		final ArrayList<Node> nodes=new ArrayList<Node>();
		final ObjectIntMap<Word> index=new ObjectIntMap<Word>(Word.class);
		final Word lookup=new Word();
		final ByteBuilder key=new ByteBuilder();
		final int[] values=new int[18];
		final ByteFile file=ByteFile.makeByteFile(input, true);
		final LineParser1 fields=new LineParser1('\t');
		boolean completed=false;
		int sampledEvery=0;
		try{
			require(file.nextLine(), "#livebranch-observation-v2");
			require(file.nextLine(), LiveBranchObserver.HEADER);
			for(byte[] line=file.nextLine(); line!=null; line=file.nextLine()){
				if(completed){throw new IllegalArgumentException("Data follows the trace completion record.");}
				fields.set(line);
				if(fields.termEquals("#complete", 0)){
					if(fields.terms()!=9 || fields.parseLong(1)!=encounters){throw new IllegalArgumentException("Missing/mismatched terminal encounter count.");}
					unique=fields.parseLong(2);
					observations=Double.parseDouble(fields.parseString(3));
					squares=Double.parseDouble(fields.parseString(4));
					controlEvery=fields.parseInt(5);
					final int recordedK=fields.parseInt(6);
					minExtend=fields.parseInt(7);
					fraction=Double.parseDouble(fields.parseString(8));
					if((k!=0 && k!=recordedK) || recordedK<24 || recordedK>127 || minExtend<1 ||
							controlEvery<1 || (sampledEvery!=0 && sampledEvery!=controlEvery) ||
							!Double.isFinite(fraction) || fraction<=0 || fraction>1){
						throw new IllegalArgumentException("Invalid or inconsistent trace K/mce/sampling/screen metadata.");
					}
					k=recordedK;
					completed=true;
					continue;
				}
				if(fields.terms()!=19){throw new IllegalArgumentException("A branch trace row must have19 fields.");}
				for(int i=0; i<values.length; i++){values[i]=fields.parseInt(i+1);}
				if(k==0){k=values[0];}
				if(k!=values[0] || k<24 || k>127 || values[1]<1){throw new IllegalArgumentException("Invalid node K/depth.");}
				for(int i=2; i<10; i++){if(values[i]<0){throw new IllegalArgumentException("Negative neighbor count.");}}
				for(int i : BINARY_FIELDS){
					if(values[i]!=0 && values[i]!=1){throw new IllegalArgumentException("Nonbinary branch/provenance flag.");}
				}
				if(values[16]<0 || values[16]>2 || values[17]<1){throw new IllegalArgumentException("Invalid screen or sampling denominator.");}
				if(values[14]==1){
					if(values[12]!=-1 || values[13]!=-1 || values[15]!=0){throw new IllegalArgumentException("Entry record invented an incoming edge or hidden stop.");}
				}else if(values[12]<0 || values[12]>1 || values[13]<0 || values[13]>3 ||
						values[2+4*values[12]+values[13]]<1){
					throw new IllegalArgumentException("Interior record has no supported incoming edge.");
				}
				final boolean stop=values[10]!=0 || values[11]!=0 || values[15]!=0;
				if(stop && values[17]!=1){throw new IllegalArgumentException("Native stops must be completely observed.");}
				if(!stop){
					if(sampledEvery!=0 && sampledEvery!=values[17]){throw new IllegalArgumentException("Control sampling changed within one trace.");}
					sampledEvery=values[17];
				}
				fields.appendTerm(key.clear(), 0);
				lookup.set(key.array, key.length());
				if(lookup.length!=k || !lookup.canonical()){throw new IllegalArgumentException("Expected one canonical A/C/G/T word of lengthK.");}
				final int found=index.get(lookup);
				final Node node;
				if(found<0){
					node=new Node(new Word(lookup), values);
					index.put(node.word, nodes.size()); nodes.add(node);
				}else{node=nodes.get(found); node.check(values);}
				node.observe(values);
				encounters=Math.addExact(encounters, 1);
			}
		}finally{if(file.close()){throw new IllegalStateException("Cannot read complete branch trace.");}}
		if(!completed){throw new IllegalArgumentException("Trace has no successful completion record.");}
		validateMoments(nodes);
		for(Node node : nodes){
			final int actual=LiveBranchPolicy.level(node.count, node.left, node.right, node.lb, node.rb, minExtend, fraction);
			if(actual!=node.level){throw new IllegalArgumentException("Stored branch screen disagrees with its bound configuration.");}
		}
		return nodes;
	}

	/** Rejects impossible global populations instead of normalizing corrupted data. */
	private void validateMoments(final ArrayList<Node> nodes){
		if(unique<0 || unique<nodes.size() || !Double.isFinite(observations) || !Double.isFinite(squares) ||
				observations<unique || squares<observations ||
				(unique==0 && (observations!=0 || squares!=0)) ||
				squares*unique<observations*observations*(1-1e-12)){
			throw new IllegalArgumentException("Invalid count-table population moments.");
		}
		double seen=0;
		for(Node node : nodes){seen+=node.count;}
		if(seen>observations*(1+1e-12)){throw new IllegalArgumentException("Observed unique-node depths exceed the whole table.");}
	}

	/** Writes one audited row per unique node, and at most one trainer row for it. */
	private void write(final ArrayList<Node> nodes, final KmerBranchOracle.Result[] truth){
		assert(nodes.size()==truth.length) : "Oracle results retain the exact sorted query order.";
		if(controlFraction<1){selectControls(nodes, truth);}
		if(!new File(output).mkdir()){throw new IllegalArgumentException("Output directory must be new: "+output);}
		final ByteStreamWriter[] writers=new ByteStreamWriter[OUTPUT_FILES.length];
		final long[][] census=new long[ROLES.length][3];
		final ByteBuilder row=new ByteBuilder();
		final float[] vector=new float[LiveBranchPolicy.INPUTS];
		final double uniqueMean=unique==0 ? 0 : observations/unique;
		final double observationMean=observations==0 ? 0 : squares/observations;
		long emitted=0, safe=0, unsafe=0, unknown=0, grey=0, controls=0;
		boolean error=false;
		try{
			for(int i=0; i<writers.length; i++){writers[i]=open(output+"/"+OUTPUT_FILES[i]);}
			row.append("word\tk\tlevel\trole\tlabel\ttarget\toccurrences\tunknown_flanks\tleft_mask\tright_mask\tpreferred_left\tpreferred_right\tencounters\tcontrol_every");
			for(int i=0; i<LiveBranchPolicy.INPUTS; i++){row.tab().append('f').append(i);}
			writers[0].print(row.nl());
			writers[1].print("#dims\t12\t1\n#schema=livebranch12-v1; one row per known grey-zone/control node; "+
					(controlFraction==1 ? "no balancing" : "controlfraction="+controlFraction+"; controlseed="+controlSeed)+"\n");
			writers[2].print("row\tword\trole\ttarget\n");
			for(int i=0; i<nodes.size(); i++){
				final Node node=nodes.get(i);
				final KmerBranchOracle.Result result=truth[i];
				final int left=LiveBranchPolicy.preferredBase(node.left), right=LiveBranchPolicy.preferredBase(node.right);
				final String label=label(result, left, right);
				final int target=label.equals("safe") ? 1 : label.equals("unsafe") ? 0 : -1;
				final int roleIndex=role(node, left, right);
				final String role=ROLES[roleIndex];
				census[roleIndex][target<0 ? 2 : target==1 ? 0 : 1]++;
				if(target<0){unknown++;}else if(target==1){safe++;}else{unsafe++;}
				LiveBranchPolicy.fill(vector, k, node.count, node.left, node.right, uniqueMean, observationMean);
				row.clear().append(node.word.bases).tab().append(k).tab().append(node.level).tab().append(role);
				row.tab().append(label).tab();
				if(target<0){row.append('.');}else{row.append(target);}
				row.tab().append(result.occurrences).tab().append(result.unknownFlanks).tab().append(result.leftMask).tab().append(result.rightMask);
				row.tab().append(left).tab().append(right).tab().append(node.encounters).tab().append(controlEvery);
				for(float value : vector){row.tab().appendSlow(value);}
				writers[0].print(row.nl());
				if(target>=0 && (roleIndex==0 || (roleIndex==2 && (controlFraction==1 || node.selectedControl)))){
					row.clear();
					for(float value : vector){row.appendSlow(value).tab();}
					writers[1].print(row.append(target).nl());
					emitted++;
					if(role.equals("grey")){grey++;}else{controls++;}
					writers[2].print(row.clear().append(emitted).tab().append(node.word.bases).tab().append(role).tab().append(target).nl());
				}
			}
			writeCensus(writers[3], census, nodes.size());
		}finally{for(ByteStreamWriter writer : writers){if(writer!=null){error=writer.poisonAndWait() | error;}}}
		if(error){throw new IllegalStateException("Incomplete labeled corpus output.");}
		final String summary="LIVE_BRANCH_CORPUS_PASS encounters="+encounters+" nodes="+nodes.size()+
				" safe="+safe+" unsafe="+unsafe+" unknown="+unknown+" training="+emitted+" grey="+grey+" controls="+controls+
				(controlFraction==1 ? "" : " controlfraction="+controlFraction+" controlseed="+controlSeed);
		final ByteStreamWriter complete=open(output+"/complete.txt");
		complete.print(summary+"\n");
		if(complete.poisonAndWait()){throw new IllegalStateException("Could not save corpus completion.");}
		outstream.println(summary);
	}

	/** One role predicate is shared by census, sampling and vector emission. */
	private int role(final Node node, final int left, final int right){
		assert(node.level>=0 && node.level<=2) : "Only validated policy levels may select training roles.";
		return node.level==LiveBranchPolicy.GREY_ZONE ? 0 : node.level==LiveBranchPolicy.HARD_STOP ? 1 :
				node.goodInterior && node.count>=minExtend && left>=0 && right>=0 &&
				node.left[left]>=minExtend && node.right[right]>=minExtend ? 2 : 3;
	}

	/**
	 * Takes the lowest seeded sequence hashes among known accepted controls.
	 * Safe and unsafe targets receive identical treatment; unknown truth is not
	 * eligible. Existing node objects are reused, with no boxed index collection.
	 * Only trainer rows are sampled; the full audit and role census remain intact.
	 */
	private void selectControls(final ArrayList<Node> nodes, final KmerBranchOracle.Result[] truth){
		assert(nodes.size()==truth.length) : "Sampling must use the same node/oracle pairing as vector emission.";
		final ArrayList<Node> candidates=new ArrayList<Node>();
		int grey=0;
		for(int i=0; i<nodes.size(); i++){
			final Node node=nodes.get(i);
			final int left=LiveBranchPolicy.preferredBase(node.left), right=LiveBranchPolicy.preferredBase(node.right);
			final String label=label(truth[i], left, right);
			if(!label.equals("safe") && !label.equals("unsafe")){continue;}
			final int role=role(node, left, right);
			if(role==0){grey++;}
			else if(role==2){
				node.sampleHash=sampleHash(node.word.bases, controlSeed);
				candidates.add(node);
			}
		}
		Collections.sort(candidates, (a, b)->{
			final int order=Long.compareUnsigned(a.sampleHash, b.sampleHash);
			return order!=0 ? order : a.word.compareTo(b.word);
		});
		final int limit=controlLimit(grey, candidates.size(), controlFraction);
		for(int i=0; i<limit; i++){candidates.get(i).selectedControl=true;}
	}

	/** Floor keeps controls at or below the requested fraction; one means retain all. */
	static int controlLimit(final int grey, final int available, final double fraction){
		if(grey<0 || available<0 || !Double.isFinite(fraction) || fraction<0 || fraction>1){
			throw new IllegalArgumentException("controlfraction must be finite and in [0,1]; population counts must be nonnegative.");
		}
		return fraction==1 ? available : (int)Math.min(available, Math.floor(grey*fraction/(1-fraction)));
	}

	/** Stable sequence-only rank; callers use a frozen per-genome seed, never a target-derived seed. */
	static long sampleHash(final byte[] word, final long seed){
		assert(word!=null && word.length>0) : "A control identity is its complete canonical junction word.";
		long hash=seed^word.length;
		for(byte base : word){hash=Tools.murmur64(hash^base);}
		return hash;
	}

	/** Separates screened candidates from controls; totals never stand in for candidate error rates. */
	private static void writeCensus(final ByteStreamWriter writer, final long[][] census, final int nodes){
		final ByteBuilder row=new ByteBuilder();
		writer.print("role\tsafe\tunsafe\tunknown\n");
		long total=0;
		for(int i=0; i<ROLES.length; i++){
			row.clear().append(ROLES[i]);
			for(long count : census[i]){row.tab().append(count); total+=count;}
			writer.print(row.nl());
		}
		assert(total==nodes) : "Every deduplicated node must contribute exactly once to the role/label census: "+total+" != "+nodes;
	}

	/** Strict local truth, not a guarantee about long-range repeat resolution or an entire assembly. */
	static String label(final KmerBranchOracle.Result result, final int left, final int right){
		if(result==null){throw new IllegalArgumentException("Reference census is required.");}
		if(left< -1 || left>3 || right< -1 || right>3){throw new IllegalArgumentException("Preferred bases must be A/C/G/T codes or -1.");}
		if(left<0 || right<0){return "no_route";}
		if(result.occurrences==0){return "absent";}
		if(result.unknownFlanks>0){return "boundary";}
		return result.unambiguous(left, right) ? "safe" : "unsafe";
	}

	/** Header comparison is outside the row hot path. */
	private static void require(final byte[] line, final String expected){
		if(line==null || !expected.equals(new String(line, StandardCharsets.US_ASCII))){throw new IllegalArgumentException("Wrong live trace header.");}
	}

	/** New native text writer; every caller must observe poisonAndWait. */
	private static ByteStreamWriter open(final String path){
		final ByteStreamWriter writer=new ByteStreamWriter(path, false, false, false);
		writer.start(); return writer;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** One persistent node; occurrence count never changes its model input or weight. */
	private static final class Node {
		Node(final Word word_, final int[] v){
			word=word_; count=v[1]; left=Arrays.copyOfRange(v, 2, 6); right=Arrays.copyOfRange(v, 6, 10);
			lb=v[10]!=0; rb=v[11]!=0; level=v[16];
		}
		/** Immutable graph evidence must agree across every encounter of one canonical word. */
		void check(final int[] v){
			if(count!=v[1] || lb!=(v[10]!=0) || rb!=(v[11]!=0) || level!=v[16]){
				throw new IllegalArgumentException("Conflicting duplicate junction evidence.");
			}
			for(int i=0; i<4; i++){
				if(left[i]!=v[2+i] || right[i]!=v[6+i]){throw new IllegalArgumentException("Conflicting duplicate neighbor counts.");}
			}
		}
		/** Only a nonhidden interior traversal of the preferred route supplies an accepted control. */
		void observe(final int[] v){
			encounters=Math.addExact(encounters, 1);
			if(v[14]==0 && v[15]==0){
				final int preferred=LiveBranchPolicy.preferredBase(v[12]==0 ? left : right);
				goodInterior|=preferred>=0 && preferred==v[13];
			}
		}
		final Word word;
		final int count, level;
		final int[] left, right;
		final boolean lb, rb;
		long encounters;
		boolean goodInterior;
		long sampleHash;
		boolean selectedControl;
	}

	/** Exact byte key; mutable lookup is never inserted, persistent copies never change. */
	private static final class Word implements Comparable<Word> {
		Word(){}
		Word(final Word other){set(Arrays.copyOf(other.bases, other.length), other.length);}
		void set(final byte[] bytes, final int len){
			assert(len>=0 && len<=bytes.length) : "Key extent must lie inside the reusable parser buffer.";
			bases=bytes; length=len; hash=1;
			for(int i=0; i<len; i++){hash=31*hash+bytes[i];}
		}
		/** Checks the supplied canonical orientation instead of silently repairing bad records. */
		boolean canonical(){
			int comparison=0;
			for(int i=0; i<length; i++){
				final int a=bases[i]<0 ? -1 : AminoAcid.baseToNumber[bases[i]];
				final int reverse=bases[length-1-i]<0 ? -1 : AminoAcid.baseToNumber[bases[length-1-i]];
				if(a<0 || reverse<0 || bases[i]!=AminoAcid.numberToBase[a]){return false;}
				if(comparison==0){comparison=Integer.compare(a, 3-reverse);}
			}
			return comparison<=0;
		}
		@Override public int hashCode(){return hash;}
		@Override public boolean equals(final Object other){
			if(!(other instanceof Word)){return false;}
			final Word word=(Word)other;
			if(length!=word.length){return false;}
			for(int i=0; i<length; i++){if(bases[i]!=word.bases[i]){return false;}}
			return true;
		}
		@Override public int compareTo(final Word other){
			for(int i=0; i<Math.min(length, other.length); i++){
				if(bases[i]!=other.bases[i]){return Byte.compare(bases[i], other.bases[i]);}
			}
			return Integer.compare(length, other.length);
		}
		byte[] bases;
		int length, hash;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final PrintStream outstream;
	private final String input, output;
	private String reference;
	private boolean circular;
	private int k, minExtend, controlEvery;
	private double fraction, observations, squares;
	private long encounters, unique;
	private double controlFraction=1;
	private long controlSeed=1;
	private static final int[] BINARY_FIELDS={10, 11, 14, 15};
	private static final String[] ROLES={"grey", "hard", "control", "excluded"};
	private static final String[] OUTPUT_FILES={"nodes.tsv", "vectors.tsv", "rows.tsv", "census.tsv"};
}
