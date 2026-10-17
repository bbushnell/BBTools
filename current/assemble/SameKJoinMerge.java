package assemble;

import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.LongHashMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Development-only reference-guided or NN-scored joining of saved raw unitigs.
 * The two decision inputs are mutually exclusive; NN mode never opens labels.
 * No read table is loaded, and no unknown label is converted into a negative.
 * @author Fischl
 */
public final class SameKJoinMerge {

	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/

	/** Uses the ordinary preparser and closes redirected diagnostic output. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, SameKJoinMerge.class, false);
		try{new SameKJoinMerge(pp.args, pp.outstream).run();}
		finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Keeps scientific choices explicit; all outputs must be new. */
	private SameKJoinMerge(final String[] args, final PrintStream outstream_){
		outstream=outstream_;
		final Parser parser=new Parser();
		for(String arg : args){
			final String[] split=arg.split("=", 2);
			String name=split[0].toLowerCase(Locale.ROOT);
			while(name.startsWith("-")){name=name.substring(1);}
			final String value=split.length>1 ? split[1] : null;
			if(name.equals("candidates")){candidates=value;}
			else if(name.equals("labels")){labels=value;}
			else if(name.equals("features")){features=value;}
			else if(name.equals("net")){network=value;}
			else if(name.equals("outscore")){scores=value;}
			else if(name.equals("cutoff")){cutoff=Float.parseFloat(value); cutoffSpecified=true;}
			else if(name.equals("outbaseline")){baseline=value;}
			else if(name.equals("outpaths")){paths=value;}
			else if(name.equals("k")){k=Parse.parseIntKMG(value);}
			else if(name.equals("mincontig")){minContig=Parse.parseIntKMG(value);}
			else if(!parser.parse(arg, name, value)){throw new IllegalArgumentException("Unknown merge argument: "+arg);}
		}
		input=parser.in1; output=parser.out1;
		final boolean modelMode=features!=null || network!=null || scores!=null || cutoffSpecified;
		if(modelMode && (labels!=null || features==null || network==null || scores==null || !Float.isFinite(cutoff) || cutoff<0 || cutoff>1)){
			throw new IllegalArgumentException("Use labels= OR features=/net=/cutoff=/outscore=; never both decision sources.");
		}
		if(input==null || output==null || candidates==null || (!modelMode && labels==null) || baseline==null || paths==null
				|| parser.in2!=null || parser.out2!=null || k<2 || minContig<1){
			throw new IllegalArgumentException("Required: in=raw.fa candidates=table labels=table OR features=/net=/cutoff=/outscore=, k=N out=merged.fa outbaseline=raw.filtered.fa outpaths=paths.tsv mincontig>=1");
		}
		final ArrayList<String> inputs=new ArrayList<String>(Arrays.asList(input, candidates));
		if(modelMode){inputs.add(features); inputs.add(network);}else{inputs.add(labels);}
		TadpoleGraph.checkPaths(input, inputs, output, baseline, paths, scores);
		if(!Tools.testInputFiles(false, true, inputs.toArray(new String[0]))
				|| !Tools.testOutputFiles(false, false, false, output, baseline, paths, scores)){
			throw new IllegalArgumentException("Inputs must exist and all merge outputs must be new writable files.");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Validates all evidence before writing, then emits matched final-length filters. */
	private void run(){
		assert(k>=2 && minContig>0) : "The overlap and final filter are explicit experimental parameters.";
		final ArrayList<Contig> contigs=TadpoleGraph.readContigs(input);
		for(int i=0; i<contigs.size(); i++){
			final String expected="unitig_"+i+" ";
			if(contigs.get(i).name==null || !contigs.get(i).name.startsWith(expected)){
				throw new IllegalArgumentException("Raw contigs must retain census order and unitig IDs; expected "+expected);
			}
		}
		final Choices choices=readChoices(contigs);
		final Result result=merge(contigs, choices.mates(), k);
		final int emitted=writeFasta(output, result.contigs, minContig);
		final int before=writeFasta(baseline, contigs, minContig);
		writePaths(result);
		outstream.println((network==null ? "SAME_K_MERGE_PASS input=" : "SAME_K_NN_MERGE_PASS input=")+contigs.size()+
				(network==null ? " supported=" : " selected=")+choices.supported+
				" joined="+result.joins+" cycle_unitigs="+result.cycleUnitigs+" output="+result.contigs.size()+
				" baseline_emitted="+before+" merged_emitted="+emitted+" mincontig="+minContig);
	}

	/** Reads consecutive paired rows and verifies the candidate against the actual FASTA. */
	private Choices readChoices(final ArrayList<Contig> contigs){
		assert(contigs!=null) : "Candidate indices refer to one fixed ordered raw-unitig set.";
		final Choices choices=new Choices(contigs.size());
		final LongHashMap seen=new LongHashMap();
		final ByteFile candidateFile=ByteFile.makeByteFile(candidates, true);
		ByteFile labelFile=null;
		SameKJoinModelGate gate=null;
		final LineParser1 cp=new LineParser1('\t'), lp=new LineParser1('\t');
		int event=0;
		boolean readError=false;
		try{
			requireHeader(candidateFile.nextLine(), SameKJoinDiagnostic.HEADER);
			if(network==null){
				labelFile=ByteFile.makeByteFile(labels, true);
				final byte[] comment=labelFile.nextLine();
				if(comment==null || comment.length<1 || comment[0]!='#'){throw new IllegalArgumentException("Missing label schema comment.");}
				requireHeader(labelFile.nextLine(), LABEL_HEADER);
			}else{gate=new SameKJoinModelGate(features, network, scores, k, cutoff);}
			for(byte[] row=candidateFile.nextLine(); row!=null; row=candidateFile.nextLine()){
				cp.set(row); event++;
				if(labelFile!=null){
					final byte[] label=labelFile.nextLine();
					if(label==null){throw new IllegalArgumentException("Missing label for candidate "+event);}
					lp.set(label);
				}
				if(cp.terms()!=36 || cp.parseInt(0)!=event || cp.parseInt(1)!=k ||
						(labelFile!=null && (lp.terms()!=17 || lp.parseInt(0)!=event || lp.parseInt(1)!=k))){
					throw new IllegalArgumentException("Mismatched candidate/label schema, K or event ID: "+event);
				}
				final int source=cp.parseInt(2), sourceSide=cp.parseInt(3), dest=cp.parseInt(4), destSide=cp.parseInt(5);
				if(source<0 || dest<0 || source>=contigs.size() || dest>=contigs.size()
						|| (sourceSide&~1)!=0 || (destSide&~1)!=0){throw new IllegalArgumentException("Invalid endpoint at event "+event);}
				final int se=2*source+sourceSide, de=2*dest+destSide;
				final long pair=((long)se<<32)|(de&0xFFFFFFFFL);
				if(se>=de || seen.contains(pair)){throw new IllegalArgumentException("Duplicate or unordered physical pair at event "+event);}
				seen.put(pair, 1);
				final Contig a=contigs.get(source), b=contigs.get(dest);
				final boolean ar=sourceSide==0, br=destSide==1;
				if(cp.parseInt(6)!=a.length() || cp.parseInt(7)!=b.length() || a.length()<k || b.length()<k){
					throw new IllegalArgumentException("Raw contig/candidate length mismatch at event "+event);
				}
				verifyText(a, ar, a.length()-k, cp.parseByteArray(12));
				verifyText(b, br, 0, cp.parseByteArray(13));
				final byte[] left=cp.parseByteArray(34), right=cp.parseByteArray(35);
				if(left.length<k || right.length<k){throw new IllegalArgumentException("Context must include both non-overlap arms.");}
				verifyText(a, ar, a.length()-left.length, left);
				verifyText(b, br, 0, right);
				verifyOverlap(a, ar, b, br, k-1);
				if(gate!=null){
					if(gate.accept(event)){choices.add(se, de);}
				}else if(lp.termEquals("supported", 2)){
					if(!lp.termEquals("1", 3) || lp.parseInt(4)!=1 || lp.parseInt(5)!=1 || lp.parseInt(6)!=1){
						throw new IllegalArgumentException("Supported label lacks unique exact context at event "+event);
					}
					choices.add(se, de);
				}else if(lp.termEquals("contradicted", 2)){
					if(!lp.termEquals("0", 3)){throw new IllegalArgumentException("Contradicted target must be zero.");}
				}else if(!lp.termEquals("unresolved", 2) || !lp.termEquals("NA", 3)){
					throw new IllegalArgumentException("Unknown label/target at event "+event);
				}
			}
			if(labelFile!=null && labelFile.nextLine()!=null){throw new IllegalArgumentException("Extra labels after final candidate.");}
			if(gate!=null){gate.finish();}
		}finally{
			readError=candidateFile.close();
			if(labelFile!=null){readError|=labelFile.close();}
			if(gate!=null){readError|=gate.close();}
		}
		if(readError){throw new RuntimeException("Error reading candidate/label evidence.");}
		return choices;
	}

	/** Keeps the full component trace, including components filtered from final FASTA. */
	private void writePaths(final Result result){
		assert(result.contigs.size()==result.paths.size()) : "Every spelled component has exactly one provenance path.";
		final ByteStreamWriter writer=new ByteStreamWriter(paths, false, false, false);
		writer.start();
		final ByteBuilder row=new ByteBuilder();
		try{
			writer.print("contig\tlength\temitted\tunitig_path\n");
			for(int i=0; i<result.contigs.size(); i++){
				final Contig c=result.contigs.get(i);
				writer.print(row.clear().append("samek_").append(i).tab().append(c.length()).tab()
						.append(c.length()>=minContig ? 1 : 0).tab().append(result.paths.get(i)).nl());
			}
		}finally{if(writer.poisonAndWait()){throw new RuntimeException("Path output failed.");}}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Spells disjoint linear chains; leaves every closed component unmerged. */
	static Result merge(final ArrayList<Contig> input, final int[] mate, final int k){
		if(k<2 || mate.length!=2L*input.size()){throw new IllegalArgumentException("Invalid merge geometry.");}
		for(int end=0; end<mate.length; end++){
			final int other=mate[end];
			if(other< -1 || other>=mate.length || (other>=0 && (other/2==end/2 || mate[other]!=end))){
				throw new IllegalArgumentException("Mates must be reciprocal connections between distinct unitigs.");
			}
		}
		final Result result=new Result();
		final boolean[] used=new boolean[input.size()];
		for(int end=0; end<mate.length; end++){
			if(mate[end]<0 && !used[end/2]){spell(input, mate, k, end, used, result);}
		}
		for(int id=0; id<input.size(); id++){
			if(used[id]){continue;}
			// A degree-two component with no free end is a cycle. Preserve original
			// unitigs rather than opening it arbitrarily or duplicating its closure.
			result.cycleUnitigs++;
			final Contig original=input.get(id);
			copyOuterEnds(result.add(original.bases.clone(), id+"+"), original, false, original, false);
			used[id]=true;
		}
		assert(result.contigs.size()+result.joins==input.size()) : "Each linear join absorbs one unitig; cycles absorb none.";
		return result;
	}

	/** Every contig is appended once; each exact K-1 overlap is emitted once. */
	private static void spell(final ArrayList<Contig> input, final int[] mate, final int k,
			final int start, final boolean[] used, final Result result){
		assert(mate[start]<0) : "Only a free physical end starts a linear component.";
		final ByteBuilder bases=new ByteBuilder(), path=new ByteBuilder();
		int entry=start, count=0;
		Contig previous=null;
		boolean previousReverse=false;
		while(entry>=0){
			final int id=entry/2;
			final Contig c=input.get(id);
			final boolean reverse=(entry&1)==1;
			if(used[id]){throw new IllegalStateException("A linear component revisited unitig "+id);}
			used[id]=true;
			final int skip=count==0 ? 0 : k-1;
			if(previous!=null){verifyOverlap(previous, previousReverse, c, reverse, skip);}
			if((long)bases.length+c.length()-skip>Integer.MAX_VALUE-8L){throw new IllegalArgumentException("Joined contig exceeds Java array capacity.");}
			for(int p=skip; p<c.length(); p++){bases.append(SameKJoinCollector.baseAt(c, reverse, p));}
			if(count>0){path.append(','); result.joins++;}
			path.append(id).append(reverse ? '-' : '+'); count++;
			previous=c; previousReverse=reverse;
			entry=mate[entry^1];
		}
		copyOuterEnds(result.add(bases.toBytes(), path.toString()), input.get(start/2), (start&1)==1,
				previous, previousReverse);
	}

	/**
	 * Exact joining removes internal ends only. Preserve the physical outer ends
	 * for later bridge/fusion eligibility; default integer zero is not topology.
	 * Coverage is intentionally recomputed from the caller's live count table.
	 */
	private static void copyOuterEnds(final Contig target, final Contig first, final boolean firstReverse,
			final Contig last, final boolean lastReverse){
		assert(target!=null && first!=null && last!=null) : "Every spelled path has two surviving physical outer ends.";
		target.leftCode=firstReverse ? first.rightCode : first.leftCode;
		target.leftRatio=firstReverse ? first.rightRatio : first.leftRatio;
		target.leftBridgeEndpoint=firstReverse ? first.rightBridgeEndpoint : first.leftBridgeEndpoint;
		target.leftFusionEndpoint=firstReverse ? first.rightFusionEndpoint : first.leftFusionEndpoint;
		target.graphLeftLimited=firstReverse ? first.graphRightLimited : first.graphLeftLimited;
		target.rightCode=lastReverse ? last.leftCode : last.rightCode;
		target.rightRatio=lastReverse ? last.leftRatio : last.rightRatio;
		target.rightBridgeEndpoint=lastReverse ? last.leftBridgeEndpoint : last.rightBridgeEndpoint;
		target.rightFusionEndpoint=lastReverse ? last.leftFusionEndpoint : last.rightFusionEndpoint;
		target.graphRightLimited=lastReverse ? last.graphLeftLimited : last.graphRightLimited;
		target.tid=first.tid==last.tid ? first.tid : -1;
	}

	/** Case/orientation-normalized sequence equality is mandatory, not merely a hash hit. */
	private static void verifyOverlap(final Contig a, final boolean ar, final Contig b, final boolean br, final int overlap){
		if(a.length()<=overlap || b.length()<=overlap){throw new IllegalArgumentException("Overlap consumes a whole unitig.");}
		for(int i=0; i<overlap; i++){
			if(SameKJoinCollector.baseAt(a, ar, a.length()-overlap+i)!=SameKJoinCollector.baseAt(b, br, i)){
				throw new IllegalArgumentException("Proposed unitigs do not share the exact K-1 overlap.");
			}
		}
	}

	/** Binds saved diagnostic contexts to the original ordered FASTA, not just IDs. */
	private static void verifyText(final Contig c, final boolean reverse, final int start, final byte[] text){
		if(start<0 || text.length<1 || (long)start+text.length>c.length()){throw new IllegalArgumentException("Context does not fit raw unitig.");}
		for(int i=0; i<text.length; i++){
			if(SameKJoinCollector.baseAt(c, reverse, start+i)!=text[i]){throw new IllegalArgumentException("Candidate context differs from raw unitig.");}
		}
	}

	/** Fixed schemas prevent a shifted column from becoming a false oracle decision. */
	private static void requireHeader(final byte[] actual, final String expected){
		if(actual==null || !expected.equals(new String(actual, StandardCharsets.US_ASCII))){throw new IllegalArgumentException("Unexpected evidence header.");}
	}

	/** Applies the same final length cutoff to both baseline and merged products. */
	private static int writeFasta(final String path, final ArrayList<Contig> contigs, final int min){
		assert(path!=null && min>0) : "Both output arms use an explicit identical minimum length.";
		final ByteStreamWriter writer=new ByteStreamWriter(path, false, false, false);
		writer.start();
		final ByteBuilder row=new ByteBuilder();
		int count=0;
		try{
			for(Contig c : contigs){
				if(c.length()<min){continue;}
				writer.print(row.clear().append('>').append(c.name).nl().append(c.bases).nl()); count++;
			}
		}finally{if(writer.poisonAndWait()){throw new RuntimeException("FASTA output failed: "+path);}}
		return count;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Counts all supported competitors before choosing any partner. */
	static final class Choices {
		Choices(final int contigs){
			if(contigs<0 || contigs>Integer.MAX_VALUE/2){throw new IllegalArgumentException("Too many physical endpoints.");}
			degree=new int[2*contigs]; partner=new int[degree.length]; Arrays.fill(partner, -1);
		}
		void add(final int a, final int b){
			assert(a>=0 && a<b && b<degree.length) : "Candidates are ordered pairs of different physical ends.";
			degree[a]=Math.addExact(degree[a], 1); degree[b]=Math.addExact(degree[b], 1);
			partner[a]=b; partner[b]=a; supported++;
		}
		/** Frozen reciprocal unique choices; self-contig proposals remain unmerged. */
		int[] mates(){
			final int[] result=new int[degree.length]; Arrays.fill(result, -1);
			for(int a=0; a<degree.length; a++){
				final int b=partner[a];
				if(degree[a]==1 && b>=0 && degree[b]==1 && partner[b]==a && a/2!=b/2){result[a]=b;}
			}
			return result;
		}
		final int[] degree, partner;
		int supported;
	}

	/** Persistent spelled products and their original-unitig paths; no input is mutated. */
	static final class Result {
		Contig add(final byte[] bases, final String path){
			assert(bases!=null && path!=null) : "Every output component retains its sequence and provenance.";
			final int id=contigs.size();
			final Contig result=new Contig(bases, "samek_"+id, id);
			contigs.add(result); paths.add(path);
			return result;
		}
		final ArrayList<Contig> contigs=new ArrayList<Contig>();
		final ArrayList<String> paths=new ArrayList<String>();
		int joins, cycleUnitigs;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final PrintStream outstream;
	private String input, output, candidates, labels, baseline, paths, features, network, scores;
	private float cutoff=Float.NaN;
	private boolean cutoffSpecified;
	private int k, minContig=500;
	private static final String LABEL_HEADER="event\tk\tlabel\ttarget\twhole_hits\tleft_hits\tright_hits\tsource_hits\tsource_unknown\tsource_left_mask\tsource_right_mask\tsource_pair_mask\tdest_hits\tdest_unknown\tdest_left_mask\tdest_right_mask\tdest_pair_mask";
}
