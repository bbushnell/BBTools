package assemble;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectIntMap;
import parse.LineParser1;
import structures.ByteBuilder;

/** Independent reference/trace fixtures for strict live-node labels and deduplication. @author Fischl */
public final class LiveBranchCorpusTest {

	/** Generates isolated test files or checks the exact expected label/feature tables. */
	public static void main(final String[] args){
		if(args.length==3 && args[0].equals("generate_sampling")){generateSampling(args[1], args[2]);}
		else if(args.length==3 && args[0].equals("check_sampling")){checkSampling(args[1], Integer.parseInt(args[2]));}
		else if(args.length==3 && args[0].equals("generate")){generate(args[1], args[2]);}
		else if(args.length==2 && args[0].equals("check")){checkOutput(args[1]);}
		else{throw new IllegalArgumentException("generate NEW_directory variant, or check output_directory");}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Nine grey nodes, twenty known controls, two unknown controls and one excluded node. */
	private static void generateSampling(final String root, final String variant){
		if(!new File(root).mkdir()){throw new IllegalArgumentException("Sampling fixture directory must be new.");}
		final String[] words=words(32);
		final ByteBuilder ref=new ByteBuilder();
		final ArrayList<String> lines=new ArrayList<String>();
		for(int i=0; i<words.length; i++){
			lines.add(row(words[i], i<9 ? 0 : i==31 ? 8 : 5, false));
			if(i==30){continue;}
			final boolean unsafe=i>=9 && i<29 && ((i%2==0)^variant.equals("flipped"));
			ref.append('>').append(i).nl().append(i==29 ? 'N' : unsafe ? 'G' : 'A').append(words[i]).append('C').nl();
		}
		if(variant.equals("reverse")){Collections.reverse(lines);}
		if(variant.equals("duplicate")){lines.add(row(words[0], 0, false)); lines.add(row(words[9], 5, false));}
		final ByteBuilder trace=new ByteBuilder("#livebranch-observation-v2\n").append(LiveBranchObserver.HEADER).nl();
		for(String line : lines){trace.append(line).nl();}
		trace.append("#complete\t").append(lines.size()).append("\t1000\t10000\t100000\t1\t24\t2\t0.5\n");
		write(root+"/trace.tsv", trace); write(root+"/reference.fa", ref);
	}

	/** Checks identities against the full audit and prints selected control words for cross-run comparisons. */
	private static void checkSampling(final String root, final int expectedControls){
		check(LiveBranchCorpus.controlLimit(233, 154759, 0.1)==25, "Control fraction must bound the emitted share, not the input population.");
		check(LiveBranchCorpus.controlLimit(366, 175464, 0.1)==40, "Middle pilot control limit differs.");
		check(LiveBranchCorpus.controlLimit(563, 184909, 0.1)==62, "Final pilot control limit differs.");
		check(LiveBranchCorpus.controlLimit(0, 20, 0.1)==0 && LiveBranchCorpus.controlLimit(9, 20, 0)==0,
				"Zero grey population or zero requested fraction must emit no sampled controls.");
		check(LiveBranchCorpus.controlLimit(9, 20, 1)==20 && LiveBranchCorpus.controlLimit(90, 2, 0.1)==2,
				"Default all-controls and scarce-control behavior changed.");
		final String[] words=words(32);
		final ObjectIntMap<String> ids=new ObjectIntMap<String>(String.class);
		for(int i=0; i<words.length; i++){ids.put(words[i], i);}
		final int[] targets=new int[words.length];
		final boolean[] seen=new boolean[words.length];
		final LineParser1 fields=new LineParser1('\t');
		final ByteFile audit=ByteFile.makeByteFile(root+"/nodes.tsv", true);
		int audited=0;
		try{
			check(audit.nextLine()!=null, "Missing audit header.");
			for(byte[] line=audit.nextLine(); line!=null; line=audit.nextLine()){
				fields.set(line);
				final int id=ids.get(fields.parseString(0));
				check(id>=0 && !seen[id], "Sampling changed the complete unique-node audit.");
				seen[id]=true; audited++;
				targets[id]=fields.termEquals(".", 5) ? -1 : fields.parseInt(5);
			}
		}finally{if(audit.close()){throw new IllegalStateException("Cannot read sampling audit.");}}
		check(audited==32, "Sampling must retain all nodes in the audit.");
		java.util.Arrays.fill(seen, false);
		final ByteFile rows=ByteFile.makeByteFile(root+"/rows.tsv", true);
		final ByteFile vectors=ByteFile.makeByteFile(root+"/vectors.tsv", true);
		int grey=0, controls=0, count=0;
		try{
			check(rows.nextLine()!=null && vectors.nextLine()!=null && vectors.nextLine()!=null, "Missing sampled output header.");
			for(byte[] line=rows.nextLine(); line!=null; line=rows.nextLine()){
				fields.set(line);
				final String word=fields.parseString(1);
				final int id=ids.get(word), target=fields.parseInt(3);
				check(id>=0 && !seen[id] && fields.parseInt(0)==++count, "Duplicate, invented or nonconsecutive sampled row.");
				seen[id]=true;
				check(target>=0 && target==targets[id], "Sampling altered truth or included an unresolved node.");
				if(id<9){check(fields.termEquals("grey", 2), "Grey node changed role."); grey++;}
				else{
					check(id<29 && fields.termEquals("control", 2), "Ineligible control entered the sample.");
					controls++; System.out.println(word);
				}
				final byte[] vector=vectors.nextLine();
				check(vector!=null, "Missing paired training vector."); fields.set(vector);
				check(fields.terms()==13 && fields.parseInt(12)==target, "Sampled vector and provenance target differ.");
			}
			check(vectors.nextLine()==null, "Extra sampled vector.");
		}finally{if(rows.close() | vectors.close()){throw new IllegalStateException("Cannot read sampled output.");}}
		check(grey==9 && controls==expectedControls, "Sampling dropped grey nodes or used the wrong control budget.");
		System.err.println("LIVE_BRANCH_SAMPLING_TEST_PASS grey="+grey+" controls="+controls);
	}

	/** Ten distinct nodes cover true/false/unknown routes and control eligibility. */
	private static void generate(final String root, final String variant){
		if(!new File(root).mkdir()){throw new IllegalArgumentException("Fixture directory must be new.");}
		final String[] words=words();
		final ByteBuilder ref=new ByteBuilder();
		for(int i=0; i<words.length; i++){
			if(i==4){continue;}// A genuinely absent query stays unresolved.
			ref.append('>').append(i).nl();
			if(i==2){ref.append('A').append(words[i]).append("CNNNA").append(words[i]).append("CNNNG").append(words[i]).append('C');}
			else{ref.append(i==1 ? 'G' : i==3 ? 'N' : 'A').append(words[i]).append('C');}
			ref.nl();
		}
		write(root+"/reference.fa", ref);
		final ArrayList<String> lines=new ArrayList<String>();
		for(int i=0; i<words.length; i++){lines.add(row(words[i], i, i==9));}
		lines.add(row(words[0], 0, false));
		lines.add(row(words[0], 0, true));
		if(variant.equals("duplicate")){lines.add(row(words[0], 0, false));}
		if(variant.equals("reverse")){Collections.reverse(lines);}
		if(variant.equals("conflict")){lines.set(10, lines.get(10).replace("\t24\t10\t", "\t24\t11\t"));}
		if(variant.equals("bad_word")){lines.set(0, "N"+lines.get(0).substring(1));}
		if(variant.equals("bad_width")){lines.set(0, lines.get(0)+"\t0");}
		if(variant.equals("bad_entry")){lines.set(11, row(words[0], 0, true).replace("\t-1\t-1\t1\t", "\t0\t0\t1\t"));}
		final ByteBuilder trace=new ByteBuilder(variant.equals("old_schema") ? "#livebranch-observation-v1\n" : "#livebranch-observation-v2\n");
		trace.append(LiveBranchObserver.HEADER).nl();
		for(String line : lines){trace.append(line).nl();}
		if(!variant.equals("missing_trailer")){
			trace.append("#complete\t").append(lines.size()+(variant.equals("bad_rows") ? 1 : 0));
			trace.append("\t100\t").append(variant.equals("bad_moments") ? "NaN" : "1.0E3");
			trace.append("\t1.0E4\t1\t24\t").append(variant.equals("bad_mce") ? 11 : 2).append("\t0.5\n");
		}
		if(variant.equals("after_trailer")){trace.append(lines.get(0)).nl();}
		write(root+"/trace.tsv", trace);
	}

	/** Literal evidence, independent of the production collector and policy implementation. */
	private static String row(final String word, final int index, final boolean entry){
		assert(word.length()==24) : "This fixture has one fixed K, unrelated to packed short/long representation.";
		final int second=index==6 ? 5 : index==7 ? 10 : index==9 ? 0 : 2;
		final boolean branch=index<5 || index==6 || index==7;
		final int level=branch ? (index==6 || index==7 ? 1 : 2) : 0;
		final ByteBuilder row=new ByteBuilder(word).append("\t24\t10\t10\t0\t").append(second);
		row.append("\t0\t0\t10\t0\t0\t").append(branch ? 1 : 0).append("\t0\t");
		row.append(entry ? -1 : 0).tab().append(entry ? -1 : index==8 ? 2 : 0);
		row.tab().append(entry ? 1 : 0).tab().append(index==8 ? 1 : 0).tab().append(level).append("\t1");
		return row.toString();
	}

	/** Requires labels from all reference occurrences, not just the first two repeat hits. */
	private static void checkOutput(final String root){
		final String[] words=words();
		final ObjectIntMap<String> indices=new ObjectIntMap<String>(String.class);
		for(int i=0; i<words.length; i++){indices.put(words[i], i);}
		final boolean[] seen=new boolean[words.length];
		final ByteFile input=ByteFile.makeByteFile(root+"/nodes.tsv", true);
		final LineParser1 fields=new LineParser1('\t');
		int nodes=0;
		try{
			check(input.nextLine()!=null, "Missing node header.");
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				fields.set(line);
				check(fields.terms()==26, "Audit output must contain14 provenance fields and12 inputs.");
				final int id=indices.get(fields.parseString(0));
				check(id>=0 && !seen[id], "A canonical node was lost, duplicated or invented.");
				seen[id]=true; nodes++;
				final String label=id==1 || id==2 ? "unsafe" : id==3 ? "boundary" : id==4 ? "absent" : id==7 ? "no_route" : "safe";
				final String role=id<5 ? "grey" : id==5 ? "control" : id==6 || id==7 ? "hard" : "excluded";
				check(fields.termEquals(role, 3) && fields.termEquals(label, 4), "Wrong label or control role at fixture"+id);
				check(fields.termEquals(id==3 || id==4 || id==7 ? "." : id==1 || id==2 ? "0" : "1", 5), "Unknown truth became a binary label.");
				if(id==2){check(fields.parseLong(6)==3 && fields.parseInt(8)==5, "Third repeat occurrence/conflicting left flank was ignored.");}
				if(id==3){check(fields.parseLong(7)==1, "N flank lost its uncertainty marker.");}
				check(fields.parseInt(1)==24, "Output K changed.");
				// Shared producer's fixed transforms are checked here independently.
				final double log10=0.125*Math.log(11)/Math.log(2);
				check(Math.abs(fields.parseDouble(14)-log10)<1e-6, "Node-depth normalization differs.");
				check(Math.abs(fields.parseDouble(23)-0.125*Math.log(24)/Math.log(2))<1e-6, "K normalization differs.");
				check(Math.abs(fields.parseDouble(24)-log10)<1e-6 && Math.abs(fields.parseDouble(25)-log10)<1e-6, "Global depth features differ.");
			}
		}finally{if(input.close()){throw new IllegalStateException("Cannot read labeled audit.");}}
		check(nodes==10, "Expected exactly ten unique junctions.");
		final ByteFile vectors=ByteFile.makeByteFile(root+"/vectors.tsv", true);
		int rows=0, positives=0;
		try{
			check("#dims\t12\t1".equals(text(vectors.nextLine())), "Wrong trainer dimensions.");
			check(vectors.nextLine()!=null, "Missing trainer schema.");
			for(byte[] line=vectors.nextLine(); line!=null; line=vectors.nextLine()){
				fields.set(line);
				check(fields.terms()==13, "Wrong trainer width.");
				for(int i=0; i<12; i++){check(Float.isFinite(fields.parseFloat(i)), "Nonfinite trainer feature.");}
				final int target=fields.parseInt(12);
				check(target==0 || target==1, "Invalid trainer target.");
				rows++; positives+=target;
			}
		}finally{if(vectors.close()){throw new IllegalStateException("Cannot read trainer output.");}}
		check(rows==4 && positives==2, "Hard, unresolved, hidden-only or entry-only nodes leaked into training, or duplicates changed weights.");
		checkCensus(root);
		System.err.println("LIVE_BRANCH_CORPUS_TEST_PASS nodes=10 training=4 positive=2 negative=2");
	}

	/** Literal role totals distinguish the candidate population from hard stops and sampled controls. */
	private static void checkCensus(final String root){
		final String[] expected={"role\tsafe\tunsafe\tunknown", "grey\t1\t2\t2",
				"hard\t1\t0\t1", "control\t1\t0\t0", "excluded\t2\t0\t0"};
		final ByteFile file=ByteFile.makeByteFile(root+"/census.tsv", true);
		try{
			for(String row : expected){check(row.equals(text(file.nextLine())), "Role census differs from the independent ten-node truth fixture.");}
			check(file.nextLine()==null, "Unexpected extra census rows.");
		}finally{if(file.close()){throw new IllegalStateException("Cannot read role census.");}}
	}

	/** Deterministic distinct canonical words; no source genome or learned model needed. */
	private static String[] words(){
		return words(10);
	}

	/** Shared generator retains the original ten-node fixture prefix. */
	private static String[] words(final int size){
		final Random random=new Random(381502);
		final String[] answer=new String[size];
		for(int i=0; i<answer.length; i++){
			final StringBuilder word=new StringBuilder(), reverse=new StringBuilder();
			for(int j=0; j<24; j++){word.append("ACGT".charAt(random.nextInt(4)));}
			for(int j=23; j>=0; j--){reverse.append("TGCA".charAt("ACGT".indexOf(word.charAt(j))));}
			final String a=word.toString(), b=reverse.toString();
			answer[i]=a.compareTo(b)<=0 ? a : b;
			for(int j=0; j<i; j++){check(!answer[i].equals(answer[j]), "Fixture must have distinct exact words.");}
		}
		return answer;
	}

	/** Native writer with explicit failure propagation for test fixtures. */
	private static void write(final String path, final ByteBuilder bytes){
		assert(bytes.length()>0) : "A fixture must not silently publish an empty input.";
		final ByteStreamWriter writer=new ByteStreamWriter(path, false, false, false);
		writer.start(); writer.print(bytes);
		if(writer.poisonAndWait()){throw new IllegalStateException("Failed to write fixture: "+path);}
	}

	/** Converts only small header records. */
	private static String text(final byte[] bytes){return bytes==null ? null : new String(bytes, StandardCharsets.US_ASCII);}

	/** Test checks remain loud independently of the assertion launch flag. */
	private static void check(final boolean condition, final String message){if(!condition){throw new AssertionError(message);}}
}
