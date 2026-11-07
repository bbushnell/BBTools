package test.prok;

import java.io.File;
import java.io.PrintStream;

import dna.AminoAcid;
import map.ObjectIntMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaSupport.Graph;
import prok.NcrnaSupportIO;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.ReduceNcrnaModels;
import shared.Shared;
import structures.ByteBuilder;

/** Native CLI fixtures for physical-window versus oriented-query scan semantics.
 * Uses the same Parser/PreParser and native byte I/O idioms as the sequence template.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaSupportScanTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaSupportScanTest.class, false);
		try{
			final Parser parser=new Parser();
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(java.util.Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown test option: "+arg);}
			}
			if(parser.out1==null){throw new IllegalArgumentException("out= must name a new fixture directory");}
			final File root=new File(parser.out1);
			if(!root.mkdir()){throw new IllegalArgumentException("Fixture directory must be new: "+root);}
			run(root, pp.outstream);
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Independent literal sequences establish three physical matches, including a reverse-only match. */
	private static void run(final File root, final PrintStream log){
		Shared.setBufferLen(1);
		final String ref=path(root, "ref.fa");
		final String queries=path(root, "windows.fa");
		final String negative=path(root, "negative.fa");
		final String hbm=path(root, "models.hbm");
		write(ref, ">m0\n"+GENE+"\n>m1\n"+PALINDROME+"\n");
		write(queries, ">forward\n"+GENE+"\n>reverse/+\n"+AminoAcid.reverseComplementBases(GENE)
			+"\n>palindrome\n"+PALINDROME+"\n>negative\n"+NEGATIVE+"\n");
		write(negative, ">negative\n"+NEGATIVE+"\n");
		writeHbm(hbm, false);
		final File both=scan(root, "both", queries, ref, "hbm="+hbm, "t=1", "reportcuts=1");
		final Graph graph=NcrnaSupportIO.read(path(both, "original.support.tsv"));
		assert(!graph.reduced && graph.queryNames.length==8 && graph.modelNames.length==2)
			: "Four physical windows require eight genuine oriented queries against both models";
		checkSummary(both, 4, 8, 16, 3);
		checkMapping(both, graph, 2);
		final ObjectIntMap<String> indices=index(graph);
		assert(graph.bestIdentities[indices.get("forward/+")]==1)
			: "The literal forward gene is identical to model m0";
		assert(graph.bestIdentities[indices.get("reverse/+/+")]<1
			&& graph.bestIdentities[indices.get("reverse/+/-")]==1)
			: "A source ID already ending in /+ stays distinct, and its reverse orientation restores the exact gene";
		assert(graph.bestIdentities[indices.get("palindrome/+")]==1
			&& graph.bestIdentities[indices.get("palindrome/-")]==1)
			: "The palindromic match passes twice but contributes only one physical-window acceptance";
		final int n=indices.get("negative/+");
		final int nr=indices.get("negative/-");
		final float negativeBest=Math.max(graph.bestIdentities[n], graph.bestIdentities[nr]);
		assert(negativeBest>=0 && negativeBest<1 && graph.starts[n]==graph.starts[n+1]
			&& graph.starts[nr]==graph.starts[nr+1])
			: "Unsupported negative orientations retain finite best identities and no exact-match edges";

		final File multi=scan(root, "multi", queries, ref, "t=3", "reportcuts=1");
		compare(graph, NcrnaSupportIO.read(path(multi, "original.support.tsv")));
		checkSummary(multi, 4, 8, 16, 3);
		checkMapping(multi, NcrnaSupportIO.read(path(multi, "original.support.tsv")), 2);
		final File forward=scan(root, "forward", queries, ref, "rcomp=f", "reportcuts=1");
		checkSummary(forward, 4, 4, 8, 2);
		checkMapping(forward, NcrnaSupportIO.read(path(forward, "original.support.tsv")), 1);
		final File allNegative=scan(root, "all_negative", negative, ref, "reportcuts=1");
		checkSummary(allNegative, 1, 2, 4, 0);

		// Native float equality is inclusive; the immediately higher float excludes this window.
		final File boundary=scan(root, "boundary", queries, ref,
			"reportcuts="+Float.toString(negativeBest)+","+Float.toString(Math.nextUp(negativeBest)));
		checkSummary(boundary, 4, 8, 16, 4, 3);
		for(final String conflict:new String[]{"reduce=t", "reduce=f", "regenerate=f", "callerverify=t",
				"target=1", "graph=missing.tsv", "genomes=missing.txt", "reportcuts=NaN", "reportcuts=0.7,,1"}){
			reject(root, "reject_"+rejected++, queries, ref, conflict);
		}
		final String badHbm=path(root, "bad.hbm");
		writeHbm(badHbm, true);
		reject(root, "bad_pair", queries, ref, "hbm="+badHbm);
		final String duplicate=path(root, "duplicate.fa");
		write(duplicate, ">same\n"+GENE+"\n>same\n"+GENE+"\n");
		reject(root, "duplicate", duplicate, ref, "reportcuts=1");
		normalReduction(root, ref);
		log.println("NCRNA_SUPPORT_SCAN_TEST_PASS orientations mapping negative denominator exact_float pairing normal_reduction");
	}

	/** Invokes the actual public parser/driver; unspecified orientation must mean both only here. */
	private static File scan(final File root, final String name, final String in, final String ref, final String... extra){
		final File out=new File(root, name);
		final String[] args=new String[6+extra.length];
		args[0]="mode=scan";
		args[1]="in="+in;
		args[2]="ref="+ref;
		args[3]="out="+out;
		args[4]="cutoff=1";
		args[5]="t=1";
		System.arraycopy(extra, 0, args, 6, extra.length);
		ReduceNcrnaModels.main(args);
		assert(out.isDirectory() && out.list().length==3)
			: "Scan mode emits only graph, window mapping and summary; no library, reduction or caller artifacts";
		return out;
	}

	/** Counts must describe physical windows even when two orientations pass. */
	private static void checkSummary(final File directory, final int windows, final int queries,
			final long alignments, final long... passing){
		final LineParser1 lp=new LineParser1('\t');
		try(Input input=new Input(path(directory, "scan.summary.tsv"))){
			lp.set(input.reader.nextLine());
			assert(lp.terms()==11 && lp.termEquals("windows", 2) && lp.termEquals("passing_windows", 5))
				: "The scan summary must expose its physical and oriented denominators";
			for(int i=0; i<passing.length; i++){
				lp.set(input.reader.nextLine());
				assert(lp.parseInt(2)==windows && lp.parseInt(3)==queries && lp.parseInt(4)==2
					&& lp.parseLong(5)==passing[i] && lp.parseLong(7)==alignments)
					: "Literal fixture denominators/work/acceptance differ at threshold row "+i;
			}
			assert(input.reader.nextLine()==null) : "Exactly one summary row is required per requested report threshold";
		}
	}

	/** Every graph query maps once to a real source, including suffix-bearing source IDs. */
	private static void checkMapping(final File directory, final Graph graph, final int orientations){
		final ObjectIntMap<String> ids=index(graph);
		final ObjectIntMap<String> windows=new ObjectIntMap<String>(String.class);
		final boolean[] seen=new boolean[graph.queryNames.length];
		final LineParser1 lp=new LineParser1('\t');
		int rows=0;
		try(Input input=new Input(path(directory, "window_queries.tsv"))){
			input.reader.nextLine();
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				lp.set(line);
				final String window=lp.parseString(0);
				final String query=lp.parseString(1);
				final int index=ids.get(query);
				assert(lp.terms()==5 && index>=0 && !seen[index]) : "Mapping must neither invent nor duplicate an oriented query";
				seen[index]=true;
				final boolean plus=lp.termEquals('+', 2);
				assert(plus || lp.termEquals('-', 2)) : "A mapped orientation must be explicit";
				assert(query.equals(orientations==1 ? window : window+(plus ? "/+" : "/-")))
					: "Orientation suffixes must preserve the complete original identifier";
				final int count=windows.get(window);
				windows.put(window, count<0 ? 1 : count+1);
				rows++;
			}
		}
		assert(rows==graph.queryNames.length) : "Every oriented graph query requires exactly one mapping row";
		for(final String window:new String[]{"forward", "reverse/+", "palindrome", "negative"}){
			assert(windows.get(window)==orientations) : "Expected orientation count for original physical window "+window;
		}
	}

	private static ObjectIntMap<String> index(final Graph graph){
		final ObjectIntMap<String> map=new ObjectIntMap<String>(String.class);
		for(int i=0; i<graph.queryNames.length; i++){map.put(graph.queryNames[i], i);}
		assert(map.size()==graph.queryNames.length) : "Graph query names must remain unique for the result join";
		return map;
	}

	/** Worker order may differ; compare native scores and supporters by the named axis. */
	private static void compare(final Graph first, final Graph second){
		final ObjectIntMap<String> ids=index(first);
		assert(first.queryNames.length==second.queryNames.length) : "Worker count cannot change the oriented population";
		for(int i=0; i<second.queryNames.length; i++){
			final int j=ids.get(second.queryNames[i]);
			assert(j>=0 && first.bestModels[j]==second.bestModels[i] && first.bestIdentities[j]==second.bestIdentities[i])
				: "Native best identity/model must not depend on worker assignment";
			final int count=first.starts[j+1]-first.starts[j];
			assert(count==second.starts[i+1]-second.starts[i]) : "Worker count cannot add or lose support edges";
			for(int edge=0; edge<count; edge++){
				final int a=first.starts[j]+edge, b=second.starts[i]+edge;
				assert(first.models[a]==second.models[b] && first.identities[a]==second.identities[b])
					: "Named support edges must preserve the original native identity";
			}
		}
	}

	private static void reject(final File root, final String name, final String in, final String ref, final String option){
		boolean failed=false;
		try{scan(root, name, in, ref, option);}catch(IllegalArgumentException expected){failed=true;}
		assert(failed && !new File(root, name).exists()) : "Invalid input/mode must fail before output creation: "+option;
	}

	private static void normalReduction(final File root, final String ref){
		final String genes=path(root, "positive.fa");
		write(genes, ">positive\n"+GENE+"\n");
		final File out=new File(root, "ordinary");
		ReduceNcrnaModels.main(new String[]{"in="+genes, "ref="+ref, "cutoff=1", "target=1", "out="+out, "t=1"});
		final Graph graph=NcrnaSupportIO.read(path(out, "original.support.tsv"));
		assert(graph.queryNames.length==1 && graph.queryNames[0].equals("positive") && new File(out, "consensus.fa").isFile()
			&& !new File(out, "scan.summary.tsv").exists()) : "Ordinary reduction keeps original names and library output";
	}

	/** Literal dense paired metadata exercises validation without invoking a model builder. */
	private static void writeHbm(final String path, final boolean wrongName){
		final ByteBuilder bb=new ByteBuilder("#HBM\t1\n");
		final String[] bases={GENE, PALINDROME};
		for(int model=0; model<2; model++){
			bb.append('>').append(wrongName && model==0 ? "wrong" : "m"+model).nl().append(bases[model]).nl();
			for(int i=0; i<bases[model].length(); i++){bb.append("1\t1\t1\t1\t0\n");}
		}
		try(Output output=new Output(path)){output.writer.print(bb);}
	}

	private static String path(final File root, final String name){return new File(root, name).getPath();}
	private static void write(final String path, final String text){
		try(Output output=new Output(path)){output.writer.print(text);}
	}

	private static int rejected=0;
	private static final String GENE="ACGTGCTAGTCAGATCGATGCTACGATCGTAGCATGCTAGTCAGCTAGCATGCTACGATCGTAC";
	private static final String PALINDROME="ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
	private static final String NEGATIVE="TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
}
