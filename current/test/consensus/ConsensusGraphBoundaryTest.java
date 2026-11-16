package test.consensus;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import consensus.BaseGraph;
import consensus.BaseGraphHelper;
import consensus.BaseNode;
import consensus.ConsensusObject;
import stream.Read;

/** Finite graph output, score-range and trimming contracts. @author Shinobu */
public final class ConsensusGraphBoundaryTest{

	/** Uses phase=baseline to identify the published defects before repairs. */
	public static void main(final String[] args) throws Exception{
		boolean baseline=false;
		for(final String arg : args){
			if(arg.equals("phase=baseline")){baseline=true;}
			else if(!arg.equals("phase=fixed")){throw new IllegalArgumentException(arg);}
		}
		final boolean oldTrim=ConsensusObject.trimNs, oldQuality=BaseGraph.qualityByMS;
		final float oldFraction=ConsensusObject.trimDepthFraction;
		try{
			ConsensusObject.trimNs=false;
			ConsensusObject.trimDepthFraction=0;
			BaseGraph.qualityByMS=true;
			text(baseline);
			scores(baseline);
			trim(baseline);
			check(Read.VALIDATE_IN_CONSTRUCTOR && Read.CHANGE_QUALITY && Read.MAX_CALLED_QUALITY()==50,
				"The public-output control requires the ordinary constructor and quality cap.");
			// The vector validator conditionally caps only when it also converts bases.
			// Select the branchless validator explicitly for the capped-output control.
			final boolean oldVector=Read.VALIDATE_VECTOR;
			try{
				Read.VALIDATE_VECTOR=false;
				qualities(baseline, false);
			}finally{Read.VALIDATE_VECTOR=oldVector;}
			// Observe graph-produced bytes before Read's independent quality capping.
			final boolean oldValidation=Read.VALIDATE_IN_CONSTRUCTOR;
			try{
				Read.VALIDATE_IN_CONSTRUCTOR=false;
				qualities(baseline, true);
			}finally{Read.VALIDATE_IN_CONSTRUCTOR=oldValidation;}
			counts();
			if(!baseline){largeWeights(); orientations(); subclassSelection();}
			System.out.println("GRAPH_BOUNDARY_PASS phase="+(baseline ? "baseline" : "fixed")+" checks="+checks);
		}finally{
			ConsensusObject.trimNs=oldTrim;
			ConsensusObject.trimDepthFraction=oldFraction;
			BaseGraph.qualityByMS=oldQuality;
		}
	}

	/** Every reference position appears exactly once, without unsolicited stderr. */
	private static void text(final boolean baseline) throws Exception{
		final BaseGraph graph=graph("text", "AC");
		final PrintStream previous=System.err;
		final ByteArrayOutputStream bytes=new ByteArrayOutputStream();
		final String actual;
		try(PrintStream capture=new PrintStream(bytes, true, "UTF-8")){
			System.setErr(capture);
			actual=graph.toText().toString();
		}finally{System.setErr(previous);}
		final String first="(Node,0,REF,1,0,0,0 ,REF:1)";
		final String second="(Node,1,REF,0,1,0,0 )";
		final String expected="text:{"+(baseline ? first : "")+first+","+second+"}";
		check(expected.equals(actual), "Each of two reference positions must have one text entry: "+actual);
		check(new String(bytes.toByteArray(), StandardCharsets.UTF_8).trim().equals(baseline ? "..." : ""),
			"Graph text conversion must not write an unsolicited diagnostic.");
	}

	/** Uniform allele probabilities have no best-to-worst range to normalize. */
	private static void scores(final boolean baseline){
		final BaseGraph flat=graph("flat", "A");
		BaseGraphHelper.setCounts(flat, 0, 1, 1, 1, 1, 0);
		BaseGraphHelper.initForScoring(flat);
		final Read read=aligned("A");
		for(final boolean old : new boolean[]{false, true}){
			for(final boolean local : new boolean[]{false, true}){
				check((old ? flat.scoreOld(read, local, false) : flat.score(read, local, false))==0.25f,
					"Raw score remains meaningful even when all four alleles tie.");
				boolean rejected=false;
				float value=0;
				try{value=old ? flat.scoreOld(read, local, true) : flat.score(read, local, true);}
				catch(ArithmeticException expected){rejected=true;}
				check(baseline ? !rejected && Float.isNaN(value) : rejected,
					"A zero score range must be rejected explicitly instead of returning NaN.");
			}
		}
		final BaseGraph definite=graph("definite", "A");
		BaseGraphHelper.initForScoring(definite);
		check(definite.score(read, false, true)==1f && definite.scoreOld(read, false, true)==1f,
			"A unique matching allele retains normalized score one in both scorers.");
	}

	/** Checks explicit N removal while retaining the shared minimum-depth rule. */
	private static void trim(final boolean baseline){
		final BaseGraph graph=graph("trim", "NACN");
		ConsensusObject.trimNs=false;
		ConsensusObject.trimDepthFraction=0;
		check(sequence(graph.traverse()).equals("NACN"), "Untrimmed consensus preserves terminal Ns.");
		ConsensusObject.trimDepthFraction=0.5f;
		check(sequence(graph.traverse()).equals(baseline ? "AC" : "NACN"),
			"Depth one meets the depth threshold; trimNs=false must preserve terminal Ns.");
		ConsensusObject.trimNs=true;
		check(sequence(graph.traverse()).equals("AC"), "Explicit N trimming removes both terminal Ns.");
		ConsensusObject.trimDepthFraction=0;
		check(sequence(graph.traverse()).equals("AC"), "N trimming also works without depth trimming.");
		ConsensusObject.trimNs=false;
	}

	/** Internal Read qualities are numeric Phred; ASCII encoding belongs to writers. */
	private static void qualities(final boolean baseline, final boolean beforeCap){
		final BaseGraph graph=graph("quality", "A");
		BaseGraphHelper.setCounts(graph, 0, 9, 1, 0, 0, 0);
		BaseGraph.qualityByMS=true;
		check(graph.traverse().quality[0]==33, "Default diversity quality remains round(39*0.85)=33.");
		BaseGraph.qualityByMS=false;
		final int q=graph.traverse().quality[0];
		check(q==(baseline ? (beforeCap ? 63 : 50) : 9),
			"Nine-to-one allele odds yield numeric Phred nine, without an ASCII offset: "+q);
		final byte[] result=new byte[2];
		final BaseNode mixed=new BaseNode('A', ConsensusObject.REF, 0);
		mixed.add((byte)'A', 9);
		mixed.add((byte)'A', 9);
		mixed.add((byte)'C', 2);
		mixed.consensus(result);
		check(result[0]=='A' && result[1]==(baseline ? 35 : 9),
			"Eighteen supporting and two opposing weight units yield numeric quality nine.");
		final BaseNode unanimous=new BaseNode('A', ConsensusObject.REF, 0);
		unanimous.add((byte)'A', 20);
		unanimous.add((byte)'A', 20);
		unanimous.consensus(result);
		check(result[1]==(baseline ? 35 : 36),
			"Forty supporting weights use the existing 0.01 denominator floor: quality36.");
		final BaseGraph insertion=graph("insertion", "A");
		BaseGraphHelper.setCounts(insertion, 0, 20, 0, 0, 0, 0);
		final BaseNode inserted=new BaseNode('A', ConsensusObject.INS, 0);
		inserted.add((byte)'A', 29);
		inserted.add((byte)'A', 30);
		insertion.ref[0].insEdge=inserted;
		final float oldMaf=ConsensusObject.MAF_ins;
		try{
			ConsensusObject.MAF_ins=0.05f;
			final Read output=insertion.traverse();
			check(sequence(output).equals("AA"), "The two supported insertion observations pass the explicit frequency threshold.");
			check(output.quality[1]==(baseline ? 36 : 4),
				"Insertion confidence uses 59/20=2.95 before logarithm, yielding numeric quality four.");
		}finally{ConsensusObject.MAF_ins=oldMaf;}
		BaseGraph.qualityByMS=true;
	}

	/** Model import replaces counts and initializes probabilities after loading. */
	private static void counts(){
		final BaseGraph graph=graph("counts", "A");
		BaseGraphHelper.setCounts(graph, 0, 1, 2, 3, 4, 5);
		check(Arrays.equals(BaseGraphHelper.getCounts(graph, 0), new int[]{1, 2, 3, 4, 5}),
			"Model import/export must retain all four base counts and the deletion count.");
		check(BaseGraphHelper.length(graph)==1, "Model coordinates retain the reference length.");
		BaseGraphHelper.initForScoring(graph);
		check(Math.abs(graph.ref[0].baseProb((byte)'T')-0.4f)<0.000001f,
			"Allele probabilities normalize over the ten base observations, excluding deletions.");
	}

	/** Reference and deletion weights can each fit an int while their sum does not. */
	private static void largeWeights(){
		final BaseGraph graph=graph("large-weights", "AA");
		BaseGraphHelper.setCounts(graph, 0, 20000000, 0, 0, 0, 20000000);
		graph.ref[0].weightSum=2000000000;
		graph.ref[0].acgtWeight[0]=2000000000;
		graph.del[0].weightSum=2000000000;
		final BaseNode insertion=new BaseNode('A', ConsensusObject.INS, 1);
		insertion.countSum=20000000;
		insertion.acgtCount[0]=20000000;
		insertion.weightSum=300000000;
		insertion.acgtWeight[0]=300000000;
		graph.ref[0].insEdge=insertion;
		check(sequence(graph.traverse()).equals("AA"),
			"Insertion weight300M is below half of the4B reference/deletion sum; int overflow must not admit it.");
	}

	/** Checks ordinary exact orientations without claiming to resolve the historical screening TODO. */
	private static void orientations(){
		final String forward="AAAACCCCGGGGTTTTACGATCGTACCA";
		final String reverse="TGGTACGATCGTAAAACCCCGGGGTTTT";
		final BaseGraph graph=graph("orientation", forward);
		for(int strand=0; strand<2; strand++){
			final byte[] bases=(strand==0 ? forward : reverse).getBytes(StandardCharsets.US_ASCII);
			final byte[] original=bases.clone();
			final Read read=new Read(bases, null, "orientation", strand);
			check(graph.findBestOrientation(read, 0.8f)==strand,
				"An exact asymmetric reference or its literal reverse complement must choose the matching orientation.");
			check(Arrays.equals(bases, original), "Orientation screening restores the caller's query bases.");
		}
	}

	/** Internal base-only selection must retain public overrides for both node kinds. */
	private static void subclassSelection(){
		final boolean oldQuality=BaseGraph.qualityByMS;
		try{
			for(final boolean insertion : new boolean[]{false, true}){
				final BaseGraph graph=graph("override", "A");
				BaseGraphHelper.setCounts(graph, 0, 2, 0, 0, 0, 0);
				final OverrideNode node=new OverrideNode(insertion ? ConsensusObject.INS : ConsensusObject.REF);
				if(insertion){graph.ref[0].insEdge=node;}else{graph.ref[0]=node;}
				int calls=0;
				for(final boolean byMS : new boolean[]{true, false}){
					BaseGraph.qualityByMS=byMS;
					final Read output=graph.traverse();
					check(sequence(output).equals(insertion ? "AC" : "C"),
						"Traversal must retain a subclass's public allele-selection override.");
					check(node.calls==++calls, "Each emitted custom node must invoke its public override once.");
					check(output.quality[output.length()-1]==(byMS ? 39 : 11),
						"Default mode uses diversity quality; alternate mode retains the custom numeric confidence.");
				}
			}
		}finally{BaseGraph.qualityByMS=oldQuality;}
	}

	/** Deliberately overrides public selection so internal shortcuts cannot bypass clients. */
	private static final class OverrideNode extends BaseNode{
		OverrideNode(final int type){
			super(type==ConsensusObject.REF ? 'A' : 'C', type, 0);
			final byte observed=(byte)(type==ConsensusObject.REF ? 'A' : 'C');
			add(observed, 20);
			add(observed, 20);
		}
		@Override
		public byte consensus(final byte[] result){
			calls++;
			result[0]='C';
			result[1]=11;
			return result[0];
		}
		int calls;
	}

	private static BaseGraph graph(final String name, final String bases){
		return new BaseGraph(name, bases.getBytes(StandardCharsets.US_ASCII), null, 0, 0);
	}

	private static Read aligned(final String bases){
		final Read read=new Read(bases.getBytes(StandardCharsets.US_ASCII), null, "read", 0);
		read.match=new byte[bases.length()];
		Arrays.fill(read.match, (byte)'m');
		read.start=0;
		read.stop=bases.length()-1;
		read.setMapped(true);
		read.setShortMatch(false);
		return read;
	}

	private static String sequence(final Read read){return new String(read.bases, StandardCharsets.US_ASCII);}
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
		checks++;
	}
	private static int checks;
}
