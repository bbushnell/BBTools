package test.clump;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;

import clump.KmerSort2;
import stream.FASTQ;

/** Direct KmerSort2 paired-output oracle for explicit out2 and out=# expansion.
 * @author Jean
 */
public final class KmerSort2Out2Test{

	/** Runs both paired-output spellings in fresh subdirectories. */
	public static void main(String[] args) throws Exception{
		check(args.length==2, "Expected phase and a fresh fixture directory");
		final boolean fixed=args[0].equals("fixed");
		check(fixed || args[0].equals("baseline"), "Unknown phase: "+args[0]);
		final Path root=Paths.get(args[1]);
		Files.createDirectory(root);

		final boolean oldForce=FASTQ.FORCE_INTERLEAVED, oldTest=FASTQ.TEST_INTERLEAVED;
		try{
			System.out.println("case\tleft_exists\tright_exists\tliteral_exists"+
					"\tleft_reads\tright_reads\tliteral_reads"+
					"\tleft_ok\tright_ok\tpair_order_ok\tliteral_ok"+
					"\tbaseline_ok\tfixed_ok");
			final Result explicit=run(root.resolve("explicit"), false, true);
			final Result hash=run(root.resolve("hash"), true, false);
			System.out.println(explicit);
			System.out.println(hash);
			check(fixed ? explicit.fixedOk && hash.fixedOk : explicit.baselineOk && hash.baselineOk,
					"Unexpected "+args[0]+" KmerSort2 paired-output behavior");
		}finally{
			FASTQ.FORCE_INTERLEAVED=oldForce;
			FASTQ.TEST_INTERLEAVED=oldTest;
		}
	}

	/** Runs one direct KmerSort2 paired-output case. */
	private static Result run(Path root, boolean hash, boolean expectExplicitLoss)
			throws Exception{
		Files.createDirectory(root);
		write(root.resolve("in_0.fq"), GROUP0);
		write(root.resolve("in_1.fq"), GROUP1);

		Throwable thrown=null;
		try{
			KmerSort2.main(new String[] {
					"in="+root.resolve("in_%.fq"),
					hash ? "out="+root.resolve("out_#.fq") : "out="+root.resolve("left.fq"),
					hash ? "out2=null" : "out2="+root.resolve("right.fq"),
					"groups=2",
					"k=7",
					"seed=1",
					"border=1",
					"hashes=4",
					"ow=t",
					"qin=33",
					"qout=33",
					"interleaved=t",
					"ecco=f",
					"dedupe=f",
					"condense=f",
					"correct=f",
					"repair=f",
					"namesort=f",
					"verbose=f"
			});
		}catch(Throwable t){thrown=t;}

		final Path left=hash ? root.resolve("out_1.fq") : root.resolve("left.fq");
		final Path right=hash ? root.resolve("out_2.fq") : root.resolve("right.fq");
		final Path literal=root.resolve("out_#.fq");
		final LinkedHashMap<String, String> observedLeft=readIfPresent(left);
		final LinkedHashMap<String, String> observedRight=readIfPresent(right);
		final LinkedHashMap<String, String> observedLiteral=readIfPresent(literal);
		final boolean leftOk=observedLeft.equals(readText(LEFT));
		final boolean rightOk=observedRight.equals(readText(RIGHT));
		final boolean pairOrderOk=pairedOrderOk(observedLeft, observedRight);
		final boolean literalOk=observedLiteral.equals(readText(LEFT+RIGHT));
		final boolean literalExists=Files.exists(literal);
		final boolean fixedOk=thrown==null && leftOk && rightOk && pairOrderOk && !literalExists;
		final boolean observedExplicitLoss=thrown==null && !fixedOk
				&& observedLeft.equals(readText(LEFT+RIGHT))
				&& observedRight.isEmpty()
				&& !literalExists;
		final boolean observedLiteralHash=thrown==null && !fixedOk
				&& observedLeft.isEmpty()
				&& observedRight.isEmpty()
				&& literalOk;
		final boolean baselineOk=expectExplicitLoss ? observedExplicitLoss : observedLiteralHash;
		return new Result(root.getFileName().toString(), Files.exists(left), Files.exists(right),
				literalExists, observedLeft.size(), observedRight.size(), observedLiteral.size(),
				leftOk, rightOk, pairOrderOk, literalOk, baselineOk, fixedOk, thrown);
	}

	/** Checks that split R1/R2 outputs preserve corresponding mate order. */
	private static boolean pairedOrderOk(final LinkedHashMap<String, String> left,
			final LinkedHashMap<String, String> right){
		if(left.size()!=right.size()){return false;}
		final Iterator<String> leftIds=left.keySet().iterator();
		final Iterator<String> rightIds=right.keySet().iterator();
		while(leftIds.hasNext()){
			final String leftId=leftIds.next(), rightId=rightIds.next();
			if(!leftId.endsWith("/1") || !rightId.endsWith("/2")){return false;}
			if(!fragmentStem(leftId).equals(fragmentStem(rightId))){return false;}
		}
		return true;
	}

	/** Returns the fragment stem before its /1 or /2 mate suffix. */
	private static String fragmentStem(final String id){
		return id.substring(0, id.length()-2);
	}

	/** Writes one grouped input shard. */
	private static void write(Path path, String text) throws Exception{
		Files.write(path, text.getBytes(StandardCharsets.US_ASCII),
				StandardOpenOption.CREATE_NEW);
	}

	/** Returns parsed FASTQ records or an empty map for a missing output. */
	private static LinkedHashMap<String, String> readIfPresent(Path path) throws Exception{
		return Files.exists(path) ? read(Files.readAllLines(path, StandardCharsets.US_ASCII)) :
				new LinkedHashMap<String, String>();
	}

	/** Parses one FASTQ text fixture. */
	private static LinkedHashMap<String, String> readText(String text) throws Exception{
		return read(java.util.Arrays.asList(text.split("\n")));
	}

	/** Strict four-line FASTQ parser keyed by read ID. */
	private static LinkedHashMap<String, String> read(List<String> lines){
		check(lines.size()%4==0, "Malformed FASTQ line count: "+lines.size());
		final LinkedHashMap<String, String> map=new LinkedHashMap<String, String>();
		for(int i=0; i<lines.size(); i+=4){
			check(lines.get(i).startsWith("@"), "Malformed header: "+lines.get(i));
			check(lines.get(i+2).equals("+"), "Malformed separator: "+lines.get(i+2));
			check(lines.get(i+1).length()==lines.get(i+3).length(),
					"Sequence/quality length mismatch");
			final String prior=map.put(lines.get(i).substring(1),
					lines.get(i+1)+'\n'+lines.get(i+3));
			check(prior==null, "Duplicate output read "+lines.get(i));
		}
		return map;
	}

	/** Lightweight assertion independent of JVM assertion flags. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Recorded result for one paired-output case. */
	private static final class Result{

		Result(String label_, boolean leftExists_, boolean rightExists_,
				boolean literalExists_, int leftReads_, int rightReads_, int literalReads_,
				boolean leftOk_, boolean rightOk_, boolean pairOrderOk_,
				boolean literalOk_, boolean baselineOk_, boolean fixedOk_, Throwable thrown_){
			label=label_;
			leftExists=leftExists_;
			rightExists=rightExists_;
			literalExists=literalExists_;
			leftReads=leftReads_;
			rightReads=rightReads_;
			literalReads=literalReads_;
			leftOk=leftOk_;
			rightOk=rightOk_;
			pairOrderOk=pairOrderOk_;
			literalOk=literalOk_;
			baselineOk=baselineOk_ && thrown_==null;
			fixedOk=fixedOk_ && thrown_==null;
			thrown=thrown_;
		}

		@Override
		public String toString(){
			return label+'\t'+leftExists+'\t'+rightExists+'\t'+literalExists+
					'\t'+leftReads+'\t'+rightReads+'\t'+literalReads+
					'\t'+leftOk+'\t'+rightOk+'\t'+pairOrderOk+
					'\t'+literalOk+'\t'+baselineOk+'\t'+fixedOk+
					(thrown==null ? "" : "\t"+thrown.getClass().getName());
		}

		final String label;
		final boolean leftExists;
		final boolean rightExists;
		final boolean literalExists;
		final int leftReads;
		final int rightReads;
		final int literalReads;
		final boolean leftOk;
		final boolean rightOk;
		final boolean pairOrderOk;
		final boolean literalOk;
		final boolean baselineOk;
		final boolean fixedOk;
		final Throwable thrown;
	}

	private static final String LEFT="@beta/1\nACGTACGTACGT\n+\nIIIIIIIIIIII\n"+
			"@alpha/1\nTTGCAACCGGTT\n+\nHHHHHHHHHHHH\n"+
			"@delta/1\nGATCGATCGATC\n+\nGGGGGGGGGGGG\n"+
			"@gamma/1\nCCCCAAAATTTT\n+\nFFFFFFFFFFFF\n";

	private static final String RIGHT="@beta/2\nTTTTCCCCAAAA\n+\nIIIIIIIIIIII\n"+
			"@alpha/2\nAAAACCCCTTTT\n+\nHHHHHHHHHHHH\n"+
			"@delta/2\nTTAACCGGTTAA\n+\nGGGGGGGGGGGG\n"+
			"@gamma/2\nGGGGTTTTAAAA\n+\nFFFFFFFFFFFF\n";

	private static final String GROUP0="@beta/1\nACGTACGTACGT\n+\nIIIIIIIIIIII\n"+
			"@beta/2\nTTTTCCCCAAAA\n+\nIIIIIIIIIIII\n"+
			"@alpha/1\nTTGCAACCGGTT\n+\nHHHHHHHHHHHH\n"+
			"@alpha/2\nAAAACCCCTTTT\n+\nHHHHHHHHHHHH\n";

	private static final String GROUP1="@delta/1\nGATCGATCGATC\n+\nGGGGGGGGGGGG\n"+
			"@delta/2\nTTAACCGGTTAA\n+\nGGGGGGGGGGGG\n"+
			"@gamma/1\nCCCCAAAATTTT\n+\nFFFFFFFFFFFF\n"+
			"@gamma/2\nGGGGTTTTAAAA\n+\nFFFFFFFFFFFF\n";

}
