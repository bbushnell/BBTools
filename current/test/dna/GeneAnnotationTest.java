package test.dna;

import dna.Exon;
import dna.Gene;
import shared.Shared;

/** Literal inclusive-exon and source-preserving merge checks for legacy Gene records.
 * Runs without gene resource files or a mapping index, in either assertion mode.
 * @author Shinobu
 */
public final class GeneAnnotationTest{

	/** Runs case=controls|utr|merge-source in a fresh JVM. */
	public static void main(String[] args){
		require(args.length==1 && args[0].startsWith("case="), "Require one case= argument");
		final String mode=args[0].substring(5);
		if(mode.equals("controls")){controls();}
		else if(mode.equals("utr")){wholeExons();}
		else if(mode.equals("merge-source")){missingSource();}
		else{throw new IllegalArgumentException(mode);}
		System.out.println("PASS "+mode);
	}

	/** Preserves spanning/coding/untranslated lengths and ordinary source retention. */
	private static void controls(){
		require(lengths("spanning-plus", gene(Shared.PLUS, new int[][]{{0, 27}}, 10, 15, "seqGene", false),
				28, 6, 10, 12), "Spanning plus-strand lengths changed");
		require(lengths("spanning-minus", gene(Shared.MINUS, new int[][]{{0, 27}}, 10, 15, "seqGene", false),
				28, 6, 12, 10), "Spanning minus-strand lengths changed");
		require(lengths("coding-only", gene(Shared.PLUS, new int[][]{{10, 15}}, 10, 15, "seqGene", false),
				6, 6, 0, 0), "Fully coding exon acquired UTR bases");
		require(lengths("zero-coding-span", gene(Shared.PLUS, new int[][]{{0, 27}}, 10, 10, "seqGene", false),
				28, 0, 10, 17), "Equal coding endpoints must retain the separate legacy zero-coding convention");
		for(byte strand : new byte[]{Shared.PLUS, Shared.MINUS}){
			require(lengths("untranslated-"+strand, gene(strand, WHOLE, 10, 10, "seqGene", true),
					19, 0, 0, 0), "Untranslated records must retain zero reported UTR lengths");
		}
		final Gene receiver=gene(Shared.PLUS, new int[][]{{0, 27}}, 10, 15, "seqGene", false);
		final Gene donor=gene(Shared.PLUS, new int[][]{{0, 27}}, 10, 15, "knownGene", false);
		final Gene merged=receiver.merge(donor);
		require(merged.primarySource==0 && sourceColumn(merged).equals("seqGene"),
				"Merge must preserve the receiver's seqGene source rather than adopt the donor's source");
		require(merged!=receiver && merged.exons==receiver.exons && receiver.primarySource==0,
				"Merge must create a record borrowing the receiver's unchanged exon array");
	}

	/** Unequal sides expose the strand swap; single-base exons expose complete loss. */
	private static void wholeExons(){
		boolean ok=true;
		ok&=lengths("whole-plus", gene(Shared.PLUS, WHOLE, 10, 15, "seqGene", false), 19, 6, 5, 8);
		ok&=lengths("whole-minus", gene(Shared.MINUS, WHOLE, 10, 15, "seqGene", false), 19, 6, 8, 5);
		ok&=lengths("one-base-utrs", gene(Shared.PLUS, new int[][]{{0, 0}, {10, 15}, {20, 20}},
				10, 15, "seqGene", false), 8, 6, 1, 1);
		require(ok, "Each whole UTR exon must contribute all of its inclusive bases");
	}

	/** Null and unrecognized receiver sources remain absent even when the donor has a source. */
	private static void missingSource(){
		boolean ok=true;
		for(String source : new String[]{null, "unlisted-source"}){
			final Gene receiver=gene(Shared.PLUS, new int[][]{{0, 27}}, 10, 15, source, false);
			final Gene donor=gene(Shared.PLUS, new int[][]{{0, 27}}, 10, 15, "seqGene", false);
			Gene merged=null;
			RuntimeException failure=null;
			try{merged=receiver.merge(donor);}
			catch(RuntimeException e){failure=e;}
			final boolean valid=failure==null && merged.primarySource==-1 && sourceColumn(merged).isEmpty() &&
					merged.exons==receiver.exons && receiver.primarySource==-1;
			System.out.println("source="+source+" failure="+(failure==null ? "none" : failure.getClass().getSimpleName())+
					" preserved="+valid);
			ok&=valid;
		}
		require(ok, "Missing source indices must not be used to index sourceCodes during merge");
	}

	/** Builds an annotation from explicit inclusive spans; expected counts are separate literals. */
	private static Gene gene(byte strand, int[][] spans, int start, int stop, String source, boolean untranslated){
		final Exon[] exons=new Exon[spans.length];
		for(int i=0; i<spans.length; i++){
			final int a=spans[i][0], b=spans[i][1];
			exons[i]=new Exon(a, b, (byte)1, strand, a<start || b>stop, start<stop && a<=stop && b>=start);
		}
		return new Gene((byte)1, strand, spans[0][0], spans[spans.length-1][1], start, stop, 7,
				"fixture", "transcript", "protein", null, null, exons, untranslated, false, true,
				source, "description", "full description");
	}

	/** Reports observed and independent expected lengths before aggregating baseline discrepancies. */
	private static boolean lengths(String label, Gene gene, int exon, int coding, int five, int three){
		final boolean ok=gene.exonLength==exon && gene.exonCodeLength==coding &&
				gene.utrLength5prime==five && gene.utrLength3prime==three;
		System.out.println(label+" exon="+gene.exonLength+" coding="+gene.exonCodeLength+
				" five="+gene.utrLength5prime+" three="+gene.utrLength3prime+" pass="+ok);
		return ok;
	}

	/** The documented legacy annotation format places source in zero-based column 14. */
	private static String sourceColumn(Gene gene){
		final String[] fields=gene.toString().split("\t", -1);
		require(fields.length==18, "Fixture descriptions must retain all 18 annotation columns");
		return fields[14];
	}

	/** Keeps independent oracle checks active when Java assertions are disabled. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Five left UTR bases, six coding bases and eight right UTR bases, separated by introns. */
	private static final int[][] WHOLE={{0, 4}, {10, 15}, {20, 27}};
}
