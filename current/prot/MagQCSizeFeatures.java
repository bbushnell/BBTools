package prot;

import fileIO.ByteFile;
import parse.LineParser1;
import structures.ByteBuilder;

/** Frozen universal marker observations and target transforms for the v1.4 size subnet.
 * Complete-source labels are accepted only by targets(); input encoding sees observed
 * family counts alone. The linear size head remains unbounded and uses identity activation.
 * @author Brian Bushnell, Yoimiya
 */
final class MagQCSizeFeatures {

	/** Loads the ordered marker list only after checking its resource and family-order bindings. */
	static MagQCSizeFeatures load(final String path, final String pin, final String requirednessPin,
			final String familyPin, final int familyCount){
		require(familyCount>=MARKERS, "Size inputs require a vocabulary with at least 100 families");
		require(DigestSuffix.file(path).equals(DigestSuffix.requireSuffix(pin, "size markers")), "Size marker bytes differ from their pin");
		DigestSuffix.requireSuffix(requirednessPin, "requiredness"); DigestSuffix.requireSuffix(familyPin, "families");
		final ByteFile input=ByteFile.makeByteFile(path, false);
		final LineParser1 parser=new LineParser1('\t');
		final int[] ranks=new int[MARKERS]; final boolean[] seen=new boolean[familyCount];
		try{
			header(input, parser, "#schema_version", "magqc_size_markers_v1");
			header(input, parser, "#requiredness_sha80", requirednessPin);
			header(input, parser, "#family_sha80", familyPin);
			header(input, parser, "#selection", "single_copy_prevalence_descending_then_family_rank_ascending");
			final byte[] columns=input.nextLine();
			require(columns!=null && java.util.Arrays.equals(columns, COLUMNS), "Size marker columns differ");
			long previous=Long.MAX_VALUE, genomes=-1; int previousRank=-1;
			for(int order=0; order<MARKERS; order++){
				final byte[] row=input.nextLine(); require(row!=null && row.length>0, "Incomplete size marker list");
				parser.set(row); require(parser.terms()==4 && number(parser, 0)==order, "Size marker ordering differs");
				final long rank=number(parser, 1), singletons=number(parser, 2), denominator=number(parser, 3);
				require(rank<familyCount && !seen[(int)rank] && denominator>0 && singletons<=denominator,
					"Invalid or repeated size marker rank/count");
				require(singletons<previous || singletons==previous && rank>previousRank, "Size markers are not in deterministic prevalence order");
				if(order==0){genomes=denominator;}
				require(denominator==genomes, "Size markers must use one complete-corpus denominator");
				ranks[order]=(int)rank; seen[(int)rank]=true; previous=singletons; previousRank=(int)rank;
			}
			require(input.nextLine()==null, "Extra size marker records");
		}finally{require(!input.close(), "Size marker read failed");}
		return new MagQCSizeFeatures(ranks, familyCount);
	}

	/** Keeps marker ordering private so all workers can safely share one frozen resource. */
	private MagQCSizeFeatures(final int[] ranks_, final int families_){ranks=ranks_.clone(); families=families_;}

	/** Appends 100 observed presence bits and returns their mean for the shared numeric formatter. */
	double appendPresence(final ByteBuilder out, final int[] observed){
		require(observed!=null && observed.length==families, "Size observations differ from the frozen family width");
		int present=0;
		for(int rank:ranks){
			final int count=observed[rank]; require(count>=0, "Negative observed size marker count");
			if(count>0){present++;}
			out.append(count>0 ? '1' : '0').tab();
		}
		return present/100.0;
	}

	/** Produces five value labels; returns bit0/bit1 when the BP/CDS log target was clamped.
	 * The linear size value remains unbounded. Dynamic residual targets belong to the trainer. */
	static int targets(final long originalBp, final long originalCds, final double geneC,
			final double geneT, final double[] out){
		require(originalBp>0 && originalCds>0 && unit(geneC) && unit(geneT) && out!=null && out.length==VALUES,
			"Size targets require positive original counts and bounded CDS quality labels");
		final double bp=(Math.log10(originalBp)-BP_LOG_OFFSET)/BP_LOG_SCALE;
		final double cds=(Math.log10(originalCds)-CDS_LOG_OFFSET)/CDS_LOG_SCALE;
		out[0]=Math.max(0,Math.min(1,bp));
		out[1]=originalBp/BP_LINEAR_SCALE;
		out[2]=Math.max(0,Math.min(1,cds));
		out[3]=geneC; out[4]=geneT;
		for(double value:out){assert(Double.isFinite(value)) : "Positive integer source sizes must produce finite affine log and linear targets";}
		return (bp<0 || bp>1 ? 1 : 0) | (cds<0 || cds>1 ? 2 : 0);
	}

	/** Uses the linear size head with a one-bp lower floor, without exponentiation or an upper biological cutoff. */
	static double logRatio(final long binBp, final double sizeLinear){
		require(binBp>0 && Double.isFinite(sizeLinear), "Size ratio requires positive observed bp and a finite linear prediction");
		// Avoid overflow even for a finite double supplied by an external consumer.
		final double predictedLog=sizeLinear<=1/BP_LINEAR_SCALE ? 0 : Math.log(sizeLinear)+Math.log(BP_LINEAR_SCALE);
		return (Math.log(binBp)-predictedLog)/Math.log(2);
	}

	/** Reads exact ordered metadata; malformed resources fail with assertions disabled as well. */
	private static void header(final ByteFile input, final LineParser1 parser, final String key, final String value){
		final byte[] row=input.nextLine(); require(row!=null, "Missing size marker metadata: "+key);
		parser.set(row); require(parser.terms()==2 && parser.termEquals(key, 0) && parser.termEquals(value, 1), "Size marker binding differs: "+key);
	}

	/** Unsigned decimal parsing rejects overflow and alternate numeric spellings. */
	private static long number(final LineParser1 parser, final int field){
		long value=0; require(parser.length(field)>0, "Missing marker integer");
		for(int i=0; i<parser.length(field); i++){
			final int digit=parser.parseByte(field, i)-'0';
			require(digit>=0 && digit<=9 && value<=(Long.MAX_VALUE-digit)/10, "Invalid marker integer"); value=value*10+digit;
		}
		return value;
	}
	/** CDS quality labels are fractions; the linear size target remains unbounded. */
	private static boolean unit(final double value){return Double.isFinite(value) && value>=0 && value<=1;}
	/** Resource and label errors must fail equally in production and fixture modes. */
	private static void require(final boolean condition, final String message){if(!condition){throw new IllegalArgumentException(message);}}

	private final int[] ranks;
	private final int families;
	static final int MARKERS=100, VALUES=5, OUTPUTS=10;
	static final double BP_LOG_OFFSET=5, BP_LOG_SCALE=3, CDS_LOG_OFFSET=1, CDS_LOG_SCALE=4, BP_LINEAR_SCALE=10000000;
	static final String OUTPUT_NAMES="size01,size_lin,cds01,gene_completeness,gene_contamination,size01_error,size_lin_error,cds01_error,gene_completeness_error,gene_contamination_error";
	private static final byte[] COLUMNS="#columns\torder\tfamily_rank\tsingle_copy_genomes\tcomplete_genomes".getBytes(java.nio.charset.StandardCharsets.US_ASCII);
}
