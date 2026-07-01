package var2;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.HashSet;

import fileIO.ByteStreamWriter;
import fileIO.TextFile;
import parse.Parse;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Converts CallVariants VCF output into training vectors for neural network training.
 * Loads truth VCF as simple (chrom,pos,ref,alt) keys — compatible with any VCF format
 * including GIAB. Loads candidates with VcfLoader for full Var objects.
 *
 * Usage: java var2.VcfToTrainingVectors in=candidates.vcf truth=truth.vcf
 *        ref=reference.fa out=vectors.tsv ploidy=2 platform=illumina
 *
 * @author UMP45
 */
public class VcfToTrainingVectors {

	public static void main(String[] args){
		String inFile=null, truthFile=null, refFile=null, outFile=null;
		int ploidy=2;

		for(String arg : args){
			String[] split=arg.split("=", 2);
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("in")){inFile=b;}
			else if(a.equals("truth") || a.equals("ref2")){truthFile=b;}
			else if(a.equals("ref")){refFile=b;}
			else if(a.equals("out")){outFile=b;}
			else if(a.equals("ploidy")){ploidy=Integer.parseInt(b);}
			else if(a.equals("includescore")){
				VectorUMP45.includeScore=Parse.parseBoolean(b);
			}else if(a.equals("platform")){
				if(b.equalsIgnoreCase("illumina")){VectorUMP45.platform=VectorUMP45.PLATFORM_ILLUMINA;}
				else if(b.equalsIgnoreCase("pacbio")){VectorUMP45.platform=VectorUMP45.PLATFORM_PACBIO;}
				else if(b.equalsIgnoreCase("nanopore") || b.equalsIgnoreCase("ont")){VectorUMP45.platform=VectorUMP45.PLATFORM_NANOPORE;}
				else if(b.equalsIgnoreCase("roche") || b.equalsIgnoreCase("sbx")){VectorUMP45.platform=VectorUMP45.PLATFORM_ROCHE;}
				else{VectorUMP45.platform=Integer.parseInt(b);}
			}
		}

		assert(inFile!=null) : "Missing in= parameter";
		assert(truthFile!=null) : "Missing truth= parameter";
		assert(outFile!=null) : "Missing out= parameter";

		PrintStream out=System.err;

		// Load reference into ScafMap (needed for homopolymer and contig end distance)
		ScafMap scafMap;
		if(refFile!=null){
			out.println("Loading reference: "+refFile);
			scafMap=ScafMap.loadReference(refFile, true);
		}else{
			out.println("Loading scaffold map from VCF header: "+inFile);
			scafMap=ScafMap.loadVcfHeader(inFile);
		}

		// Load truth VCF as simple variant keys (works with ANY VCF format)
		out.println("Loading truth: "+truthFile);
		HashSet<String> truthSet=loadTruthKeys(truthFile, scafMap);
		out.println("Truth variant keys loaded: "+truthSet.size());

		// Load candidate VCF with full extended info (BBTools format)
		out.println("Loading candidates: "+inFile);
		VarMap candidateMap=VcfLoader.loadVcfFile(inFile, scafMap, true, true);
		long candidateSize=candidateMap.size();
		out.println("Candidate variants loaded: "+candidateSize);

		double pairingRate=candidateMap.properPairRate;
		double totalQualityAvg=candidateMap.totalQualityAvg;
		double totalMapqAvg=candidateMap.totalMapqAvg;
		double readLengthAvg=candidateMap.readLengthAvg;

		out.println("Dataset stats: pairingRate="+pairingRate+
				" qualAvg="+totalQualityAvg+
				" mapqAvg="+totalMapqAvg+
				" readLenAvg="+readLengthAvg);

		// Write training vectors
		ByteStreamWriter bsw=new ByteStreamWriter(outFile, true, false, false);
		bsw.start();

		ByteBuilder bb=new ByteBuilder();
		bb.append("#dims\t").append(VectorUMP45.DIMS).append('\t').append(1).nl();
		bsw.print(bb);

		long written=0, excluded=0, positive=0, negative=0;

		for(Var v : candidateMap){
			int count=v.alleleCount();
			if(count<1){continue;}

			v.calcCoverage(scafMap);

			// Compute phred score for labeling
			double phred=v.phredScore(pairingRate, totalQualityAvg, totalMapqAvg,
					readLengthAvg, 1, ploidy, scafMap, null);

			// Check truth membership using key
			String key=makeKey(v, scafMap);
			//[var2/VcfToTrainingVectors#001 FIXED] DELs now match: makeKeyFromVcfFields was made anchor-inclusive
			//to agree with makeKey/fromVCF. See the FIXED note in makeKeyFromVcfFields for the full rationale.
			boolean inTruth=truthSet.contains(key);

			// Assign continuous label
			// Deadzone (25-30) only applies to NIST-NEGATIVE variants
			double label;
			if(inTruth){
				label=0.7+Tools.mid(0, phred, 30)/100.0;
			}else if(phred>=30){
				label=0.5+Tools.mid(30, phred, 50)/100.0;
			}else if(phred>25){
				excluded++;
				continue; // Ambiguous zone — NIST-negative only
			}else{
				label=phred/100.0;
			}
			label=Tools.mid(0, label, 1);

			// Compute feature vector
			float[] vec=VectorUMP45.makeVector(v, pairingRate, totalQualityAvg,
					totalMapqAvg, readLengthAvg, ploidy, scafMap);

			// Write TSV line
			bb.clear();
			for(int i=0; i<vec.length; i++){
				if(i>0){bb.tab();}
				bb.append(vec[i], 6);
			}
			bb.tab().append(label, 6).nl();
			bsw.print(bb);

			written++;
			if(label>0.5){positive++;}
			else{negative++;}
		}

		bsw.poisonAndWait();

		out.println("Written: "+written+" (positive: "+positive+", negative: "+negative+")");
		out.println("Excluded (ambiguous 25-30): "+excluded);
		out.println("Output: "+outFile);
	}

	/**
	 * Loads truth VCF as a set of canonical variant keys.
	 * Parses only CHROM, POS, REF, ALT — works with any VCF format (GIAB, GATK, etc.).
	 */
	private static HashSet<String> loadTruthKeys(String fname, ScafMap scafMap){
		HashSet<String> set=new HashSet<String>();
		TextFile tf=new TextFile(fname);
		for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
			if(line.startsWith("#")){continue;}
			String[] fields=line.split("\t", 6);
			if(fields.length<5){continue;}

			String chrom=fields[0];
			int pos=Integer.parseInt(fields[1]);
			String ref=fields[3];
			String alt=fields[4];

			// Handle multi-allelic (comma-separated ALT)
			String[] alts=alt.split(",");
			for(String a : alts){
				String key=makeKeyFromVcfFields(chrom, pos, ref, a.trim(), scafMap);
				if(key!=null){set.add(key);}
			}
		}
		tf.close();
		return set;
	}

	/**
	 * Creates canonical variant key from VCF fields for matching.
	 * Converts VCF coordinates to internal 0-based format.
	 */
	private static String makeKeyFromVcfFields(String chrom, int pos, String ref, String alt, ScafMap scafMap){
		int scafNum=scafMap.getNumber(chrom);
		if(scafNum<0){return null;}

		int refLen=ref.length();
		int altLen=alt.length();

		// Determine type and normalize coordinates (same logic as VcfToVar)
		int start, stop;
		String normalizedAlt;

		if(refLen!=altLen && altLen>0 && refLen>0){
			// Indel: strip leading base
			normalizedAlt=alt.substring(1);
			start=pos; // stays 1-based then converted
			//FIXED [var2/VcfToTrainingVectors#001] (was HIGH; G11 2026-06-29, UMP45-cleared): mirror
			//VcfToVar.fromVCF's DEL convention. The old `refLen--` (unconditional) made DEL truth keys
			//anchor-EXCLUSIVE (stop=pos+len), but candidate keys (from fromVCF) are anchor-INCLUSIVE
			//(stop=pos+len+1), so truthSet.contains() was ALWAYS false for deletions -> true DELs mislabeled
			//as non-truth -> wrong NN training targets. Now `if(refLen==1){refLen--;}` decrements only for the
			//insertion case (matching fromVCF), so DEL keys line up. SUB/INS unchanged; degenerate guard below
			//unaffected. (No shipped net used this tool, per UMP45 — low-risk fix, no retrain needed.)
			if(refLen==1){refLen--;}
			if(refLen==0 && normalizedAlt.isEmpty()){
				// Pure insertion of 1 base that got stripped — degenerate
				return null;
			}
		}else{
			normalizedAlt=alt;
			start=pos-1; // convert to 0-based
		}
		stop=start+refLen;

		return scafNum+":"+start+":"+stop+":"+normalizedAlt.toUpperCase();
	}

	/** Creates matching key from internal Var object */
	private static String makeKey(Var v, ScafMap scafMap){
		return v.scafnum+":"+v.start+":"+v.stop+":"+new String(v.allele).toUpperCase();
	}
}
